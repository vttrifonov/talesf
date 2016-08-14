/*

Copyright (c) 2011-2012, Daniel S. Standage <daniel.standage@gmail.com> and
Erin Doyle <edoyle@iastate.edu> with modifications by Nick Booher <njbooher@gmail.com>.
See README for license details.

*/

// System libraries
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <bcutils/Hashmap.h>
#include <bcutils/Array.h>
#include <bcutils/bcutils.h>

// Initialize the kseq library
#include <bcutils/kseq.h>
KSEQ_INIT(gzFile, gzread)

typedef struct
{
  int strand;
  char *sequence;
  char *sequence_name;
  unsigned long index;
  double score;
} BindingSite;

Hashmap *talesf_kwargs;

/*
 * Utility
 */

void create_options_string(char *options_str, char *rvd_str) {

  char cutoff_str[32];
  char rvds_eq_str[256];

  double cutoff = *((double *) hashmap_get(talesf_kwargs, "cutoff"));
  int forward_only = *((int *) hashmap_get(talesf_kwargs, "forward_only"));
  int c_upstream = *((int *) hashmap_get(talesf_kwargs, "c_upstream"));

  strcat(options_str, "options_used:");

  if (!forward_only) {
    strcat(options_str, "search reverse complement, ");
  }

  strcat(options_str, "upstream_base = ");

  if (c_upstream != 1) {
    strcat(options_str, "T ");
  }

  if (c_upstream != 0) {
    strcat(options_str, "C ");
  }

  sprintf(cutoff_str, ", cutoff = %.2lf, ", cutoff);
  strcat(options_str, cutoff_str);

  sprintf(rvds_eq_str, "rvd_sequence = %s", rvd_str);
  strcat(options_str, rvds_eq_str);

  strcat(options_str, "\n");

}

double *create_lookahead_array(unsigned int *rvd_seq, unsigned int rvd_seq_len, double cutoff, double best_score, double **scoring_matrix) {
  
  double *lookahead_array = calloc(rvd_seq_len, sizeof(double));
  
  lookahead_array[rvd_seq_len - 1] = cutoff * best_score;
  
  for (int i = rvd_seq_len - 2; i >= 0; i--) {
    
    double *scores = scoring_matrix[rvd_seq[i+1]];
    
    double min = scores[0];
    
    for (int j = 0; j < 4; j++) {
      if (scores[j] < min) {
        min = scores[j];
      }
    }
    
    lookahead_array[i] = lookahead_array[i + 1] - (min);
      
  }
  
  return lookahead_array;

}

/*
 * Core
 */

int binding_site_compare_score(const void * a, const void * b)
{
  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  double real_a_score = floorf(real_a->score * 100 + 0.5) / 100;
  double real_b_score = floorf(real_b->score * 100 + 0.5) / 100;

  double score_diff = (real_a_score - real_b_score);

  if(score_diff < 0) {
    return -1;
  } else if(score_diff > 0) {
    return 1;
  } else {
    return 0;
  }

}

int print_results(Array *results, FILE *log_file) {
  
  char *output_filepath = hashmap_get(talesf_kwargs, "output_filepath");
  char *organism_name = hashmap_get(talesf_kwargs, "organism_name");
  unsigned int num_rvds = *((unsigned int *) hashmap_get(talesf_kwargs, "rvd_seq_len"));
  double best_score = *((double *) hashmap_get(talesf_kwargs, "best_score"));
  int forward_only = *((int *) hashmap_get(talesf_kwargs, "forward_only"));
  
  char *source_str = "TALESF";
  char options_str[512];
  
  // strcat doesn't seem to work unless you do this
  *options_str = '\0';

  char *plus_strand_sequence;
  FILE *gff_out_file = NULL;
  FILE *tab_out_file = NULL;
  FILE *genome_browser_file = NULL;

  int is_genome = (*organism_name != '\0');
  int genome_using_gbrowse = (is_genome && (strcmp(organism_name, "oryza_sativa") == 0 || strcmp(organism_name, "arabidopsis_thaliana") == 0));
  
  size_t output_filepath_length;
  char* temp_output_filepath;
  
  char *rvd_string_printable = strdup(hashmap_get(talesf_kwargs, "rvd_string"));
  
  char *pos = strstr(rvd_string_printable, " ");
  
  while (pos != NULL) {

    strncpy(pos, "_", 1);
    pos = strstr(rvd_string_printable, " ");

  }

  create_options_string(options_str, rvd_string_printable);

  output_filepath_length = strlen(output_filepath) + 5;
  temp_output_filepath = calloc(output_filepath_length + 1, sizeof(char));
  
  sprintf(temp_output_filepath, "%s.txt", output_filepath);
  tab_out_file = fopen(temp_output_filepath, "w");
  memset(temp_output_filepath, '\0', output_filepath_length);
  sprintf(temp_output_filepath, "%s.gff3", output_filepath);
  gff_out_file = fopen(temp_output_filepath, "w");
  
  if(is_genome) {
    
    memset(temp_output_filepath, '\0', output_filepath_length);
    
    if(genome_using_gbrowse) {
      sprintf(temp_output_filepath, "%s.gff", output_filepath);
    } else {
      sprintf(temp_output_filepath, "%s.bed", output_filepath);
    }
    
    genome_browser_file = fopen(temp_output_filepath, "w");
    
    
  }
  
  free(temp_output_filepath);

  if(!gff_out_file || !tab_out_file || (is_genome && !genome_browser_file)) {
    
    fprintf(log_file, "Error: unable to open output file '%s'\n", output_filepath);
    
    free(rvd_string_printable);
    
    return 1;
    
  }

  // Tab file header
  if (forward_only) {
    fprintf(tab_out_file, "table_ignores:Plus strand sequence\n");
  }

  fprintf(tab_out_file, options_str);

  fprintf(tab_out_file, "Best Possible Score:%.2lf\n", best_score);
  fprintf(tab_out_file, "Sequence Name\tStrand\tScore\tStart Position\tTarget Sequence\tPlus strand sequence\n");

  // GFF file header
  fprintf(gff_out_file, "##gff-version 3\n");

  if (forward_only) {
    fprintf(gff_out_file, "#table_display_tags:target_sequence\n");
  } else {
    fprintf(gff_out_file, "#table_display_tags:target_sequence,plus_strand_sequence\n");
  }

  fprintf(gff_out_file, "#%s", options_str);

  fprintf(gff_out_file, "#Best Possible Score:%.2lf\n", best_score);

  // bed file header

  if(genome_using_gbrowse) {
    fprintf(genome_browser_file, "##gff-version 3\n");
  }
  else if(is_genome) {
    fprintf(genome_browser_file, "track name=\"TAL Targets\" description=\"Targets for RVD sequence %s\" visibility=2 useScore=1\n", rvd_string_printable);
  }

  for(int i = 0; i < array_size(results); i++) {

    BindingSite *site = (BindingSite *)array_get(results, i);
    char *sequence = site->sequence;
    char strand = '+';
    char *tab_strand = "Plus";

    if(site->strand > 0)
      plus_strand_sequence = sequence;
    else {
      
      int seq_len = num_rvds + 2;

      plus_strand_sequence = sequence;
      sequence = malloc(sizeof(char)*(seq_len+1));
      sequence[seq_len] = '\0';

      for(int j = 0; j < seq_len; j++) {
        
        char base = site->sequence[seq_len - j - 1];
        
        if(base == 'A' || base == 'a')
          sequence[j] = 'T';
        else if(base == 'C' || base == 'c')
          sequence[j] = 'G';
        else if(base == 'G' || base == 'g')
          sequence[j] = 'C';
        else if(base == 'T' || base == 't')
          sequence[j] = 'A';
        else if(base == ' ')
          sequence[j] = ' ';
        else {
          fprintf(stderr, "Error: unexpected character '%c'\n", base);
          free(sequence);
          free(rvd_string_printable);
          fclose(gff_out_file);
          fclose(tab_out_file);
          if (genome_browser_file) {
            fclose(genome_browser_file);
          }
          return 1;
        }
        
      }
      
      strand = '-';
      tab_strand = "Minus";
      
    }

    fprintf( tab_out_file, "%s\t%s\t%.2lf\t%lu\t%s\t%s\n",
             site->sequence_name, tab_strand, site->score, site->index + 1, sequence, plus_strand_sequence);

    fprintf( gff_out_file, "%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\trvd_sequence=%s;target_sequence=%s;plus_strand_sequence=%s;\n",
             site->sequence_name, source_str, "TAL_effector_binding_site", site->index + 1,
             site->index + num_rvds, site->score, strand, rvd_string_printable, sequence, plus_strand_sequence);

    if(is_genome && i < 10000) {

      if(genome_using_gbrowse) {

        fprintf( genome_browser_file, "chr%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\tName=site%d;\n",
                 site->sequence_name, source_str, "TAL_effector_binding_site", site->index + 1,
                 site->index + num_rvds, site->score, strand, i);

      } else {

        int bed_score = floorf((best_score / site->score * 1000) + 0.5);
        fprintf( genome_browser_file,"%s\t%lu\t%lu\tsite%d\t%d\t%c\n",
                 site->sequence_name, site->index, site->index + num_rvds - 1, i, bed_score, strand);

      }

    }

    if(plus_strand_sequence != sequence) {
      free(sequence);
    }

  }

  free(rvd_string_printable);
  fclose(gff_out_file);
  fclose(tab_out_file);

  if(is_genome) {
    fclose(genome_browser_file);
  }

  return 0;

}

double score_binding_site(kseq_t *seq, unsigned long i, unsigned int *rvd_seq, unsigned int rvd_seq_len, double **scoring_matrix, double *lookahead_array, int reverse) {

  double total_score = 0.0;
  int num_rvds = rvd_seq_len;

  if (!reverse) {

    for (unsigned long j = 0; j < rvd_seq_len; j++) {

      double *scores = scoring_matrix[rvd_seq[j]];

      if (seq->seq.s[i+j] == 'A' || seq->seq.s[i+j] == 'a')
        total_score += scores[0];
      else if (seq->seq.s[i+j] == 'C' || seq->seq.s[i+j] == 'c')
        total_score += scores[1];
      else if (seq->seq.s[i+j] == 'G' || seq->seq.s[i+j] == 'g')
        total_score += scores[2];
      else if (seq->seq.s[i+j] == 'T' || seq->seq.s[i+j] == 't')
        total_score += scores[3];
      else
        total_score += lookahead_array[num_rvds - 1] + 1;

      if (total_score > lookahead_array[j])
        return -1;

    }

  } else {

    for (unsigned long j = 0; j < rvd_seq_len; j++) {

      double *scores = scoring_matrix[rvd_seq[j]];

      unsigned long k = i + rvd_seq_len - j - 2;

      if (seq->seq.s[k] == 'A' || seq->seq.s[k] == 'a')
        total_score += scores[3];
      else if (seq->seq.s[k] == 'C' || seq->seq.s[k] == 'c')
        total_score += scores[2];
      else if (seq->seq.s[k] == 'G' || seq->seq.s[k] == 'g')
        total_score += scores[1];
      else if (seq->seq.s[k] == 'T' || seq->seq.s[k] == 't')
        total_score += scores[0];
      else
        total_score += lookahead_array[num_rvds - 1] + 1;

      if (total_score > lookahead_array[j])
        return -1;
    }

  }

  return total_score;

}

BindingSite *create_binding_site(kseq_t *seq, unsigned long i, int num_rvds, double score, int reverse) {

  int seq_name_len = strlen(seq->name.s);
  
  BindingSite *site = malloc(sizeof(BindingSite));
  
  site->sequence = calloc(num_rvds + 2 + 1, sizeof(char));
  site->sequence[num_rvds + 2] = '\0';
  
  site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
  site->sequence_name[seq_name_len] = '\0';
  strncpy(site->sequence_name, seq->name.s, seq_name_len);
  
  site->score = score;
  
  if (!(reverse == 1)) {
    
    site->strand = 1;
    site->index = i;
    
    strncpy(site->sequence, seq->seq.s + site->index - 1, 1);
    site->sequence[1] = ' ';
    strncpy(site->sequence + 2, seq->seq.s + site->index, num_rvds);
    
  } else {
    
    site->strand = -1;
    site->index = i - 1;
    
    strncpy(site->sequence, seq->seq.s + site->index, num_rvds);
    site->sequence[num_rvds] = ' ';
    strncpy(site->sequence + num_rvds + 1, seq->seq.s + site->index + num_rvds, 1);
    
  }
  
  for(int j = 0; j < num_rvds + 2 + 1; j++) {
    site->sequence[j] = toupper(site->sequence[j]);
  }
  
  return site;

}

void cpu_the_whole_shebang(kseq_t *seq, double *lookahead_array, Array *results) {
  
  int c_upstream = *((int *) hashmap_get(talesf_kwargs, "c_upstream"));
  int forward_only = *((int *) hashmap_get(talesf_kwargs, "forward_only"));
  unsigned int *rvd_seq = hashmap_get(talesf_kwargs, "rvd_seq");
  unsigned int num_rvds = *((unsigned int *) hashmap_get(talesf_kwargs, "rvd_seq_len"));
  double **scoring_matrix = hashmap_get(talesf_kwargs, "scoring_matrix");
  
  #pragma omp parallel for schedule(static)
  for(unsigned long i = 1; i <= seq->seq.l - num_rvds; i++) {
    
    if((c_upstream != 0 && (seq->seq.s[i-1] == 'C' || seq->seq.s[i-1] == 'c')) || (c_upstream != 1 && (seq->seq.s[i-1] == 'T' || seq->seq.s[i-1] == 't'))) {
      
      double score = score_binding_site(seq, i, rvd_seq, num_rvds, scoring_matrix, lookahead_array, 0);
      
      if(score != -1) {
        
        BindingSite *site = create_binding_site(seq, i, num_rvds, score, 0);
        
        #pragma omp critical (add_result)
        array_add(results, site);
        
      }
      
    }
    
    if(!forward_only) {
      
      if((c_upstream != 0 && (seq->seq.s[i + num_rvds - 1] == 'G' || seq->seq.s[i + num_rvds - 1] == 'g')) || (c_upstream != 1 && (seq->seq.s[i + num_rvds - 1] == 'A' || seq->seq.s[i + num_rvds - 1] == 'a'))) {
        
        double score = score_binding_site(seq, i, rvd_seq, num_rvds, scoring_matrix, lookahead_array, 1);
        
        if(score != -1) {
          
          BindingSite *site = create_binding_site(seq, i, num_rvds, score, 1);
          
          #pragma omp critical (add_result)
          array_add(results, site);
          
        }
        
      }
    }
    
  }
  
}

// Identify and print out TAL effector binding sites
void find_binding_sites(FILE *log_file, kseq_t *seq, double *lookahead_array, Array *results) {
  
  unsigned int num_rvds = *((unsigned int *) hashmap_get(talesf_kwargs, "rvd_seq_len"));
  
  if(num_rvds > seq->seq.l) {
    logger(log_file, "Warning: skipping sequence '%s' since it is shorter than the RVD sequence\n", seq->seq.s);
    return;
  }
  
  logger(log_file, "Scanning %s for binding sites (length %ld)", seq->name.s, seq->seq.l);

  cpu_the_whole_shebang(seq, lookahead_array, results);

}

int run_talesf_task(Hashmap *kwargs) {
  
  talesf_kwargs = kwargs;
  
  // Options
  char *seq_filename = hashmap_get(kwargs, "seq_filename");
  char *log_filepath = hashmap_get(kwargs, "log_filepath");
  unsigned int *rvd_seq = hashmap_get(kwargs, "rvd_seq");
  unsigned int rvd_seq_len = *((unsigned int *) hashmap_get(kwargs, "rvd_seq_len"));
  double best_score = *((double *) hashmap_get(kwargs, "best_score"));
  double cutoff = *((double *) hashmap_get(kwargs, "cutoff"));
  int numprocs = *((int *) hashmap_get(kwargs, "num_procs"));
  double **scoring_matrix = hashmap_get(kwargs, "scoring_matrix");
  
  // Setup the logger

  FILE *log_file = stdout;

  if (log_filepath && strcmp(log_filepath, "NA") != 0) {
    log_file = fopen(log_filepath, "a");
  }
  
  // Open sequence file
  
  gzFile seqfile;
  
  seqfile = gzopen(seq_filename, "r");
  
  if (!seqfile) {
    logger(log_file, "Error: unable to open sequence '%s'", seq_filename);
    if (log_file != stdout) {
      fclose(log_file);
    }
    return 1;
  }
  
  Array *results = array_new( sizeof(BindingSite *) );
  
  // Define score cutoffs for match sites
  
  double *lookahead_array = create_lookahead_array(rvd_seq, rvd_seq_len, cutoff, best_score, scoring_matrix);
  
  // Begin processing

  int abort = 0;
  
  omp_set_num_threads(numprocs);
  
  kseq_t *seq = kseq_init(seqfile);
  
  int result;
  
  while ((result = kseq_read(seq)) >= 0) {
    find_binding_sites(log_file, seq, lookahead_array, results);
  }
  
  kseq_destroy(seq);
  gzclose(seqfile);

  if(!abort) {

    qsort(results->data, array_size(results), sizeof(BindingSite *), binding_site_compare_score);

    abort = print_results(results, log_file);

    logger(log_file, "Finished");

  }

  // Free memory

  if(results) {

    for(int i = 0; i < array_size(results); i++) {
      
      BindingSite *site = (BindingSite *)array_get(results, i);

      free(site->sequence);
      free(site->sequence_name);
      free(site);

    }

    array_delete(results, NULL);

  }
  
  if (lookahead_array) {
    free(lookahead_array);
  }

  if(log_file != stdout) {
    fclose(log_file);
  }

  return abort;
  
}
