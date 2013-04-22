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

typedef struct {
  char *sequence[2];
  char *sequence_name;
  unsigned long indexes[2];
  double scores[2];
  int spacer_length;
  int f_idx;
  int r_idx;
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

double *create_lookahead_array(Array *rvd_seq, double cutoff, double best_score, Hashmap *diresidue_scores) {
  
  double *lookahead_array = calloc(array_size(rvd_seq), sizeof(double));
  
  lookahead_array[array_size(rvd_seq) - 1] = cutoff * best_score;
  
  for (int i = array_size(rvd_seq) - 2; i >= 0; i--) {
    
    double *scores = hashmap_get(diresidue_scores, array_get(rvd_seq, i + 1));
    
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

int binding_site_compare_score(const void * a, const void * b) {
  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  int str_cmp_result = strcmp(real_a->sequence_name, real_b->sequence_name);

  if (str_cmp_result != 0) {

    return str_cmp_result;

  } else {

    double real_a_score = floorf(real_a->scores[0] * 100 + 0.5) / 100;
    double real_b_score = floorf(real_b->scores[0] * 100 + 0.5) / 100;

    double score_diff = (real_a_score - real_b_score);

    if (score_diff < 0) {
      return -1;
    } else if (score_diff > 0) {
      return 1;
    } else {

      double real_a_score2 = floorf(real_a->scores[1] * 100 + 0.5) / 100;
      double real_b_score2 = floorf(real_b->scores[1] * 100 + 0.5) / 100;
      double score_diff2 = (real_a_score2 - real_b_score2);

      if (score_diff2 < 0) {
        return -1;
      } else if (score_diff2 > 0) {
        return 1;
      } else {
        return 0;
      }

    }

  }

}

int binding_site_compare_pos(const void * a, const void * b) {

  BindingSite *real_a = *((BindingSite **)a);
  BindingSite *real_b = *((BindingSite **)b);

  int str_cmp_result = strcmp(real_a->sequence_name, real_b->sequence_name);

  if (str_cmp_result != 0) {

    return str_cmp_result;

  } else if (real_a->f_idx != real_b->f_idx) {

    return ((real_a->f_idx < real_b->f_idx) ? -1 : 1);

  } else if (real_a->r_idx != real_b->r_idx) {

    return ((real_a->r_idx < real_b->r_idx) ? -1 : 1);

  } else if (real_a->indexes[0] != real_b->indexes[0]) {

    return ((real_a->indexes[0] < real_b->indexes[0]) ? -1 : 1);

  } else if (real_a->indexes[1] != real_b->indexes[1]) {

    return ((real_a->indexes[1] < real_b->indexes[1]) ? -1 : 1);

  } else {

    return 0;

  }

}

int print_results(Array *results, Array **rvd_seqs, double best_score, double best_score2, FILE *log_file) {

  char *output_filepath = hashmap_get(talesf_kwargs, "output_filepath");
  char *organism_name = hashmap_get(talesf_kwargs, "organism_name");

  double combined_best_score = best_score + best_score2;

  char *source_str = "TALESF";
  char options_str[512];

  // strcat doesn't seem to work unless you do this
  *options_str = '\0';

  FILE *tab_out_file = NULL;
  FILE *gff_out_file = NULL;
  FILE *genome_browser_file = NULL;

  int is_genome = (*organism_name != '\0');
  int genome_using_gbrowse = (is_genome && (strcmp(organism_name, "oryza_sativa") == 0 || strcmp(organism_name, "arabidopsis_thaliana") == 0));

  size_t output_filepath_length;
  char* temp_output_filepath;

  char *rvd_string_printable = strdup(hashmap_get(talesf_kwargs, "rvd_string"));
  str_replace(rvd_string_printable, " ", "_");

  char *rvd_string2_printable = strdup(hashmap_get(talesf_kwargs, "rvd_string2"));
  str_replace(rvd_string2_printable, " ", "_");

  //create_options_string(options_str, rvd_str);

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

  if (!gff_out_file || !tab_out_file || (is_genome && !genome_browser_file)) {

    fprintf(log_file, "Error: unable to open output file '%s'\n", output_filepath);

    // TODO: find a better way to do this
    free(rvd_string_printable);
    free(rvd_string2_printable);

    return 1;

  }

  //fprintf(tab_out_file, options_str);


  // Tab File Header
  fprintf(tab_out_file, "RVD1 Best Possible Score:%.2lf\n", best_score);
  fprintf(tab_out_file, "RVD2 Best Possible Score:%.2lf\n", best_score2);
  fprintf(tab_out_file, "Sequence Name\tTAL 1\tTAL 2\tTAL 1 Score\tTAL 2 Score\tTAL 1 Start\tTAL 2 Start\tSpacer Length\tTAL 1 Target\tTAL 2 Target\n");


  // GFF file header
  fprintf(gff_out_file, "##gff-version 3\n");

  fprintf(gff_out_file, "#table_display_tags:tal1,tal2,tal1_target,tal2_target,tal1_score,tal2_score,spacer_length\n");

  fprintf(gff_out_file, "#RVD1 Best Possible Score:%.2lf\n", best_score);
  fprintf(gff_out_file, "#RVD2 Best Possible Score:%.2lf\n", best_score2);

  // Browser file header

  if(genome_using_gbrowse) {
    fprintf(genome_browser_file, "##gff-version 3\n");
  } else if(is_genome) {
    fprintf(genome_browser_file, "track name=\"TAL Targets\" description=\"Targets for TALEN Pair '%s' '%s'\" visibility=2 useScore=1\n", rvd_string_printable, rvd_string2_printable);
  }


  for (unsigned long i = 0; i < array_size(results); i++) {

    BindingSite *site = (BindingSite *) array_get(results, i);

    char tal1_name[16];
    char tal2_name[16];

    int tal2_seq_len = array_size(rvd_seqs[site->r_idx]) + 2;
    char *tal2_sequence = calloc(tal2_seq_len + 1, sizeof(char));

    sprintf(tal1_name, "RVD%d", site->f_idx + 1);
    sprintf(tal2_name, "RVD%d", site->r_idx + 1);

    for (int j = 0; j < tal2_seq_len; j++) {
      char base = site->sequence[1][tal2_seq_len - j - 1];
      if (base == 'A' || base == 'a')
        tal2_sequence[j] = 'T';
      else if (base == 'C' || base == 'c')
        tal2_sequence[j] = 'G';
      else if (base == 'G' || base == 'g')
        tal2_sequence[j] = 'C';
      else if (base == 'T' || base == 't')
        tal2_sequence[j] = 'A';
      else if (base == ' ')
        tal2_sequence[j] = ' ';
      else {
        fprintf(stderr, "Error: unexpected character '%d'\n", base);

        // TODO: find a better way to do this
        free(tal2_sequence);
        free(rvd_string_printable);
        free(rvd_string2_printable);
        fclose(tab_out_file);

        return 1;

      }

    }

    unsigned long end_pos = site->indexes[0] + array_size(rvd_seqs[0]) + site->spacer_length + array_size(rvd_seqs[1]);
    double combined_score = site->scores[0] + site->scores[1];

    fprintf( tab_out_file, "%s\t%s\t%s\t%.2lf\t%.2lf\t%lu\t%lu\t%d\t%s\t%s\n",
             site->sequence_name, tal1_name, tal2_name, site->scores[0], site->scores[1], site->indexes[0], site->indexes[1], site->spacer_length, site->sequence[0], tal2_sequence);

    fprintf( gff_out_file, "%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\ttal1=%s;tal2=%s;tal1_target=%s;tal2_target=%s;tal1_score=%.2lf;tal2_score=%.2lf;spacer_length=%d;\n",
             site->sequence_name, source_str, "TAL_effector_binding_site", site->indexes[0] + 1,
             end_pos, combined_score, '.', tal1_name, tal2_name, site->sequence[0], tal2_sequence, site->scores[0], site->scores[1], site->spacer_length);

    if(is_genome && i < 10000) {

      if(genome_using_gbrowse) {

        fprintf( genome_browser_file, "chr%s\t%s\t%s\t%lu\t%lu\t%.2lf\t%c\t.\tName=site%lu;\n",
                 site->sequence_name, source_str, "TAL_effector_binding_site", site->indexes[0] + 1,
                 end_pos, combined_score, '.', i);

      } else {

        int bed_score = floorf((combined_best_score / combined_score * 1000) + 0.5);
        fprintf( genome_browser_file,"%s\t%lu\t%lu\tsite%lu\t%d\t%c\n",
                 site->sequence_name, site->indexes[0], end_pos - 1, i, bed_score, '+');

      }

    }

    free(tal2_sequence);

  }

  free(rvd_string_printable);
  free(rvd_string2_printable);
  fclose(tab_out_file);
  fclose(gff_out_file);
  if(is_genome) {
    fclose(genome_browser_file);
  }

  return 0;

}

double score_binding_site(kseq_t *seq, unsigned long i, Array *rvd_seq, Hashmap *diresidue_scores, double *lookahead_array, int reverse) {

  double total_score = 0.0;
  int num_rvds = array_size(rvd_seq);

  if (!reverse) {

    for (unsigned long j = 0; j < array_size(rvd_seq); j++) {

      char *rvd = array_get(rvd_seq, j);
      double *scores = hashmap_get(diresidue_scores, rvd);

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

    for (unsigned long j = 0; j < array_size(rvd_seq); j++) {

      char *rvd = array_get(rvd_seq, j);
      double *scores = hashmap_get(diresidue_scores, rvd);

      unsigned long k = i + array_size(rvd_seq) - j - 1;

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

BindingSite *create_binding_site(kseq_t *seq, unsigned long i, unsigned long j, int num_forward_rvds, double forward_score, int num_reverse_rvds, double reverse_score, int spacer_size, int f_idx, int r_idx) {

  int seq_name_len = strlen(seq->name.s);

  BindingSite *site = malloc(sizeof(BindingSite));

  site->sequence_name = calloc(seq_name_len + 1, sizeof(char));
  site->sequence_name[seq_name_len] = '\0';
  strncpy(site->sequence_name, seq->name.s, seq_name_len);

  site->spacer_length = spacer_size;

  site->f_idx = f_idx;
  site->r_idx = r_idx;

  // Plus

  site->indexes[0] = i;

  site->sequence[0] = calloc(num_forward_rvds + 2 + 1, sizeof(char));
  site->sequence[0][num_forward_rvds + 2] = '\0';

  strncpy(site->sequence[0], seq->seq.s + i - 1, 1);

  // Upstream
  site->sequence[0][1] = ' ';
  strncpy(site->sequence[0] + 2, seq->seq.s + i, num_forward_rvds);
  
  for(int k = 0; k < num_forward_rvds + 2 + 1; k++) {
    site->sequence[0][k] = toupper(site->sequence[0][k]);
  }

  site->scores[0] = forward_score;

  // Minus

  site->indexes[1] = j;

  site->sequence[1] = calloc(num_reverse_rvds + 2 + 1, sizeof(char));
  site->sequence[1][num_reverse_rvds + 2] = '\0';

  strncpy(site->sequence[1], seq->seq.s + j - num_reverse_rvds + 1, num_reverse_rvds);

  // Upstream
  site->sequence[1][num_reverse_rvds] = ' ';
  strncpy(site->sequence[1] + num_reverse_rvds + 1, seq->seq.s + j + 1, 1);
  
  for(int k = 0; k < num_reverse_rvds + 2 + 1; k++) {
    site->sequence[1][k] = toupper(site->sequence[1][k]);
  }
  
  site->scores[1] = reverse_score;

  return site;

}

// Identify and print out TAL effector binding sites
void find_binding_sites(FILE *log_file, kseq_t *seq, Array **rvd_seqs, Hashmap *diresidue_scores, double **lookahead_arrays, Array *results) {

  int c_upstream = *((int *) hashmap_get(talesf_kwargs, "c_upstream"));
  int spacer_min = *((int *) hashmap_get(talesf_kwargs, "spacer_min"));
  int spacer_max = *((int *) hashmap_get(talesf_kwargs, "spacer_max"));
  
  int count_only = *((int *) hashmap_get(talesf_kwargs, "count_only"));

  int **count_results_array = NULL;
  
  if (count_only) {
    count_results_array = hashmap_get(talesf_kwargs, "count_results_array");
  }

  if (array_size(rvd_seqs[0]) + array_size(rvd_seqs[0]) + spacer_min > seq->seq.l) {
    logger(log_file, "Warning: skipping sequence '%s' since it is shorter than the RVD sequence\n", seq->seq.s);
    return;
  }

  logger(log_file, "Scanning %s for binding sites (length %ld)", seq->name.s, seq->seq.l);

  for (int f_idx = 0; f_idx < 2; f_idx++) {
    
    for (int r_idx = 0; r_idx < 2; r_idx++) {
        
      Array *forward_rvd_seq = rvd_seqs[f_idx];
      Array *reverse_rvd_seq = rvd_seqs[r_idx];
      
      double *forward_lookahead = lookahead_arrays[f_idx];
      double *reverse_lookahead = lookahead_arrays[r_idx];
        
      int num_forward_rvds = array_size(forward_rvd_seq);
      int num_reverse_rvds = array_size(reverse_rvd_seq);

      for (unsigned long i = 1; i <= seq->seq.l - num_forward_rvds; i++) {
        
        char forward_upstream = seq->seq.s[i-1];
        
        int forward_upstream_is_c = (forward_upstream == 'C' || forward_upstream == 'c');
        int forward_upstream_is_t = (forward_upstream == 'T' || forward_upstream == 't');

        if ((c_upstream != 0 && forward_upstream_is_c) || (c_upstream != 1 && forward_upstream_is_t)) {
          
          double forward_score = score_binding_site(seq, i, forward_rvd_seq, diresidue_scores, forward_lookahead, 0);
          
          if (forward_score != -1) {
              
            for (int spacer_size = spacer_min; spacer_size < spacer_max + 1; spacer_size++) {
                
               unsigned long j = i + num_forward_rvds + spacer_size + num_reverse_rvds - 1;
               
               if (j >= seq->seq.l - 2) continue;
               
               char reverse_upstream = seq->seq.s[j + 1];
               
               if ((c_upstream != 0 && (forward_upstream_is_c && (reverse_upstream == 'G' || reverse_upstream == 'g'))) || (c_upstream != 1 && (forward_upstream_is_t && (reverse_upstream == 'A' || reverse_upstream == 'a')))) {
                   
                 double reverse_score = score_binding_site(seq, j - num_reverse_rvds + 1, reverse_rvd_seq, diresidue_scores, reverse_lookahead, 1);
                 
                 if (reverse_score != -1) {
                   
                   if (count_only) {
                     
                     count_results_array[f_idx][r_idx]++;

//                     if (count_results_array[f_idx][r_idx] >= 5) {
//                       return;
//                     }
                     
                   } else {
                   
                     BindingSite *site = create_binding_site(seq, i, j, num_forward_rvds, forward_score, num_reverse_rvds, reverse_score, spacer_size, f_idx, r_idx);
                   
                     #pragma omp critical (add_result)
                     array_add(results, site);
                   
                   }
                   
                 }
                   
               }
                
            }
            
          }
          
        }
        
      }
      
    }
    
  }

}

int run_paired_talesf_task(Hashmap *kwargs) {

  talesf_kwargs = kwargs;

  // Options
  char *seq_filename = hashmap_get(kwargs, "seq_filename");
  char *rvd_string = hashmap_get(kwargs, "rvd_string");
  char *rvd_string2 = hashmap_get(kwargs, "rvd_string2");
  char *log_filepath = hashmap_get(kwargs, "log_filepath");

  double weight = *((double *) hashmap_get(kwargs, "weight"));
  double cutoff = *((double *) hashmap_get(kwargs, "cutoff"));

  int numprocs = *((int *) hashmap_get(kwargs, "num_procs"));

  int count_only = *((int *) hashmap_get(kwargs, "count_only"));

  // Setup the logger

  FILE *log_file = stdout;

  if (log_filepath && strcmp(log_filepath, "NA") != 0) {
    log_file = fopen(log_filepath, "a");
  }

  // Determine number of sequences in the input file

  int seq_num;
  char cmd[256], line[32];

  sprintf(cmd, "grep '^>' %s | wc -l", seq_filename);
  FILE *fasta_filesize_in = popen(cmd, "r");

  if (!fasta_filesize_in) {
    perror("Error: unable to check fasta file size\n");
    logger(log_file, "Error: unable to check fasta file size");
    if (log_file != stdout) fclose(log_file);
    return 1;
  }

  fgets(line, sizeof(line), fasta_filesize_in);
  pclose(fasta_filesize_in);
  seq_num = atoi(line);

  // Process RVD sequences

  Array *results = array_new( sizeof(BindingSite *) );

  Array *rvd_seq = rvd_string_to_array(rvd_string);
  Array *rvd_seq2 = rvd_string_to_array(rvd_string2);

  Array *rvd_seqs[2];

  rvd_seqs[0] = rvd_seq;
  rvd_seqs[1] = rvd_seq2;

  Array *joined_rvd_seq = array_concat(rvd_seq, rvd_seq2);

  // Get RVD/bp matching scores

  Hashmap *diresidue_probabilities = get_diresidue_probabilities(joined_rvd_seq, weight);
  Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities);
  hashmap_delete(diresidue_probabilities, NULL);

  // Compute optimal score for the RVD sequences

  double best_score = get_best_score(rvd_seq, diresidue_scores);
  double best_score2 = get_best_score(rvd_seq2, diresidue_scores);

  // Define score cutoffs for match sites

  double *lookahead_arrays[2];

  lookahead_arrays[0] = create_lookahead_array(rvd_seq, cutoff, best_score, diresidue_scores);
  lookahead_arrays[1] = create_lookahead_array(rvd_seq2, cutoff, best_score2, diresidue_scores);

  // Begin processing

  int abort = 0;

  omp_set_num_threads(numprocs);

  #pragma omp parallel
  {

    // Open sequence file
    gzFile seqfile = gzopen(seq_filename, "r");

    if (!seqfile) {

      logger(log_file, "Error: unable to open sequence '%s'", seq_filename);
      abort = 1;

    } else {

      kseq_t *seq = kseq_init(seqfile);

      int j = 0;

      #pragma omp for schedule(static)
      for (int i = 0; i < seq_num; i++) {

        #pragma omp flush (abort)
        if (!abort) {

          while (j <= i) {

            int result = kseq_read(seq);

            if (result < 0) {
              logger(log_file, "Error: problem parsing data from '%s'", seq_filename);
              abort = 1;
            }

            j++;

          }

          if (!abort) {
            find_binding_sites(log_file, seq, rvd_seqs, diresidue_scores, lookahead_arrays, results);
          }

        }

      }

      kseq_destroy(seq);
      gzclose(seqfile);

    }

  }

  if (!abort) {

    if (!count_only) {

      qsort(results->data, array_size(results), sizeof(BindingSite *), binding_site_compare_pos);

      abort = print_results(results, rvd_seqs, best_score, best_score2, log_file);

    }

    logger(log_file, "Finished");

  }

  // Free memory

  if (results) {

    for (int i = 0; i < array_size(results); i++) {

      BindingSite *site = (BindingSite *) array_get(results, i);

      free(site->sequence[0]);
      free(site->sequence[1]);
      free(site->sequence_name);
      free(site);

    }

    array_delete(results, NULL);

  }

  free(lookahead_arrays[0]);
  free(lookahead_arrays[1]);

  if (rvd_seq) {
    array_delete(rvd_seq, free);
  }

  if (rvd_seq2) {
    array_delete(rvd_seq2, free);
  }

  if (joined_rvd_seq) {
    array_delete(joined_rvd_seq, NULL);
  }

  if (diresidue_scores) {
    hashmap_delete(diresidue_scores, free);
  }

  if (log_file != stdout) {
    fclose(log_file);
  }

  return abort;

}
