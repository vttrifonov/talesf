#include <getopt.h>
#include <stdio.h>
#include <omp.h>

#include "pairedtalesf.h"

#include <bcutils/Hashmap.h>
#include <bcutils/Array.h>
#include <bcutils/bcutils.h>

#define BIGGEST_RVD_SCORE_EVER 100

// Print usage statement
void print_usage(FILE *out_stream, char *prog_name)
{
  fprintf( out_stream, "\nUsage: %s [options] sequence_file_path \"rvdseq\" \"rvdseq2\"\n"
           "  Options:\n"
           "    -c|--cupstream        sets the allowed upstream bases; 0 for T only, 1 for C only, 2 for either\n"
           "    -h|--help             print this help message and exit\n"
           "    -m|--min              the minimum allowed spacer size; default is 15\n"
           "    -n|--numprocs         the number of processors to use; default is 1\n"
           "    -o|--outfile          template filename to which output will be written; both a tab-delimited file "
           "                          and gff3 file will be produced \n"
           "    -w|--weight           user-defined weight; default is 0.9\n"
           "    -x|--max              the maximum allowed spacer size; default is 30\n"
           "    -t|--cutoffmult       multiple of best score at which potential sites will be\n"
           "                          filtered; default is 3.0\n\n", prog_name );
}

int main(int argc, char **argv)
{

  char *prog_name;
  char *seq_filepath;
  char *rvd_string;
  char *rvd_string2;
  char *log_filepath;
  int num_procs;
  char out_filepath[256];
  int c_upstream;
  double weight;
  double cutoff;
  int min;
  int max;

  // Set Defaults

  num_procs = 1;
  weight = 0.9;
  cutoff = 3.0;
  log_filepath = "NA";
  c_upstream = 0;
  min = 15;
  max = 30;

  prog_name = argv[0];

  int opt, opt_index;
  const char *opt_str = "c:hn:o:s:w:t:m:x:";
  const struct option otsf_options[] =
  {
    { "help", no_argument, NULL, 'h' },
    { "numprocs", required_argument, NULL, 'n' },
    { "outfile", required_argument, NULL, 'o' },
    { "weight", required_argument, NULL, 'w' },
    { "cutoffmult", required_argument, NULL, 't' },
    { "cupstream", required_argument, NULL, 'c' },
    { "min", required_argument, NULL, 'm' },
    { "max", required_argument, NULL, 'x' },
    { NULL, no_argument, NULL, 0 },
  };

  for( opt = getopt_long(argc, argv + 0, opt_str, otsf_options, &opt_index);
       opt != -1;
       opt = getopt_long(argc, argv + 0, opt_str, otsf_options, &opt_index) )
  {
    switch(opt)
    {

      case 'c':
        if( sscanf(optarg, "%d", &c_upstream) != 1 )
        {
          fprintf(stderr, "Error: unable to convert cupstream '%s' to an integer\n", optarg);
          return 1;
        }

        if ( c_upstream != 0 && c_upstream != 1 && c_upstream != 2) {
          fprintf(stderr, "Error: cupstream must be 0, 1, or 2\n");
          return 1;
        }

        break;

      case 'h':
        print_usage(stdout, prog_name);
        return 0;

      case 'n':
        if( sscanf(optarg, "%d", &num_procs) != 1 )
        {
          fprintf(stderr, "Error: unable to convert numprocs '%s' to an integer\n", optarg);
          return 1;
        }
        if( num_procs > omp_get_num_procs())
        {
          fprintf(stderr, "Error: numprocs was %d but only %d are available\n", num_procs, omp_get_num_procs());
          return 1;
        }
        break;

      case 'm':
        if( sscanf(optarg, "%d", &min) != 1 )
        {
          fprintf(stderr, "Error: unable to convert min '%s' to an integer\n", optarg);
          return 1;
        }
        break;

      case 'x':
        if( sscanf(optarg, "%d", &max) != 1 )
        {
          fprintf(stderr, "Error: unable to convert max '%s' to an integer\n", optarg);
          return 1;
        }
        break;

      case 'o':
        strcpy(out_filepath, optarg);
        break;

      case 'w':
        if( sscanf(optarg, "%lf", &weight) != 1 )
        {
          fprintf(stderr, "Error: unable to convert weight '%s' to a double\n", optarg);
          return 1;
        }
        break;

      case 't':
        if( sscanf(optarg, "%lf", &cutoff) != 1 )
        {
          fprintf(stderr, "Error: unable to convert cutoff multiple '%s' to a double\n", optarg);
          return 1;
        }
        break;
    }
  }

  // Parse arguments
  if(argc - optind != 3)
  {
    fputs("Error: must provide sequence (file path) and 2 RVD sequences (strings)\n", stderr);
    print_usage(stderr, prog_name);
    return 1;
  }

  seq_filepath = argv[optind];
  rvd_string = argv[optind + 1];
  rvd_string2 = argv[optind + 2];

  Hashmap *talesf_kwargs = hashmap_new(32);

  Array *rvd_array = rvd_string_to_array(rvd_string);
  Array *rvd_array2 = rvd_string_to_array(rvd_string2);

  Array *joined_rvd_array = array_concat(rvd_array, rvd_array2);

  // Get RVD/bp matching scores

  Hashmap *diresidue_probabilities = get_diresidue_probabilities(joined_rvd_array, weight);
  Hashmap *diresidue_scores = convert_probabilities_to_scores(diresidue_probabilities);
  hashmap_delete(diresidue_probabilities, NULL);
  
  // Convert hashmap to int map
  
  hashmap_add(diresidue_scores, "XX", double_array(0, 0, 0, 0, BIGGEST_RVD_SCORE_EVER));
  
  double **scoring_matrix = calloc(hashmap_size(diresidue_scores), sizeof(double*));
  
  Hashmap *rvd_to_int = hashmap_new(hashmap_size(diresidue_scores));
  unsigned int *rvd_ints = calloc(hashmap_size(diresidue_scores), sizeof(unsigned int));
  
  char **diresidues = hashmap_keys(diresidue_scores);
  
  for (unsigned int i = 0; i < hashmap_size(diresidue_scores); i++) {
  
    rvd_ints[i] = i;
    hashmap_add(rvd_to_int, diresidues[i], rvd_ints + i);
  
    scoring_matrix[i] = hashmap_get(diresidue_scores, diresidues[i]);
    scoring_matrix[i][4] = BIGGEST_RVD_SCORE_EVER;
  
  }
  
  // Convert RVD seqs to int seqs
  unsigned int *rvd_seqs[2];
  rvd_seqs[0] = (unsigned int*) calloc(array_size(rvd_array), sizeof(unsigned int));
  rvd_seqs[1] = (unsigned int*) calloc(array_size(rvd_array2), sizeof(unsigned int));
  
  for (unsigned int i = 0; i < array_size(rvd_array); i++) {
    rvd_seqs[0][i] = *(unsigned int *)(hashmap_get(rvd_to_int, array_get(rvd_array, i)));
  }
  
  for (unsigned int i = 0; i < array_size(rvd_array2); i++) {
    rvd_seqs[1][i] = *(unsigned int *)(hashmap_get(rvd_to_int, array_get(rvd_array2, i)));
  }
  
  unsigned int rvd_seqs_lens[2];
  rvd_seqs_lens[0] = array_size(rvd_array);
  rvd_seqs_lens[1] = array_size(rvd_array2);
  
  // Compute optimal scores for the RVD sequences
  double best_scores[2];
  best_scores[0] = get_best_score(rvd_array, diresidue_scores);
  best_scores[1] = get_best_score(rvd_array2, diresidue_scores);
  
  hashmap_add(talesf_kwargs, "seq_filename", seq_filepath);
  hashmap_add(talesf_kwargs, "rvd_seqs", rvd_seqs);
  hashmap_add(talesf_kwargs, "rvd_seqs_lens", rvd_seqs_lens);
  hashmap_add(talesf_kwargs, "rvd_string", rvd_string);
  hashmap_add(talesf_kwargs, "rvd_string2", rvd_string2);
  hashmap_add(talesf_kwargs, "best_scores", best_scores);
  hashmap_add(talesf_kwargs, "scoring_matrix", scoring_matrix);
  hashmap_add(talesf_kwargs, "output_filepath", out_filepath);
  hashmap_add(talesf_kwargs, "log_filepath", log_filepath);
  hashmap_add(talesf_kwargs, "weight", &weight);
  hashmap_add(talesf_kwargs, "cutoff", &cutoff);
  hashmap_add(talesf_kwargs, "c_upstream", &c_upstream);
  hashmap_add(talesf_kwargs, "spacer_min", &min);
  hashmap_add(talesf_kwargs, "spacer_max", &max);
  hashmap_add(talesf_kwargs, "num_procs", &num_procs);
  hashmap_add(talesf_kwargs, "organism_name", "");
  
  int count_only = 0;
  hashmap_add(talesf_kwargs, "count_only", &count_only);

  int task_result = run_paired_talesf_task(talesf_kwargs);

  hashmap_delete(talesf_kwargs, NULL);
  
  if (rvd_seqs[0]) {
    free(rvd_seqs[0]);
  }
  
  if (rvd_seqs[1]) {
    free(rvd_seqs[1]);
  }
  
  if (scoring_matrix) {
    free(scoring_matrix);
  }
  
  if (rvd_to_int) {
    hashmap_delete(rvd_to_int, NULL);
  }
  
  if (rvd_ints) {
    free(rvd_ints);
  }
  
  if (diresidues) {
    free(diresidues);
  }
  
  if (rvd_array) {
    array_delete(rvd_array, free);
  }

  if (rvd_array2) {
    array_delete(rvd_array2, free);
  }

  if (joined_rvd_array) {
    array_delete(joined_rvd_array, NULL);
  }

  if (diresidue_scores) {
    hashmap_delete(diresidue_scores, free);
  }

  return task_result;

}
