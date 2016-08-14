#include <getopt.h>
#include <stdio.h>
#include <omp.h>

#include "talesf.h"

#include <bcutils/Hashmap.h>
#include <bcutils/Array.h>
#include <bcutils/bcutils.h>

#define BIGGEST_RVD_SCORE_EVER 100

// Print usage statement
void print_usage(FILE *out_stream, char *prog_name)
{
  fprintf( out_stream, "\nUsage: %s [options] sequence_file_path \"rvdseq\"\n"
           "  Options:\n"
           "    -c|--cupstream        sets the allowed upstream bases; 0 for T only, 1 for C only, 2 for either\n"
           "    -f|--forwardonly      only search the forward strand of the sequence\n"
           "    -h|--help             print this help message and exit\n"
           "    -n|--numprocs         the number of processors to use; default is 1\n"
           "    -o|--outfile          template filename to which output will be written; both a tab-delimited file "
           "                          and gff3 file will be produced \n"
           "    -w|--weight           user-defined weight; default is 0.9\n"
           "    -x|--cutoffmult       multiple of best score at which potential sites will be\n"
           "                          filtered; default is 3.0\n\n", prog_name );
}

int main(int argc, char **argv)
{

  char *prog_name;
  char *seq_filepath;
  char *rvd_string;
  char *log_filepath;
  int forward_only;
  int num_procs;
  char out_filepath[256];
  int c_upstream;
  double weight;
  double cutoff;

  // Set Defaults

  forward_only = 0;
  num_procs = 1;
  weight = 0.9;
  cutoff = 3.0;
  log_filepath = "NA";
  c_upstream = 0;

  prog_name = argv[0];

  int opt, opt_index;
  const char *opt_str = "c:gfhn:o:s:w:x:";
  const struct option otsf_options[] =
  {
    { "forwardonly", no_argument, NULL, 'f' },
    { "help", no_argument, NULL, 'h' },
    { "numprocs", required_argument, NULL, 'n' },
    { "outfile", required_argument, NULL, 'o' },
    { "weight", required_argument, NULL, 'w' },
    { "cutoffmult", required_argument, NULL, 'x' },
    { "cupstream", required_argument, NULL, 'c' },
    { NULL, no_argument, NULL, 0 },
  };

  for( opt = getopt_long(argc, argv + 0, opt_str, otsf_options, &opt_index);
       opt != -1;
       opt = getopt_long(argc, argv + 0, opt_str, otsf_options, &opt_index) )
  {
    switch(opt)
    {
      case 'f':
        forward_only = 1;
        break;

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

      case 'x':
        if( sscanf(optarg, "%lf", &cutoff) != 1 )
        {
          fprintf(stderr, "Error: unable to convert cutoff multiple '%s' to a double\n", optarg);
          return 1;
        }
        break;
    }
  }

  // Parse arguments
  if(argc - optind != 2)
  {
    fputs("Error: must provide sequence (file) and RVD sequence (string)\n", stderr);
    print_usage(stderr, prog_name);
    return 1;
  }

  seq_filepath = argv[optind];
  rvd_string = argv[optind + 1];
  
  Hashmap *talesf_kwargs = hashmap_new(32);
  
  Array *rvd_array = rvd_string_to_array(rvd_string);
  
  // Get RVD/bp matching scores

  Hashmap *diresidue_probabilities = get_diresidue_probabilities(rvd_array, weight);
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
  
  unsigned int *rvd_seq = (unsigned int*) calloc(array_size(rvd_array), sizeof(unsigned int));
  
  for (unsigned int i = 0; i < array_size(rvd_array); i++) {
    rvd_seq[i] = *(unsigned int *)(hashmap_get(rvd_to_int, array_get(rvd_array, i)));
  }
  
  unsigned int rvd_seq_len = array_size(rvd_array);
  
  double best_score = get_best_score(rvd_array, diresidue_scores);
  
  int scoring_matrix_length = hashmap_size(diresidue_scores);
  
  hashmap_add(talesf_kwargs, "seq_filename", seq_filepath);
  hashmap_add(talesf_kwargs, "rvd_seq", rvd_seq);
  hashmap_add(talesf_kwargs, "rvd_seq_len", &rvd_seq_len);
  hashmap_add(talesf_kwargs, "rvd_string", rvd_string);
  hashmap_add(talesf_kwargs, "best_score", &best_score);
  hashmap_add(talesf_kwargs, "scoring_matrix", scoring_matrix);
  hashmap_add(talesf_kwargs, "scoring_matrix_length", &scoring_matrix_length);
  hashmap_add(talesf_kwargs, "output_filepath", out_filepath);
  hashmap_add(talesf_kwargs, "log_filepath", log_filepath);
  hashmap_add(talesf_kwargs, "weight", &weight);
  hashmap_add(talesf_kwargs, "cutoff", &cutoff);
  hashmap_add(talesf_kwargs, "c_upstream", &c_upstream);
  hashmap_add(talesf_kwargs, "num_procs", &num_procs);
  hashmap_add(talesf_kwargs, "organism_name", "");
  
  hashmap_add(talesf_kwargs, "forward_only", &forward_only);

  int task_result = run_talesf_task(talesf_kwargs);

  hashmap_delete(talesf_kwargs, NULL);
  
  if (rvd_seq) {
    free(rvd_seq);
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

  if (diresidue_scores) {
    hashmap_delete(diresidue_scores, free);
  }

  return task_result;

}
