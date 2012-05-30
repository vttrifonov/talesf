#include <getopt.h>
#include <stdio.h>
#include <omp.h>

#include "pairedtalesf.h"
#include "Hashmap.h"

// Print usage statement
void print_usage(FILE *out_stream, char *prog_name)
{
  fprintf( out_stream, "\nUsage: %s [options] sequence_file_path \"rvdseq\"\n"
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

  hashmap_add(talesf_kwargs, "seq_filename", seq_filepath);
  hashmap_add(talesf_kwargs, "rvd_string", rvd_string);
  hashmap_add(talesf_kwargs, "rvd_string2", rvd_string2);
  hashmap_add(talesf_kwargs, "output_filepath", out_filepath);
  hashmap_add(talesf_kwargs, "log_filepath", log_filepath);
  hashmap_add(talesf_kwargs, "weight", &weight);
  hashmap_add(talesf_kwargs, "cutoff", &cutoff);
  hashmap_add(talesf_kwargs, "c_upstream", &c_upstream);
  hashmap_add(talesf_kwargs, "num_procs", &num_procs);
  hashmap_add(talesf_kwargs, "spacer_min", &min);
  hashmap_add(talesf_kwargs, "spacer_max", &max);
  hashmap_add(talesf_kwargs, "organism_name", "");

  int task_result = run_paired_talesf_task(talesf_kwargs);

  hashmap_delete(talesf_kwargs, NULL);

  return task_result;

}
