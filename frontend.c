#include <getopt.h>
#include <stdio.h>

#include "talesf.h"
#include "Hashmap.h"

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
  const char *opt_str = "c:fhn:o:s:w:x:";
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

  hashmap_add(talesf_kwargs, "seq_filename", seq_filepath);
  hashmap_add(talesf_kwargs, "rvd_string", rvd_string);
  hashmap_add(talesf_kwargs, "output_filepath", out_filepath);
  hashmap_add(talesf_kwargs, "log_filepath", log_filepath);
  hashmap_add(talesf_kwargs, "weight", &weight);
  hashmap_add(talesf_kwargs, "cutoff", &cutoff);
  hashmap_add(talesf_kwargs, "forward_only", &forward_only);
  hashmap_add(talesf_kwargs, "c_upstream", &c_upstream);
  hashmap_add(talesf_kwargs, "num_procs", &num_procs);
  hashmap_add(talesf_kwargs, "organism_name", "");

  int task_result = run_talesf_task(talesf_kwargs);

  hashmap_delete(talesf_kwargs, NULL);

  return task_result;

}
