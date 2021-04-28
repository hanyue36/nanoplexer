#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include "demultiplex.h"

#define VERSION "0.73939133"

static inline int64_t mm_parse_num(const char *str)
{
  double x;
  char *p;
  x = strtod(str, &p);
  if (*p == 'G' || *p == 'g') x *= 1e9;
  else if (*p == 'M' || *p == 'm') x *= 1e6;
  else if (*p == 'K' || *p == 'k') x *= 1e3;
  return (int64_t)(x + .499);
}

static int usage()
{
  fprintf(stderr, "Usage: nanoplexer [options] input.fastq\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -b  FILE    barcode file\n");
  fprintf(stderr, " -s  FILE    barcode combination file\n");
  fprintf(stderr, " -p  CHAR    output path\n");
  fprintf(stderr, " -l  FILE    output log file\n");
  fprintf(stderr, " -B  NUM     batch size [default 400M]\n");
  fprintf(stderr, " -t  INT     number of threads [default 3]\n");
  fprintf(stderr, " -m  INT     match score [default 2]\n");
  fprintf(stderr, " -x  INT     mismatch score [default 2]\n");
  fprintf(stderr, " -o  INT     gap open score [default 3]\n");
  fprintf(stderr, " -e  INT     gap extension score [default 1]\n");
  fprintf(stderr, " -E  FLOAT   error rate [recommend 0.08 for pb, 0.15 for ont]\n");
  fprintf(stderr, " -h          help information\n");
  fprintf(stderr, " -v          show version number\n\n");
  fprintf(stderr, "Example:\n");
  fprintf(stderr, "nanoplexer -b barcode.fa -s sample.txt -p output_path -E 0.08 input.fastq\n\n");
  return 0;
}

int main(int argc, char *argv[])
{
  const char *optstring = "b:s:p:l:B:t:m:x:o:e:E:hv";
  static struct option long_options[] = {
    {"barcode", required_argument, 0, 'b'},
    {"sample", required_argument, 0, 's'},
    {"path", required_argument, 0, 'p'},
    {"log", optional_argument, 0, 'l'},
    {"batch", optional_argument, 0, 'B'},
    {"threads", optional_argument, 0, 't'},
    {"match", optional_argument, 0, 'm'},
    {"mismatch", optional_argument, 0, 'x'},
    {"open", optional_argument, 0, 'o'},
    {"extension", optional_argument, 0, 'e'},
    {"error", required_argument, 0, 'E'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {0, 0, 0, 0}
  };

  int c;
  opt_t opt;
  opt_init(&opt);

  while((c = getopt_long(argc, argv, optstring, long_options, NULL)) != -1) {
    switch (c) {
      case 'b': opt.fb = strdup(optarg); break;
      case 's': opt.fs = strdup(optarg); break;
      case 'p': opt.path = strdup(optarg); break;
      case 'l': opt.log = fopen(optarg, "w"); break;
      case 'B': opt.chunk_size = (int)mm_parse_num(optarg); break;
      case 't': opt.n_threads = atoi(optarg); break;
      case 'm': opt.match = atoi(optarg); break;
      case 'x': opt.mismatch = atoi(optarg); break;
      case 'o': opt.open = atoi(optarg); break;
      case 'e': opt.ext = atoi(optarg); break;
      case 'E': opt.err = atof(optarg); break;
      case 'h': return usage();
      case 'v': fprintf(stdout, "Version: %s\n", VERSION); return 0;
      default:  return usage();
    }
  }

  if (optind + 1 > argc) return usage();
  if (!opt.fb) {
    print_log("Please specify barcode file.\n");
    return usage();
  }
  if (!opt.fs) {
    print_log("Please specify barcode combination file.\n");
    return usage();
  }
  if (!opt.path) {
    print_log("Please specify output path.\n");
    return usage();
  }
  if (!opt.err) {
    print_log("Please specify estimated error rate for each read.\n");
    return usage();
  }

  print_log("Start demultiplexing sequence data.");

  opt_set(&opt, optind < argc ? argv[optind] : 0);
  demultiplex_data(&opt);
  opt_free(&opt);

  print_log("Finished demultiplexing sequence data.");

  return 0;
}