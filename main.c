#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include "demultiplex.h"

#define VERSION "0.1"

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
  fprintf(stderr, " -d  FILE    dual barcode pair file\n");
  fprintf(stderr, " -p  CHAR    output path\n");
  fprintf(stderr, " -l  FILE    output log file\n");
  fprintf(stderr, " -M  CHAR    output mode, fastq or fasta [default fastq]\n");
  fprintf(stderr, " -B  NUM     batch size [default 400M]\n");
  fprintf(stderr, " -t  INT     number of threads [default 3]\n");
  fprintf(stderr, " -L  INT     target length for detecting barcode [default 150]\n");
  fprintf(stderr, " -m  INT     match score [default 2]\n");
  fprintf(stderr, " -x  INT     mismatch score [default 2]\n");
  fprintf(stderr, " -o  INT     gap open score [default 3]\n");
  fprintf(stderr, " -e  INT     gap extension score [default 1]\n");
  fprintf(stderr, " -s  INT     minimal alignment score for demultiplexing\n");
  fprintf(stderr, " -i          ignore parameter estimation\n");
  fprintf(stderr, " -h          help information\n");
  fprintf(stderr, " -v          show version number\n\n");
  fprintf(stderr, "-b -p must be specified.\n\n");
  fprintf(stderr, "Example:\n");
  fprintf(stderr, "nanoplexer -b barcode.fa -p /ouput/ input.fastq\n\n");
  return 0;
}

int main(int argc, char *argv[])
{
  const char *optstring = "b:d:p:l:M:B:t:L:m:x:o:e:s:ihv";
  static struct option long_options[] = {
    {"barcode", required_argument, 0, 'b'},
    {"dual", optional_argument, 0, 'd'},
    {"path", required_argument, 0, 'p'},
    {"log", optional_argument, 0, 'l'},
    {"mode", optional_argument, 0, 'M'},
    {"batch", optional_argument, 0, 'B'},
    {"threads", optional_argument, 0, 't'},
    {"length", optional_argument, 0, 'L'},
    {"match", optional_argument, 0, 'm'},
    {"mismatch", optional_argument, 0, 'x'},
    {"open", optional_argument, 0, 'o'},
    {"extension", optional_argument, 0, 'e'},
    {"score", optional_argument, 0, 's'},
    {"ignore", no_argument, 0, 'i'},
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
      case 'd': opt.dual = strdup(optarg); break;
      case 'p': opt.path = strdup(optarg); break;
      case 'l': opt.log = fopen(optarg, "w"); break;
      case 'M':
        if (strcmp(optarg, "fasta") == 0 || strcmp(optarg, "fa") == 0)  opt.mode = 1;
        break;
      case 'B': opt.chunk_size = (int)mm_parse_num(optarg); break;
      case 't': opt.n_threads = atoi(optarg); break;
      case 'L': opt.LEN = atoi(optarg); opt.MASK_LEN = opt.LEN / 2; break;
      case 'm': opt.match = atoi(optarg); break;
      case 'x': opt.mismatch = atoi(optarg); break;
      case 'o': opt.open = atoi(optarg); break;
      case 'e': opt.ext = atoi(optarg); break;
      case 's': opt.score = atoi(optarg); break;
      case 'i': opt.flag = 1; break;
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
  if (!opt.path) {
    print_log("Please specify output path.\n");
    return usage();
  }
  if (opt.flag && !opt.score) {
    print_log("Please specify alignment score if ignoring parameter estimation.\n");
    return usage();
  };

  opt_set(&opt, optind < argc ? argv[optind] : 0);
  demultiplex_data(&opt);
  opt_free(&opt);

  print_log("Finished demultiplexing sequence data.");

  return 0;
}