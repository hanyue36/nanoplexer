#ifndef DEMULTIPLEX_H
#define DEMULTIPLEX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(dual, int)

#define ERR 0.15

#define BUFSIZE 110000000U
#define WRITESIZE 100000000U

typedef struct 
{
  char **name;
  int8_t **nt, **nt_rc, *len;
  int idx, capacity;

  FILE **ptr;
  char **buffer;
  int32_t *offset;
  int file_num;

  khash_t(dual) *h;
} bc_t;

typedef struct 
{
  char *fb, *path, *dual;

  int n_threads, chunk_size, mode;
  int LEN, MASK_LEN;

  int match, mismatch;
  int open, ext;
  int8_t mat[25];
  int score, flag;

  FILE *log;
  kseq_t *ks;
  bc_t *bc;
} opt_t;

void opt_init(opt_t *opt);

void opt_set(opt_t *opt, char *fn);

void opt_free(opt_t *opt);

void print_log(char *format, ...);

void rc_seq(char *seq, int len);

int8_t *seq_to_nt(char *seq, int len, int flag);

void demultiplex_data(opt_t *opt);

#endif