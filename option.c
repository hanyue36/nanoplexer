#include "demultiplex.h"
#include <unistd.h>
#include <sys/stat.h>

void opt_init(opt_t *opt)
{
  memset(opt, 0, sizeof(opt_t));

  opt->match = 2; opt->mismatch = 2;
  opt->open = 3; opt->ext = 1;

  opt->n_threads = 3;
  opt->chunk_size = 400000000;

  opt->bc = (bc_t *)malloc(sizeof(bc_t));
  memset(opt->bc, 0, sizeof(bc_t));
}

void opt_set_mat(int match, int mismatch, int8_t mat[25])
{
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      mat[k++] = i == j? match : -mismatch;
      mat[k++] = 0;
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
}

void opt_allocate_bc(opt_t *opt)
{
  bc_t *bc = opt->bc;
  bc->capacity = 16;
  bc->name = (char **)malloc(sizeof(char *) * 16);
  bc->len = (int8_t *)malloc(sizeof(int8_t) * 16);
  bc->nt = (int8_t **)malloc(sizeof(int8_t *) * 16);
  bc->nt_rc = (int8_t **)malloc(sizeof(int8_t *) * 16);
  gzFile fp = gzopen(opt->fb, "r");
  kseq_t *seq = kseq_init(fp);
  while(kseq_read(seq) >= 0) {
    if (bc->idx + 1 >= bc->capacity) {
      bc->capacity <<= 1;
      bc->name = (char **)realloc(bc->name, sizeof(char *) * bc->capacity);
      bc->len = (int8_t *)realloc(bc->len, sizeof(int8_t) * bc->capacity);
      bc->nt = (int8_t **)realloc(bc->nt, sizeof(int8_t *) * bc->capacity);
      bc->nt_rc = (int8_t **)realloc(bc->nt_rc, sizeof(int8_t *) * bc->capacity);
    }
    bc->name[bc->idx] = strdup(seq->name.s);
    bc->len[bc->idx] = seq->seq.l;
    bc->nt[bc->idx] = seq_to_nt(seq->seq.s, seq->seq.l, 0);
    bc->nt_rc[bc->idx] = seq_to_nt(seq->seq.s, seq->seq.l, 1);
    bc->idx += 1;
  }
  kseq_destroy(seq);
  gzclose(fp);

  bc->score = min_score(opt);
}

void opt_hash_bc(opt_t *opt)
{
  if (access(opt->path, F_OK)) {
    int status = mkdir(opt->path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status) print_log("Failed to create output dir %s.", opt->path), exit(1);
  }

  int i = 0; char fname[1024];
  khint_t k; int absent;
  bc_t *bc = opt->bc;
  bc->h = kh_init(sample);

  bc->file_num = bc->idx * bc->idx + 1;
  bc->ptr = (FILE **)malloc(sizeof(FILE *) * (bc->file_num));
  bc->buffer = (char **)malloc(sizeof(char *) * (bc->file_num));
  bc->offset = (int32_t *)calloc((bc->file_num), sizeof(int32_t));

  FILE *fp = fopen(opt->fs, "r");
  char *line = NULL; size_t len;
  char fn[128], bc1[128], bc2[128], tmp[256];
  while(getline(&line, &len, fp) != -1) {
    sscanf(line, "%s\t%s\t%s", fn, bc1, bc2);
    sprintf(tmp, "%s%s", bc1, bc2);
    tmp[strlen(bc1) + strlen(bc2)] = '\0';
    k = kh_put(sample, bc->h, tmp, &absent);
    if (absent) {
      sprintf(fname, "%s/%s.fastq", opt->path, fn);
      bc->ptr[i] = fopen(fname, "w");
      bc->buffer[i] = (char *)malloc(sizeof(char) * BUFSIZE);
      kh_key(bc->h, k) = strdup(tmp);
      kh_val(bc->h, k) = i++;
    }
  }
  sprintf(fname, "%s/unclassified.fastq", opt->path);
  bc->ptr[i] = fopen(fname, "w");
  bc->buffer[i] = (char *)malloc(sizeof(char) * BUFSIZE);
  bc->file_num = i;
  free(line);
  fclose(fp);
}

void opt_set_bc(opt_t *opt)
{
  opt_allocate_bc(opt);

  opt_hash_bc(opt);

  print_log("Minimal demultiplex score is %.2f.", opt->bc->score);
}

void opt_set(opt_t *opt, char *fn)
{
  opt_set_mat(opt->match, opt->mismatch, opt->mat);

  gzFile fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
  opt->ks = kseq_init(fp);

  opt_set_bc(opt);
}

void opt_free(opt_t *opt)
{
  int i;
  bc_t *bc = opt->bc;

  for(i = 0; i <= bc->file_num; i++) {
    if (bc->offset[i]) {
      bc->buffer[i][bc->offset[i]] = '\0';
      fprintf(bc->ptr[i], "%s", bc->buffer[i]);
    }
    fclose(bc->ptr[i]);
    free(bc->buffer[i]);
  }
  for(i = 0; i < bc->idx; i++) {
    free(bc->name[i]);
    free(bc->nt[i]); free(bc->nt_rc[i]);
  }

  khint_t k;
  for(k = kh_begin(bc->h); k != kh_end(bc->h); ++k)
    if (kh_exist(bc->h, k)) free((char*)kh_key(bc->h, k));
  kh_destroy(sample, bc->h);
  free(opt->fs);

  free(bc->ptr); free(bc->buffer); free(bc->offset);
  free(bc->name); free(bc->len);
  free(bc->nt); free(bc->nt_rc);
  free(bc);

  gzFile fp = opt->ks->f->f;
  kseq_destroy(opt->ks);
  gzclose(fp);
  free(opt->fb); free(opt->path);
  if (opt->log) fclose(opt->log);
}