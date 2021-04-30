#include "ssw.h"
#include "demultiplex.h"
#include "ksort.h"

typedef struct 
{
  int idx, strand, pos, score;
  float ratio;
} aln_t;

#define pair_gt(a, b) ((a).ratio > (b).ratio)
KSORT_INIT(pair, aln_t, pair_gt)
KSORT_INIT_GENERIC(float)

typedef struct 
{
  int idx;
  char *name, *seq, *qual;
  int l_name, l_seq, l_qual;
  aln_t aln1[64], aln2[64], aln3[64];
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
  int size = 0, m, n;
  bseq1_t *seqs;
  m = n = 0; seqs = 0;
  while (kseq_read(ks) >= 0) {
    bseq1_t *s;
    if (n >= m) {
      m = m? m<<1 : 256;
      seqs = realloc(seqs, m * sizeof(bseq1_t));
    }
    s = &seqs[n];
    s->name = strdup(ks->name.s); s->seq = strdup(ks->seq.s);
    s->l_name = ks->name.l; s->l_seq = ks->seq.l;
    if (ks->qual.s)  s->qual = strdup(ks->qual.s), s->l_qual = ks->qual.l;
    size += seqs[n++].l_seq;
    if (size >= chunk_size) break;
  }
  *n_ = n;
  return seqs;
}

void aln_core(char *seq, opt_t *opt, aln_t *aln)
{
  int i, num, len;
  bc_t *bc = opt->bc;

  s_profile *profile;
  s_align *result;

  num = bc->idx;
  len = strlen(seq);
  int8_t *nt = seq_to_nt(seq, len, 0);
  profile = ssw_init(nt, len, opt->mat, 5, 2);

  for(i = 0; i < num; i++) {
    result = ssw_align(profile, bc->nt[i], bc->len[i], opt->open, opt->ext, 0, 0, 0, len / 2);
    aln[i].idx = i; aln[i].strand = 0;
    aln[i].pos = result->read_end1; aln[i].score = result->score1;
    aln[i].ratio = (float)result->score1 / bc->score[i];
    align_destroy(result);
  }
  for(i = 0; i < num; i++) {
    result = ssw_align(profile, bc->nt_rc[i], bc->len[i], opt->open, opt->ext, 0, 0, 0, len / 2);
    aln[num + i].idx = i; aln[num + i].strand = 1;
    aln[num + i].pos = result->read_end1; aln[num + i].score = result->score1;
    aln[num + i].ratio = (float)result->score1 / bc->score[i];
    align_destroy(result);
  }

  ks_mergesort(pair, num * 2, aln, 0);

  init_destroy(profile);
  free(nt);
}

void aln_barcode(opt_t *opt, bseq1_t *seq)
{
  bc_t *bc = opt->bc;
  seq->idx = bc->file_num;
  if (seq->l_seq <= 450) return;

  char *front = (char *)calloc(151, sizeof(char));
  char *middle = (char *)calloc(seq->l_seq - 299, sizeof(char));
  char *rear = (char *)calloc(151, sizeof(char));
  memcpy(front, seq->seq, 150);
  memcpy(middle, seq->seq + 150, seq->l_seq - 300);
  memcpy(rear, seq->seq + seq->l_seq - 150, 150);
  aln_core(front, opt, seq->aln1);
  aln_core(middle, opt, seq->aln2);
  aln_core(rear, opt, seq->aln3);
  free(front); free(middle); free(rear);

  //TODO: signal decrease
  if ((seq->aln1[0].ratio >= 1 && seq->aln3[0].ratio >= 1) \
  && (seq->aln2[0].ratio <= 1)) {
    khint_t k; char tmp[1024];
    sprintf(tmp, "%s%s", bc->name[seq->aln1[0].idx], bc->name[seq->aln3[0].idx]);
    tmp[bc->len[seq->aln1[0].idx] + bc->len[seq->aln3[0].idx]] = '\0';
    k = kh_get(sample, bc->h, tmp);
    if (k != kh_end(bc->h)) seq->idx = kh_val(bc->h, k);
    sprintf(tmp, "%s%s", bc->name[seq->aln3[0].idx], bc->name[seq->aln1[0].idx]);
    tmp[bc->len[seq->aln3[0].idx] + bc->len[seq->aln1[0].idx]] = '\0';
    k = kh_get(sample, bc->h, tmp);
    if (k != kh_end(bc->h)) seq->idx = kh_val(bc->h, k);
  }
}

void log_output(opt_t *opt, aln_t *aln, int flag)
{
  fprintf(opt->log, "%s\t%c\t%d\t%d\t%.2f%c", \
  opt->bc->name[aln->idx], "+-"[aln->strand], \
  aln->pos, aln->score, aln->ratio, "\t\n"[flag]);
}


void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct 
{
  int n_seqs;
  bseq1_t *seqs;
  opt_t *opt;
} data_for_t;

static void worker_for(void *_data, long i, int tid)
{
  data_for_t *data = (data_for_t*)_data;
  aln_barcode(data->opt, &data->seqs[i]);
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
  int i;
  opt_t *opt = (opt_t *)shared;
  if (step == 0) {
    data_for_t *ret;
    ret = calloc(1, sizeof(data_for_t));
    ret->seqs = bseq_read(opt->ks, opt->chunk_size, &ret->n_seqs);
    ret->opt = opt;
    if (ret->seqs) return ret;
    else free(ret);
  } else if (step == 1) {
    data_for_t *data = (data_for_t*)_data;
    kt_for(opt->n_threads, worker_for, data, data->n_seqs);
    return data;
  } else if (step == 2) {
    data_for_t *data = (data_for_t*)_data;
    for (i = 0; i < data->n_seqs; ++i) {
      bseq1_t *s = &data->seqs[i];
      bc_t *bc = opt->bc;
      if (opt->log && (s->l_seq > 450)) {
        fprintf(opt->log, "%s\t%d\t", s->name, s->l_seq);
        log_output(opt, &(s->aln1[0]), 0); log_output(opt, &(s->aln1[1]), 0);
        log_output(opt, &(s->aln2[0]), 0); log_output(opt, &(s->aln2[1]), 0);
        log_output(opt, &(s->aln3[0]), 0); log_output(opt, &(s->aln3[1]), 1);
      }
      sprintf(bc->buffer[s->idx] + bc->offset[s->idx], "@%s\n%s\n+\n%s\n", \
      s->name, s->seq, s->qual);
      bc->offset[s->idx] += (s->l_name + s->l_seq + s->l_qual + 6);
      if (bc->offset[s->idx] > WRITESIZE) {
        bc->buffer[s->idx][bc->offset[s->idx]] = '\0';
        fprintf(bc->ptr[s->idx], "%s", bc->buffer[s->idx]);
        bc->offset[s->idx] = 0U;
      }
      free(s->name); free(s->seq);
      if (s->qual)  free(s->qual);
    }
    free(data->seqs); free(data);
  }
  return 0;
}

void demultiplex_data(opt_t *opt)
{
  kt_pipeline(opt->n_threads, worker_pipeline, opt, 3);
}
