#include "ssw.h"
#include "demultiplex.h"

typedef struct 
{
  char *name, *seq, *qual, *comment;
  int l_name, l_seq, l_qual, l_comment;
  int idx1, strand1, score1;
  int idx2, strand2, score2;
  int idx;
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
    if (ks->comment.s) s->comment = strdup(ks->comment.s), s->l_comment = ks->comment.l;
    s->score1 = s->score2 = 0;
    size += seqs[n++].l_seq;
    if (size >= chunk_size) break;
  }
  *n_ = n;
  return seqs;
}

void aln_core(char *seq, opt_t *opt, int *idx, int *strand, int *score)
{
  int i, num;
  bc_t *bc = opt->bc;

  s_profile *profile;
  s_align *result;

  num = bc->idx;
  int8_t *nt = seq_to_nt(seq, opt->LEN, 0);
  profile = ssw_init(nt, opt->LEN, opt->mat, 5, 2);

  for(i = 0; i < num; i++) {
    result = ssw_align(profile, bc->nt[i], bc->len[i], opt->open, opt->ext, 0, 0, 0, opt->MASK_LEN);
    if (result->score1 > *score)  *idx = i, *strand = 0, *score = result->score1;
    align_destroy(result);
  }
  for(i = 0; i < num; i++) {
    result = ssw_align(profile, bc->nt_rc[i], bc->len[i], opt->open, opt->ext, 0, 0, 0, opt->MASK_LEN);
    if (result->score1 > *score)  *idx = i, *strand = 1, *score = result->score1;
    align_destroy(result);
  }

  init_destroy(profile);
  free(nt);
}

void aln_barcode(opt_t *opt, bseq1_t *seq)
{
  seq->idx = opt->bc->file_num;
  if (seq->l_seq <= opt->LEN) return;

  char *front = (char *)calloc(opt->LEN + 1, sizeof(char));
  char *rear = (char *)calloc(opt->LEN + 1, sizeof(char));
  memcpy(front, seq->seq, opt->LEN);
  memcpy(rear, seq->seq + seq->l_seq - opt->LEN, opt->LEN);
  aln_core(front, opt, &seq->idx1, &seq->strand1, &seq->score1);
  aln_core(rear, opt, &seq->idx2, &seq->strand2, &seq->score2);
  free(front); free(rear);

  if (opt->dual) {
    if (seq->score1 > opt->score && seq->score2 > opt->score) {
      bc_t *bc = opt->bc;
      khint_t k; char tmp[1024];
      if (seq->strand1 == 0 && seq->strand2 == 0) {
        sprintf(tmp, "%s%s", bc->name[seq->idx1], bc->name[seq->idx2]);
        tmp[bc->len[seq->idx1] + bc->len[seq->idx2]] = '\0';
        k = kh_get(dual, bc->h, tmp);
        if (k != kh_end(bc->h)) seq->idx = kh_val(bc->h, k);
      } else if (seq->strand1 == 1 && seq->strand2 == 1) {
        sprintf(tmp, "%s%s", bc->name[seq->idx2], bc->name[seq->idx1]);
        tmp[bc->len[seq->idx2] + bc->len[seq->idx1]] = '\0';
        k = kh_get(dual, bc->h, tmp);
        if (k != kh_end(bc->h)) seq->idx = kh_val(bc->h, k);
      }
    }
  } else {
    if ((seq->score1 > seq->score2) && (seq->score1 > opt->score)) seq->idx = seq->idx1;
    else if ((seq->score1 < seq->score2) && (seq->score2 > opt->score)) seq->idx = seq->idx2;
    else if (((seq->score1 == seq->score2) && (seq->score1 > opt->score)) && (seq->idx1 == seq->idx2)) seq->idx = seq->idx1;
  }
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

      if (opt->log && (s->l_seq > opt->LEN)) {
        fprintf(opt->log, "%s\t%s\t%c\t%d\t%s\t%c\t%d\n", \
        s->name, bc->name[s->idx1], "+-"[s->strand1], s->score1, \
        bc->name[s->idx2], "+-"[s->strand2], s->score2);
      }

      if (opt->mode) {
        sprintf(bc->buffer[s->idx] + bc->offset[s->idx], ">%s\n%s\n", s->name, s->seq);
        bc->offset[s->idx] += (s->l_name + s->l_seq + 3);
      } else {
        sprintf(bc->buffer[s->idx] + bc->offset[s->idx], "@%s %s\n%s\n+\n%s\n", \
        s->name, s->comment, s->seq, s->qual);
        bc->offset[s->idx] += (s->l_name + s->l_seq + s->l_qual + s->l_comment + 7);
      }
      if (bc->offset[s->idx] > WRITESIZE) {
        bc->buffer[s->idx][bc->offset[s->idx]] = '\0';
        fprintf(bc->ptr[s->idx], "%s", bc->buffer[s->idx]);
        bc->offset[s->idx] = 0U;
      }
      free(s->name); free(s->seq);
      if (s->comment) free(s->comment);
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
