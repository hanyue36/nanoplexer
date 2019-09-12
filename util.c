#include <stdarg.h>
#include <time.h>
#include "demultiplex.h"

void print_log(char *format, ...)
{
  time_t t;
  struct tm *lt;
  time(&t);
  lt = localtime(&t);

  va_list ap;
  va_start(ap, format);
  fprintf(stderr, "%02d-%02d-%02d %02d:%02d:%02d ", \
  lt->tm_year+1900, lt->tm_mon+1, lt->tm_mday, lt->tm_hour, \
  lt->tm_min, lt->tm_sec);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr, "\n");
}

int8_t nt_table[128] = {
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void rc_seq(char *seq, int len)
{
  int i;
  char c1, c2;
  for(i = 0; i < len >> 1; ++i) {
    c1 = "ACGT"[3 - nt_table[(int8_t)seq[i]]];
    c2 = "ACGT"[3 - nt_table[(int8_t)seq[len - 1 - i]]];
    seq[i] = c2; seq[len - 1 - i] = c1;
  }
  if (len & 1)
    seq[len >> 1] = "ACGT"[3 - nt_table[(int8_t)seq[len >> 1]]];
}

int8_t *seq_to_nt(char *seq, int len, int flag)
{
  if (flag) rc_seq(seq, len);

  int i;
  int8_t *ret = (int8_t *)malloc(sizeof(int8_t) * len);
  for(i = 0; i < len; i++)  ret[i] = nt_table[(int8_t)seq[i]];

  return ret;
}