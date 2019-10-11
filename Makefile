CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function

PROGS=nanoplexer

all:$(PROGS)

nanoplexer: main.c option.c barcode.c util.c ssw.c kthread.c demultiplex.h kseq.h khash.h ssw.h
	$(CC) $(CFLAGS) -o $@ $^ -lz -lm -pthread

clean:
	rm -rf *.o a.out $(PROGS)