# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**nanoplexer** is a C tool for demultiplexing Nanopore long-read sequencing data. It extracts front and rear 150bp sequences from reads, aligns them against barcode sequences using Smith-Waterman alignment, and assigns each read to a sample barcode.

## Build Commands

```sh
make        # Build the nanoplexer binary
make clean  # Remove build artifacts
```

Dependencies: zlib (`-lz`), math library (`-lm`), and pthread (`-pthread`).

## Running

```sh
./nanoplexer -b barcode.fa -p output_path -t 8 input.fastq
./nanoplexer -b barcode.fa -p output_path -l log input.fastq   # with alignment log
cat sequence_id*.fastq | ./nanoplexer -b barcode.fa -p output_path -  # stdin
./nanoplexer -b barcode.fa -d dual_barcode_pair.txt -p output_path input.fastq  # dual barcodes
```

## Architecture

- **main.c**: Entry point, CLI argument parsing via getopt_long
- **option.c**: Option initialization (`opt_init`), validation (`opt_set`), cleanup (`opt_free`)
- **barcode.c**: Barcode data structures and file handling
- **util.c**: Utility functions including reverse complement (`rc_seq`) and sequence conversion (`seq_to_nt`)
- **ssw.c/ssw.h**: Smith-Waterman alignment implementation (from SSW library)
- **kthread.c**: Threading utilities from klib
- **demultiplex.h**: Core data structures:
  - `bc_t`: Barcode information (names, sequences, file pointers, buffers)
  - `opt_t`: Program options (barcode file, output path, threads, alignment scores)
  - `demultiplex_data()`: Main demultiplexing logic
- **kseq.h**: FASTQ/FASTA parsing (from klib)
- **khash.h**: Hash table implementation (from klib)

The demultiplexing process reads FASTQ records, extracts front/rear barcode regions, performs SSW alignment against all barcodes, and writes each read to the appropriate output file based on the best match.

## Key Parameters

- `-t`: Number of threads (default: 3)
- `-L`: Target length for barcode detection (default: 150bp)
- `-m`/`-x`: Match/mismatch scores (default: 2)
- `-o`/`-e`: Gap open/extension penalties (default: 3/1)
- `-s`: Minimum alignment score threshold
- `-d`: Dual barcode pair file (tab-delimited: output_name, 5'_barcode, 3'_barcode)
