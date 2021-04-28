# nanoplexer

#### Introduction
***nanoplexer*** is a standard tool to demultiplex Nanopore long read sequencing data. It extracts front and rear 150bp sequences to align aginst barcode sequences and identify the best hit. To install nanoplexer,
```sh
git clone https://github.com/hanyue36/nanoplexer.git
cd nanoplexer && make
```
The only library dependency is zlib.
#### Usage examples
- demultiplex pacbio data according to barcode file
```
./nanoplexer -b barcode.fa -s sample.txt -p output_path -E 0.08 input.fastq
```
- demultiplex nanopore data and output alignment information
```
./nanoplexer -b barcode.fa -s sample.txt -p output_path -E 0.15 -l log input.fastq
```
- demultiplex pacbio data from stdin stream
```
cat sequence_id*.fastq | ./nanoplexer -b barcode.fa -s sample.txt -p output_path -E 0.08 -
```
#### FAQ
- Format for barcode combination file
<br/>Tab-delimited for each line: output file name, 5' barcode name, 3' barcode name

#### Getting Help
Please use the [Issues page][issue] if you have questions. You may also directly contact Yue Han at hanyue89tj@gmail.com.

[issue]: https://github.com/hanyue36/nanoplexer/issues
