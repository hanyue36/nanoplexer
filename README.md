# nanoplexer

#### Introduction
***nanoplexer*** is a standard tool to demultiplex Nanopore long read sequencing data. It extracts front and rear 150bp sequences to align aginst barcode sequences and identify the best hit. To install nanoplexer,
```sh
git clone https://github.com/hanyue36/nanoplexer.git
cd nanoplexer && make
```
The only library dependency is zlib


##### Install with conda
![conda_badge](
https://anaconda.org/bioconda/nanoplexer/badges/version.svg)
![conda_badge](https://anaconda.org/bioconda/nanoplexer/badges/latest_release_date.svg
)
![conda_badge](https://anaconda.org/bioconda/nanoplexer/badges/platforms.svg
)
![conda_badge](https://anaconda.org/bioconda/nanoplexer/badges/downloads.svg)
```sh
conda install -c bioconda nanoplexer
```

#### Usage examples
- demultiplex data according to barcode file
```
./nanoplexer -b barcode.fa -p output_path -t 8 input.fastq
```
- demultiplex data and output alignment information
```
./nanoplexer -b barcode.fa -p output_path -l log input.fastq
```
- demultiplex data from stdin stream
```
cat sequence_id*.fastq | ./nanoplexer -b barcode.fa -p output_path -
```
- demultiplex data according to dual barcode file
```
./nanoplexer -b barcode.fa -d dual_barcode_pair.txt -p output_path input.fastq
```
#### FAQ
- Format for dual barcode pair file
<br/>Tab-delimited for each line: output file name, 5' barcode name, 3' barcode name

#### Getting Help
Please use the [Issues page][issue] if you have questions. You may also directly contact Yue Han at hanyue89tj@gmail.com.

[issue]: https://github.com/hanyue36/nanoplexer/issues
