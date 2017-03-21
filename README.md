![pizzly](logo.jpg)

# pizzly 

Fast fusion detection using kallisto

### About

pizzly is a program for detecting gene fusions from RNA-Seq data of cancer samples. 
It requires running [kallisto](https://pachterlab.github.io/kallisto) with the `--fusion`
parameter (available in version `0.43.1` or later). 

### Building

pizzly uses the [SeqAn](http://www.seqan.de/) and requires a recent compiler with `c++14` enabled.
This means `gcc 4.9` or later on linux or `clang 3.5` or later. Since these compilers are not always
awailable on cluster systems we provide a precompiled binary under [releases](https://github.com/pmelsted/pizzly/releases). Compiling pizzly also
requires `cmake 3.0.0` or later.

In order to compile just run the following in the source directory
```
mkdir build
cd build
cmake ..
make
make install
```

### Ingredients

Pizzly requires the reference transcriptome in FASTA format as well as a GTF file describing the transcriptome.
We recommend using the [Ensembl](http://www.ensembl.org/index.html) transcriptomes.

The example below assumes you have your transcriptome in FASTA format as `index.fa.gz`, the GTF file `index.gtf.gz` 
and your paired-end RNA-Seq data sets in `r1.fastq.g` and `r2.fastq.gz`

### Running

First we create the kallisto index

```
kallisto index -i index.idx -k 31 index.fa.gz
```

Next we quantify using kallisto with fusion detection enabled

```
kallisto quant -i index.idx --fusion -o output r1.fastq.gz r2.fastq.gz
```

This creates the file `output/fusion.txt` which is used by pizzly, finally we run pizzly

```
pizzly -k 31 --gtf index.gtf --cache index.cache.txt --align-score 2 \
        --insert-size 250 --fasta index.fa.gz --output test output/fusion.txt
```

The parameters to set are 

* `--insert-size`, which should be the maximum insert size of the paired-end library (kallisto will estimate this from the reads)
* `--align-score`, the number of mismatches allowed when aligning reads to a reference transcript


### Output

The `--output test` parameter is used as a prefix and two files are created `test.fusions.fasta` and `test.json`


### License

GPL v3.0
