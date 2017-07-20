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

The example below assumes you have your transcriptome in FASTA format as `transcripts.fa.gz`, the GTF file `transcripts.gtf.gz` 
and your paired-end RNA-Seq data sets in `r1.fastq.gz` and `r2.fastq.gz`

### Running

First we create the kallisto index

```
kallisto index -i index.idx -k 31 transcripts.fa.gz
```

Next we quantify using kallisto with fusion detection enabled

```
kallisto quant -i index.idx --fusion -o output r1.fastq.gz r2.fastq.gz
```

This creates the file `output/fusion.txt` which is used by pizzly, finally we run pizzly

```
pizzly -k 31 --gtf transcripts.gtf --cache index.cache.txt --align-score 2 \
        --insert-size 400 --fasta transcripts.fa.gz --output test output/fusion.txt
```

The parameters to set are 

* `--insert-size`, which should be the maximum insert size of the paired-end library (kallisto will estimate this from the reads, default used is 400)
* `--align-score`, the number of mismatches allowed when aligning reads to a reference transcript (default used is 2)
`--ignore-protein`, ignore any information about protein coding in the annotation, **warning** this will probably lead to an increase in the number of false positives reported.
* `--cache`, if this file does not exist, pizzly will parse the GTF (which can take up to a minute or two) and store the required data in the cached file. If the cache file exists (if you've run pizzly previously on the same GTF file), pizzly will parse this file instead, which is much faster than parsing the GTF. 


A more sophisticated example is in the `test` directory which contains a `snakemake` workflow to index, quantify, call fusions and requantify using `kallisto` and `pizzly`.



### Output

The `--output test` parameter is used as a prefix and two files are created `test.fusions.fasta` and `test.json`, this contains the filtered fusion calls. For unfiltered fusion calls, use the corresponding `.unfiltered` files.

### Scripts

The `scripts` subfolder contains useful Python scripts

- `get_fragment_length.py` examines an `abundance.h5` produced by `kallisto` and finds the 95th percentile of the fragment length distribution
- `flatten_json.py` reads the `.json` output and converts to a simple gene table

### Annotations

pizzly has been tested on [Ensembl](http://www.ensembl.org/) (versions 75+) and [Gencode](http://www.gencodegenes.org/) (version 19+) annotations. We recommend using the latest Ensembl annotations (version 87 [GTF](ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/), [FASTA](ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cdna/)) for running with pizzly.

Note that for gencode you will need to modify the FASTA file to remove pipe symbols (`|`) from the target names. The following should work (use `gzcat` on macosx)

```
zcat gencode.v26.transcripts.fa.gz  | tr '|' ' ' | gzip -1 >  gencode.v26.transcripts.fixed.fa.gz
```

The FASTA file used must be the same one that was used to build the kallisto index.


### License

GPL v3.0
