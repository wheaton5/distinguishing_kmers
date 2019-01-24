# distinguishing_kmers
Memory efficient (approximate) thresholded kmer set subtraction using a mod hash minimizer.

To Install
Requirements: Rust ver 1.3 or later
```
git clone git@github.com:wheaton5/distinguishing_kmers.git
cd distinguishing_kmers
cargo build --release
```
```
./target/release/distinguishing_kmers -h
distinguishing_kmers
Haynes Heaton <whheaton@gmail.com>
Thresholded kmer set subtraction using counting bloom filters. Designed for finding haplotype distinguishing kmers with
pedigree illumina data.

USAGE:
    distinguishing_kmers [OPTIONS] --kmers_in <kmers_in>... --kmers_subtract <kmers_subtract>... --min_source_count <min_source_count> --min_subtract_count <min_subtract_count>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --difference_threshold <difference_threshold>
            threshold on difference between counts of infile and subtract files

        --estimated_kmers <estimated_kmers>
            estimated number of unique kmers. Good rule of thumb is genome size * 2

    -i, --kmers_in <kmers_in>...                         fastq files from which kmers are taken
    -s, --kmers_subtract <kmers_subtract>...
            fastq files of kmers to be threshold subtracted from kmers_in files.

        --min_source_count <min_source_count>            min kmer count in input files to report in output
        --max_subtract_count <max_subtract_count>        max kmer count in subtract files before subtracting from kmers_in
```        

Example using tiny made up data
```
./target/release/distinguishing_kmers --kmers_in test/data/kmers_in.fastq.gz --kmers_subtract test/data/kmers_subtract.fastq.gz --min_source_count 5 --max_subtract_count 4

This program will take kmers from files
        test/data/kmers_in.fastq.gz
with at least 5 count and subtract kmers from files
        test/data/kmers_subtract.fastq.gz
with at least 3 count
Reminder! This program uses counting bloom filters and thus counts are approximate but can only be correct or over-counts
CTCGACTCTGCCGCAGCTACC   6       0

```
