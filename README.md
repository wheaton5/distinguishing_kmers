# distinguishing_kmers
Memory efficient (approximate) thresholded kmer set subtraction using counting bloom filters.

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
        --min_subtract_count <min_subtract_count>        min kmer count in subtract files to subtract from kmers_in
```        
