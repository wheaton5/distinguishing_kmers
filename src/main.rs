extern crate clap;
extern crate flate2;
extern crate bloom;
extern crate debruijn;

use flate2::read::GzDecoder;
use std::io;
use std::io::prelude::*;
use std::fs::File;

use std::cmp::min;

use bloom::{ASMS,CountingBloomFilter};

use clap::{Arg, App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;

fn main() {
    let matches = App::new("distinguishing_kmers")
        .author("Haynes Heaton <whheaton@gmail.com>")
        .about("Thresholded kmer set subtraction using counting bloom filters. Designed for finding haplotype distinguishing kmers with pedigree illumina data.")
        .arg(Arg::with_name("kmers_in")
                 .short("i")
                 .long("kmers_in")
                 .takes_value(true)
                .multiple(true)
                .required(true)
                 .help("fastq files from which kmers are taken"))
        .arg(Arg::with_name("kmers_subtract")
                 .short("s")
                 .long("kmers_subtract")
                 .takes_value(true)
                .required(true)
                .multiple(true)
                 .help("fastq files of kmers to be threshold subtracted from kmers_in files."))
        .arg(Arg::with_name("min_source_count")
                .long("min_source_count")
                .takes_value(true)
                .help("min kmer count in input files to report in output")
                .required(true))
        .arg(Arg::with_name("min_subtract_count")
                .long("min_subtract_count")
                .takes_value(true)
                .help("min kmer count in subtract files to subtract from kmers_in")
                .required(true))
        .arg(Arg::with_name("difference_threshold")
                .long("difference_threshold")
                .takes_value(true)
                .help("threshold on difference between counts of infile and subtract files"))
        .arg(Arg::with_name("estimated_kmers")
                .long("estimated_kmers")
                .help("estimated number of unique kmers. Good rule of thumb is genome size * 2")
                .takes_value(true))
        .get_matches();

    eprintln!("");
    let min_source_count = matches.value_of("min_source_count").unwrap_or("5");
    let min_source_count: u32 = min_source_count.to_string().parse::<u32>().unwrap();
    let min_subtract_count = matches.value_of("min_subtract_count").unwrap_or("3");
    let min_subtract_count: u32 = min_subtract_count.to_string().parse::<u32>().unwrap();
    let difference_threshold = matches.value_of("difference_threshold").unwrap_or("0");
    let difference_threshold: u32 = difference_threshold.to_string().parse::<u32>().unwrap();
    let estimated_kmers = matches.value_of("estimated_kmers").unwrap_or("6000000000");
    let estimated_kmers: u32 = estimated_kmers.to_string().parse::<u32>().unwrap();
    let kmers_in: Vec<_> = matches.values_of("kmers_in").unwrap().collect();
    let kmers_subtract: Vec<_> = matches.values_of("kmers_subtract").unwrap().collect();
    eprintln!("This program will take kmers from files");
    for f in &kmers_in {
        eprintln!("\t{}",f);
    }
    eprintln!("with at least {} count and subtract kmers from files",min_source_count);
    for f in &kmers_subtract {
        eprintln!("\t{}",f);
    }
    eprintln!("with at least {} count",min_subtract_count);
    eprintln!("Reminder! This program uses counting bloom filters and thus counts are approximate but can only be correct or over-counts");
    let k_size: usize = 21; // THIS CANNOT BE CHANGED WITHOUT CHANGING Kmer datatype
    let source_counting_bits = 6;
    let subtract_counting_bits = 4;
    let bloom_kmer_in_counts = count_kmers_fastq(&kmers_in, source_counting_bits, estimated_kmers, k_size);
    let bloom_kmer_subtract_counts = count_kmers_fastq(&kmers_subtract, subtract_counting_bits, estimated_kmers, k_size);

    subtract_kmers(kmers_in, bloom_kmer_in_counts, bloom_kmer_subtract_counts, min_source_count, min_subtract_count, difference_threshold, k_size);
}


fn subtract_kmers(kmers_in: Vec<&str>, kmer_counts: CountingBloomFilter, subtract_counts: CountingBloomFilter, 
                min_source_count: u32, min_subtract_count: u32, difference_threshold: u32, k_size: usize) {
    for kmer_file in kmers_in {
        let file = match File::open(kmer_file) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening file {} with error {}", kmer_file,error)
            },
        };
        let gz = GzDecoder::new(file);
        for (line_number, line) in io::BufReader::new(gz).lines().enumerate() {
            if line_number % 4 == 1 {
                let dna: DnaString = DnaString::from_dna_string(&line.unwrap());
                for kmer_start in 0..(dna.len()-k_size) {
                    let k: Kmer21 = dna.get_kmer(kmer_start);
                    let krc = k.rc();
                    let to_hash = min(k.to_u64(), krc.to_u64());
                    let source_count = kmer_counts.estimate_count(&to_hash);
                    if source_count >= min_source_count {
                        let subtract_count = subtract_counts.estimate_count(&to_hash);
                        if subtract_count <= min_subtract_count && source_count - subtract_count >= difference_threshold {
                            if k.to_u64() < krc.to_u64() {
                                println!("{}\t{}\t{}",k.to_string(),source_count,subtract_count);
                            } else {
                                println!("{}\t{}\t{}",krc.to_string(),source_count,subtract_count);
                            }
                        }
                    }
                }
            }
        }
    }
}

fn count_kmers_fastq(kmers_in: &Vec<&str>, counting_bits: usize, estimated_kmers: u32, k_size: usize) -> CountingBloomFilter {
    let mut kmer_counts: CountingBloomFilter = CountingBloomFilter::with_rate(counting_bits, 0.05, estimated_kmers);
    for kmer_file in kmers_in {
        let file = match File::open(kmer_file) {
            Ok(file) => file,
            Err(error) => {
                panic!("There was a problem opening the file {} with error {}", kmer_file, error)
            },
        };
        let gz = GzDecoder::new(file);
        for (line_number, line) in io::BufReader::new(gz).lines().enumerate() {
            if line_number % 4 == 1 {
                let dna: DnaString = DnaString::from_dna_string(&line.unwrap()); //Vec<u64> 2bit encoded
                for kmer_start in 0..(dna.len() - k_size) {
                    let k: Kmer21 = dna.get_kmer(kmer_start); //statically typed kmer size
                    kmer_counts.insert(&min(k.to_u64(), k.rc().to_u64()));
                }
            }
        }
    }
    kmer_counts
}

type Kmer21 = VarIntKmer<u64, K21>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct K21;

impl KmerSize for K21 {
    fn K() -> usize {
        21
    }
}
