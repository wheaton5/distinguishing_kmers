#[macro_use]
extern crate clap;
extern crate flate2;
extern crate bloom;
extern crate debruijn;
extern crate dna_io;

use flate2::read::GzDecoder;
use std::io;
use std::io::prelude::*;
use std::fs::File;

use std::cmp::min;

use bloom::{ASMS,CountingBloomFilter,BloomFilter};
use std::collections::HashMap;

use clap::{App};

use debruijn::*;
use debruijn::dna_string::*;
use debruijn::kmer::*;
static mut KMER_SIZE: usize = 21;

fn main() {
    let (min_source_count, 
        max_subtract_count, 
        difference_threshold,
        modimizer,
        mod_index,
        kmers_in,
        kmers_subtract) = load_params();
    eprintln!("");
    eprintln!("This program will take kmers from files");
    for f in &kmers_in {
        eprintln!("\t{}",f);
    }
    eprintln!("with at least {} count and subtract kmers from files",min_source_count);
    for f in &kmers_subtract {
        eprintln!("\t{}",f);
    }
    eprintln!("with at least {} count",max_subtract_count);
    unsafe {
        eprintln!("Only considers kmers where the 2bit representation of the kmer % 9 == 0 to improve speed and memory usage. This is reasonable because otherwise for every difference you get {} overlapping {}mers but with this you get on average just over 2 overlapping {}mers.",KMER_SIZE,KMER_SIZE,KMER_SIZE);
    }
    //let source_counting_bits = 6;
    //let subtract_counting_bits = 4;
    let kmer_in_counts = count_kmers_fastq_exact(&kmers_in, modimizer, mod_index);
    let kmer_subtract_counts = count_kmers_fastq_exact(&kmers_subtract, modimizer, mod_index);
    subtract_kmers_exact(kmer_in_counts, kmer_subtract_counts, min_source_count, max_subtract_count, difference_threshold);
}


fn subtract_kmers_exact(kmers_in: HashMap<u64,u16>, subtract_counts: HashMap<u64,u16>, 
                min_source_count: u16, min_subtract_count: u16, difference_threshold: u16) {
    for (kmeru64, source_count) in &kmers_in {
        if *source_count >= min_source_count {
            let subtract_count = match subtract_counts.get(&kmeru64) {Some(x) => *x, None => 0,};
            if subtract_count <= min_subtract_count && source_count - subtract_count >= difference_threshold {
                let kmer = KmerX::from_u64(*kmeru64);
                println!("{}\t{}\t{}",kmer.to_string(), source_count, subtract_count);
            }
        }
    }
}

fn count_kmers_fastq_exact(kmers_in: &Vec<String>, modimizer: i64, mod_index: i64) -> HashMap<u64,u16> {
    let mut kmer_counts: HashMap<u64, u16> = HashMap::new();

    for kmer_file in kmers_in {
        let mut reader = dna_io::DnaReader::from_path(kmer_file);

        'recloop: loop {
            let record = match reader.next() {
                None => break 'recloop,
                Some(x) => x,
            };
                
            for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                let to_hash = min(k.to_u64(), k.rc().to_u64());
                if to_hash % modimizer as u64 == mod_index as u64 {
                    let mut count = kmer_counts.entry(to_hash).or_insert(0);
                    if *count < 65535 {
                        *count += 1;
                    }
                }
            }
        }
    }
    kmer_counts
}

#[allow(dead_code)]
fn subtract_kmers(kmers_in: Vec<&str>, kmer_counts: CountingBloomFilter, subtract_counts: CountingBloomFilter, 
                min_source_count: u32, min_subtract_count: u32, difference_threshold: u32, k_size: usize, estimated_kmers: u32) {
    let mut visited_kmer = BloomFilter::with_rate(0.03, estimated_kmers);
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
                for kmer_start in 0..(dna.len() - k_size + 1) {
                    let k: KmerX = dna.get_kmer(kmer_start);
                    let krc = k.rc();
                    let to_hash = min(k.to_u64(), krc.to_u64());
                    if visited_kmer.contains(&to_hash) {
                        continue;
                    }
                    visited_kmer.insert(&to_hash);
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

#[allow(dead_code)]
fn count_kmers_fastq(kmers_in: &Vec<String>, counting_bits: usize, estimated_kmers: u32, k_size: usize, 
        _modimizer: i64, _mod_index: i64) -> CountingBloomFilter {
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
                for kmer_start in 0..(dna.len() - k_size + 1) {
                    let k: KmerX = dna.get_kmer(kmer_start); //statically typed kmer size
                    println!("{}",k.to_string());
                    kmer_counts.insert(&min(k.to_u64(), k.rc().to_u64()));
                }
            }
        }
    }
    kmer_counts
}

type KmerX = VarIntKmer<u64, KX>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct KX;

impl KmerSize for KX {
    fn K() -> usize {
        unsafe {
            KMER_SIZE
        }
    }
}

fn load_params() -> (u16, u16, u16, i64, i64, Vec<String>, Vec<String>) {
    let yaml_params = load_yaml!("params.yml");
    let params = App::from_yaml(yaml_params).get_matches();
    let kmer_size = params.value_of("kmer_size").unwrap_or("21");
    let kmer_size: usize = kmer_size.to_string().parse::<usize>().unwrap();
    if kmer_size > 32 { panic!("kmer size too large, only support k < 32. found {}", kmer_size); }
    unsafe {
        KMER_SIZE = kmer_size;
    }
    let min_source_count = params.value_of("min_source_count").unwrap_or("5");
    let min_source_count: u16 = min_source_count.to_string().parse::<u16>().unwrap();
    let max_subtract_count = params.value_of("max_subtract_count").unwrap_or("0");
    let max_subtract_count: u16 = max_subtract_count.to_string().parse::<u16>().unwrap();
    let difference_threshold = params.value_of("difference_threshold").unwrap_or("0");
    let difference_threshold: u16 = difference_threshold.to_string().parse::<u16>().unwrap();
    let modimizer = params.value_of("modimizer").unwrap_or("9");
    let modimizer: i64 = modimizer.to_string().parse::<i64>().unwrap();
    let mod_index = params.value_of("mod_index").unwrap_or("0");
    let mod_index: i64 = mod_index.to_string().parse::<i64>().unwrap();
    let mut kmers_in: Vec<String> = Vec::new();
    for kmer_in in params.values_of("kmers_in").unwrap() {
        kmers_in.push(kmer_in.to_string());
    }
    let mut kmers_subtract: Vec<String> = Vec::new();
    for kmer_subtract in params.values_of("kmers_subtract").unwrap() {
        kmers_subtract.push(kmer_subtract.to_string());
    }
    (min_source_count,
    max_subtract_count,
    difference_threshold,
    modimizer,
    mod_index,
    kmers_in,
    kmers_subtract)
}
