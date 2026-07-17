use std::{fs::File, io::{BufReader, BufWriter}};

use flexmap::keys::KHashEntry;

use crate::flexalign::time;


#[derive(Clone, Savefile)]
struct Container {
    pub a: u32,
    pub b: u16,
}


fn test() {


    let mut test_a = vec![0u64; 100_000_000];
    let mut test_b = vec![KHashEntry { key: 0, range_len: 0, range_start: 0 }; 50_000_000];

    for i in 0..50_000_000 {
        test_a[i] = i as u64;
        test_b[i].range_len = i as u32;
        // test_b[i].b = i as u16 % 200;
    }

    let mut file_a = BufWriter::new(File::create("data/large_data/fasta/test_a.bin").unwrap());
    let mut file_b = BufWriter::new(File::create("data/large_data/fasta/test_b.bin").unwrap());

    let config = bincode::config::standard();
    let _ = bincode::encode_into_std_write(&test_a, &mut file_a, config);
    let _ = bincode::encode_into_std_write(&test_b, &mut file_b, config);


    let mut file_a = BufReader::new(File::open("data/large_data/fasta/test_a.bin").unwrap());
    let mut file_b = BufReader::new(File::open("data/large_data/fasta/test_b.bin").unwrap());

    let config = bincode::config::standard();
    let (duration, result) = time(|| bincode::decode_from_std_read::<Vec<KHashEntry>, _, _>(&mut file_b, config));
    eprintln!("Load B took {:?} : {}", duration, result.unwrap().len());

    let config = bincode::config::standard();
    let (duration, result) = time(|| bincode::decode_from_std_read::<Vec<u64>, _, _>(&mut file_a, config));
    eprintln!("Load A took {:?} : {}", duration, result.unwrap().len());

    let mut file_a = File::open("data/large_data/fasta/test_a.bin").unwrap();
    let mut file_b = File::open("data/large_data/fasta/test_b.bin").unwrap();

    let config = bincode::config::standard();
    let (duration, result) = time(|| bincode::decode_from_std_read::<Vec<u64>, _, _>(&mut file_a, config));
    eprintln!("Load A took {:?} : {}", duration, result.unwrap().len());

    let config = bincode::config::standard();
    let (duration, result) = time(|| bincode::decode_from_std_read::<Vec<KHashEntry>, _, _>(&mut file_b, config));
    eprintln!("Load B took {:?} : {}", duration, result.unwrap().len());

}


pub trait TestTrait {
    fn test();
}
pub struct TestStruct;

impl TestTrait for TestStruct {
    fn test() {
        eprintln!("Test successful");
    }
}

pub fn test2() {
    test2_worker::<TestStruct>();
}

pub fn test2_worker<T: TestTrait>() {
    T::test();
}
