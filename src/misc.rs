
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

    let _ = save(&mut file_a, 1, &test_a);
    let _ = save(&mut file_b, 1, &test_b);


    let mut file_a = BufReader::new(File::open("data/large_data/fasta/test_a.bin").unwrap());
    let mut file_b = BufReader::new(File::open("data/large_data/fasta/test_b.bin").unwrap());

    let (duration, result) = time(|| load::<Vec<KHashEntry>>(&mut file_b, 1));
    eprintln!("Load B took {:?} : {}", duration, result.unwrap().len());

    let (duration, result) = time(|| load::<Vec<u64>>(&mut file_a, 1));
    eprintln!("Load A took {:?} : {}", duration, result.unwrap().len());

    let mut file_a = File::open("data/large_data/fasta/test_a.bin").unwrap();
    let mut file_b = File::open("data/large_data/fasta/test_b.bin").unwrap();

    let (duration, result) = time(|| load::<Vec<u64>>(&mut file_a, 1));
    eprintln!("Load A took {:?} : {}", duration, result.unwrap().len());

    let (duration, result) = time(|| load::<Vec<KHashEntry>>(&mut file_b, 1));
    eprintln!("Load B took {:?} : {}", duration, result.unwrap().len());

}