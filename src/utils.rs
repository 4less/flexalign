use std::{cmp::max, collections::HashSet, fs::File, io::{self, BufRead}, path::Path};



fn read_lines_from_file(filename: &str) -> io::Result<Vec<String>> {
    let path = Path::new(filename);
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let lines: Vec<String> = reader
        .lines()
        .collect::<Result<_, _>>()?;

    Ok(lines)
}


pub fn infer_output_prefix(input: &[String]) -> Vec::<String> {
    let tokens: Vec<Vec<_>> = input.iter()
        .map(|s| { 
            let s = s.strip_suffix(".gz").unwrap_or(s);
            let s = s.strip_suffix(".bz").unwrap_or(s);
            let s = s.strip_suffix(".bz2").unwrap_or(s);     // Remove .gz if present
            let s = s.rsplit_once('.').map_or(s, |(left, _)| left); // Split by rightmost '.' and keep left part
            s.split("/").collect::<Vec<&str>>()
        })
        .collect();

    let longest = tokens.iter().fold(0, |acc,vec| return max(acc, vec.len()));

    let duplications = (0..longest)
        .map(|i| {
            tokens.iter()
                .map(|t| t.get(i).unwrap_or(&""))
                .collect::<HashSet<_>>().len()
        }).collect::<Vec<_>>();

    let indices = duplications.iter().enumerate().filter(|(_, &x)| x > 1).map(|(i, _)| i).collect::<Vec<_>>();

    let output_prefixes = tokens.iter().map(|token| {
        indices.iter()
            .filter(|&&i| i < token.len())
            .map(|&i| token[i].to_string())
            .collect::<Vec<String>>().join("_")
    }).collect::<Vec<_>>();

    output_prefixes
}
