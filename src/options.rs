use std::{path::PathBuf, str::FromStr};

use clap::Parser;
use clap_derive::Args;

use crate::utils::infer_output_prefix;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)] // term_width sets it fixed, max term_width can be smaller
pub struct Args {
    /// Forward read of pair (.fastq, .fq)
    #[arg(num_args(0..), short = '1', long, default_values_t = ["".to_string()], action = clap::ArgAction::Append)]
    pub fwd: Vec<String>,

    /// Reverse read of pair (.fastq, .fq)
    #[arg(num_args(0..), short = '2', long, default_values_t = ["".to_string()], action = clap::ArgAction::Append)]
    pub rev: Vec<String>,

    /// Output file
    #[arg(short = '0', long)] // String::default()
    pub output: Option<String>,

    /// Database reference
    #[arg(short = 'r', long = "reference", default_value_t = String::default())]
    pub reference: String,

    /// Input map file 
    #[arg(short, long, default_value_t = String::default())]
    pub map: String,

    /// threads 
    #[arg(short, long, default_value_t = 1)]
    pub threads: u32,

    /// How many minimizers should be looked at
    #[arg(short = 'a', long = "ranges", default_value_t = 15)]
    pub ranges: u32,

    /// For a single minimizer, how many occurances may there be at max.
    #[arg(short = 'b', long = "max-range-size", default_value_t = 256)]
    pub max_range_size: usize,

    /// For all occurrences of a key, flexalign only takes the seeds with the highest matching flanking region.
    /// This limits the number of values to be retrieved in this scenario. 
    #[arg(short = 'f', long = "max-best-flex", default_value_t = 16)]
    pub max_best_flex: usize,

    /// After the seeds are grouped into anchors, the top x will be extended with the use of hamming distance.
    /// This affects speed negatively but sensitivity and precision positively
    #[arg(short = 'x', long = "extend-top-x", default_value_t = 4)]
    pub extend_top_x: usize,

    /// align the top y anchors. This happens after anchor extension 
    #[arg(short = 'y', long = "align-top-y", default_value_t = 4)]
    pub align_top_y: usize,

    /// Minimum number of ranges for lookup. With max-best-flex defines, none of the ranges might actually yield any seeds.
    #[arg(long = "min-ranges", default_value_t = 4)]
    pub min_ranges: usize,

    /// force_build
    #[arg(long = "force-build", action)]
    pub force_build: bool,

    /// force_build
    #[arg(long = "debug", action)]
    pub debug: bool,
}

#[derive(Debug)]
pub struct Options {
    pub fwd: Vec<PathBuf>,
    pub rev: Vec<Option<PathBuf>>,
    pub output_prefix: Option<Vec<PathBuf>>,
    pub reference: PathBuf,
    pub reference_database: PathBuf,
    
    pub args: Args,
}


impl Options {
    pub fn from_args(args: Args) -> Self {
        let mut options = Options {
            fwd: vec![PathBuf::default(); 0],
            rev: vec![None; 0],
            reference: PathBuf::default(),
            reference_database: PathBuf::default(),
            output_prefix: None,
            args: args,
        };
        Self::init(&mut options);
        options
    }

    pub fn init(&mut self) {
        self.fwd.extend(self.args.fwd.iter().map(|x| x.into()));
        self.rev.extend(self.args.rev.iter().map(|x| Some(x.into())));

        if self.fwd.len() > 1 {
            if self.args.output.is_none() {
                panic!("When processing multiple files in one run you need to provide an output folder for the results to be stored. (--output FOLDER)")
            }

            let inputs = self.fwd.iter().map(|x| x.to_string_lossy().into_owned()).collect::<Vec<String>>();
            self.output_prefix = Some(infer_output_prefix(&inputs)
                .iter()
                .map(|s| { 
                    let mut p = PathBuf::from_str(&self.args.output.as_ref().unwrap()).expect("Cannot turn string into path");
                    p.push(s);
                    p
                })
                .collect::<Vec<_>>());
        } else if self.fwd.len() == 1 && self.args.output.is_some() {
            let s = self.fwd.first().unwrap().to_str().unwrap();
            let s = s.strip_suffix(".gz").unwrap_or(s);
            let s = s.strip_suffix(".bz").unwrap_or(s);
            let s = s.strip_suffix(".bz2").unwrap_or(s);     // Remove .gz if present
            let s = s.rsplit_once('.').map_or(s, |(left, _)| left);
            
            self.output_prefix = Some(vec![PathBuf::from(s); 0]);
        }
        
        if self.output_prefix.is_some() {
            for s in self.output_prefix.as_ref().unwrap() {
                println!("{:?}", s);
            }
        }

        self.reference.push(self.args.reference.clone());
    }
}

pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}