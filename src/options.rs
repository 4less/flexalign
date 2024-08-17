use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)] // term_width sets it fixed, max term_width can be smaller
pub struct Args {
    /// Forward read of pair (.fastq, .fq)
    #[arg(short = '1', long, default_value_t = String::default())]
    pub fwd: String,

    /// Reverse read of pair (.fastq, .fq)
    #[arg(short = '2', long, default_value_t = String::default())]
    pub rev: String,

    /// Output file
    #[arg(short = '0', long, default_value_t = String::default())]
    pub output: String,

    /// Database reference
    #[arg(short = 'r', long = "reference", default_value_t = String::default())]
    pub reference: String,

    /// Input map file 
    #[arg(short, long, default_value_t = String::default())]
    pub map: String,

    /// threads 
    #[arg(short, long, default_value_t = 1)]
    pub threads: u32,

    /// ranges 
    #[arg(short = 'a', long = "ranges", default_value_t = 15)]
    pub ranges: u32,

    /// max_range_size 
    #[arg(short = 'b', long = "max-range-size", default_value_t = 256)]
    pub max_range_size: usize,

    /// max_range_size 
    #[arg(short = 'f', long = "max-best-flex", default_value_t = 16)]
    pub max_best_flex: usize,

    /// Minimum number of ranges for lookup. With max-best-flex defines, none of the ranges might actually yield any seeds.
    #[arg(long = "min-ranges", default_value_t = 4)]
    pub min_ranges: usize,

    /// force_build
    #[arg(long = "force-build", action)]
    pub force_build: bool,
}

#[derive(Debug)]
pub struct Options {
    pub fwd: Vec<PathBuf>,
    pub rev: Vec<Option<PathBuf>>,
    pub reference: PathBuf,
    pub reference_database: PathBuf,
    output_prefix: Vec<PathBuf>,
    
    pub args: Args,
}


impl Options {
    pub fn from_args(args: Args) -> Self {
        let mut options = Options {
            fwd: vec![PathBuf::default(); 0],
            rev: vec![None; 0],
            reference: PathBuf::default(),
            reference_database: PathBuf::default(),
            output_prefix: vec![PathBuf::default(); 0],
            args: args,
        };
        Self::init(&mut options);
        options
    }

    pub fn init(&mut self) {
        self.fwd.push(self.args.fwd.clone().into());

        self.rev.push(if self.args.rev.is_empty() { 
            None 
        } else {
            Some(self.args.rev.clone().into())
        });
        
        self.reference.push(self.args.reference.clone());
        self.output_prefix.push("test.out".into());
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