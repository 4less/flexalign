
#![feature(map_try_insert)]
#![feature(const_option)]
#![feature(const_trait_impl)]
// #![feature(effects)]

use std::fs::File;
use std::io::{BufReader, BufWriter};

use clap::Parser;
#[allow(unused_parens)]
pub mod options;
pub mod utils;
use colored::control::SHOULD_COLORIZE;
use colored::Colorize;
use flexalign::flexalign::{run, time};
use flexalign::misc::test2;
use flexalign::options::Args;
use flexmap::keys::{FMKeysHash, KHashEntry};
use savefile::{load, save};
use savefile_derive::Savefile;



fn logo() -> String {
    format!(
        "{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}", 
        "  __ _           ".yellow(), "     _ _             ".red().bold(),
        " / _| | _____  __".yellow(), "__ _| (_) __ _ _ __  ".red().bold(),
        "| |_| |/ _ \\ \\/ /".yellow(), " _` | | |/ _` | '_ \\ ".red().bold(),
        "|  _| |  __/>  <".yellow(), " (_| | | | (_| | | | |".red().bold(),
        "|_| |_|\\___/_/\\_\\".yellow(), "__,_|_|_|\\__, |_| |_|".red().bold(),
        "                          |___/       ".red().bold(),
    )
}

fn logo2() -> String {
    format!("\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}\n{}\n{}\n{}\n",
        " .d888 888                    ","        888 d8b                   ".red().bold(),
        "d88P\"  888                  ","          888 Y8P                   ".red().bold(),
        "888    888                   ","         888                       ".red().bold(),
        "888888 888  .d88b.  888  888 "," 8888b.  888 888  .d88b.  88888b.  ".red().bold(),
        "888    888 d8P  Y8b `Y8bd8P' ","    \"88b 888 888 d88P\"88b 888 \"88b ".red().bold(),
        "888    888 88888888   X88K   ",".d888888 888 888 888  888 888  888 ".red().bold(),
        "888    888 Y8b.     .d8\"\"8b."," 888  888 888 888 Y88b 888 888  888 ".red().bold(),
        "888    888  \"Y8888  888  888 ","\"Y888888 888 888  \"Y88888 888  888 ".red().bold(),
        "                                                   888          ".red().bold(),
        "                                              Y8b d88P          ".red().bold(),
        "                                               \"Y88P\"   ".red().bold())
}

fn main() {
    // CAUTION: do not colorize anything that goes into stdout
    // otherwise the resulting sam files will be broken.
    SHOULD_COLORIZE.set_override(true);

    eprintln!("{}", logo());

    let args: Args = Args::parse();
    let (duration, _) = time(|| run(args));

    eprintln!("Flexalign took {:?}", duration);
}