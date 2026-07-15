use clap::Parser;
use xenofilters::{config::Cli, run, Error};

fn main() -> Result<(), Error> {
    let cli = Cli::parse();
    run(cli)
}
