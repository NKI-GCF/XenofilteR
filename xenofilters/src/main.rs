use clap::Parser;
use xenofilters::{Error, config::Cli, run};

fn main() -> Result<(), Error> {
    let cli = Cli::parse();
    run(cli)
}
