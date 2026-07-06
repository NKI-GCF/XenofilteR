use xenofilters::{Error,config::Config, run};
use clap::Parser;

fn main() -> Result<(), Error> {
    let config = Config::parse();
    run(config)
}
