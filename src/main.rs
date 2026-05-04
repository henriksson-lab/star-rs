mod cli;

pub mod generated {
    pub mod functions;
    pub mod structs;
}

fn main() {
    let args = cli::cli_args();
    match cli::run_cli(&args) {
        Ok(result) => {
            cli::print_result(&result);
            std::process::exit(result.exit_code);
        }
        Err(err) => {
            eprint!("{}", err);
            if !err.ends_with('\n') {
                eprintln!();
            }
            std::process::exit(1);
        }
    }
}
