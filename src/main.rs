fn main() {
    let args = star_rs::cli::cli_args();
    match star_rs::cli::run_cli(&args) {
        Ok(result) => {
            star_rs::cli::print_result(&result);
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
