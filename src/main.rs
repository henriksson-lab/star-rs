mod cli;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[used]
#[cfg_attr(target_os = "linux", unsafe(link_section = ".init_array"))]
static MIMALLOC_OPTIONS_INIT: extern "C" fn() = init_mimalloc_options;

extern "C" fn init_mimalloc_options() {
    unsafe {
        // mimalloc v3 option ids: eager commit, arena eager commit, purge delay.
        mi_option_set_default(3, 0);
        mi_option_set_default(4, 0);
        mi_option_set_default(15, 0);
    }
}

unsafe extern "C" {
    fn mi_option_set_default(option: i32, value: libc::c_long);
}

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
