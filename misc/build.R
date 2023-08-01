fmesher_install <- function(repo = ".", debug = FALSE) {
  flags <- pkgbuild::compiler_flags(debug = debug)
  flag_names <- names(flags)
  if (debug) {
    flags <- paste0(flags, " -DFMESHER_DEBUG")
    names(flags) <- flag_names
  }
  flags <- as.list(flags)
  if (identical(repo, ".")) {
    pkgbuild::with_debug(
      devtools::install_local(repo, force = TRUE),
      CFLAGS = flags[["CFLAGS"]],
      CXXFLAGS = flags[["CXXFLAGS"]],
      FFLAGS = flags[["FFLAGS"]],
      FCFLAGS = flags[["FCFLAGS"]],
      debug = debug)
  } else {
    pkgbuild::with_debug(
      devtools::install_github(repo, force = TRUE),
      CFLAGS = flags[["CFLAGS"]],
      CXXFLAGS = flags[["CXXFLAGS"]],
      FFLAGS = flags[["FFLAGS"]],
      FCFLAGS = flags[["FCFLAGS"]],
      debug = debug)
  }
}
