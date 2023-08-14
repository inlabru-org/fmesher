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



fmesher_clang_tidy <- function() {
  CPPFLAGS <- "-I/home/flindgre/local/R-4.3.1/lib/R/include -DNDEBUG -DFMESHER_WITH_R -I/home/flindgre/R/x86_64-pc-linux-gnu-library/4.3/Rcpp/include -I/usr/local/include"
  SOURCE <- paste0(
    "src/",
    c(
      "basis.cc",
      "fmesher.cc",
      "fmesher_helpers.cc",
      "ioutils.cc",
      "locator.cc",
      "mesh.cc",
      "Rcpp_interface.cc",
      "trees.cc",
      "vector.cc"
    ),
    collapse = " "
  )
  cmd <- paste("clang-tidy", SOURCE, "--", CPPFLAGS)
  print(cmd)
  system(cmd, intern = TRUE)
}
