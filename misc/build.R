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



fmesher_clang_tidy <- function(files = NULL) {
  if (is.null(files)) {
    files <-
      c(
        "basis.cc",
        "fmesher_helpers.cc",
        "ioutils.cc",
        "locator.cc",
        "mesh.cc",
        "meshc.cc",
        "Rcpp_interface.cc",
        "trees.cc",
        "vector.cc"
      )
  }
  CPPFLAGS <- paste0("-I", R.home("include"),
                     " -I", file.path(system.file(package = "Rcpp"), "include"),
                     " -I/usr/local/include",
                     " -DNDEBUG -DFMESHER_WITH_R")
  SOURCE <- paste0(file.path("src", files), collapse = " ")
  cmd <- paste("clang-tidy", SOURCE, "--", CPPFLAGS)
  print(cmd)
  system(cmd, intern = TRUE)
}
