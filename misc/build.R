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



fmesher_clang_tidy <- function(files = NULL, standalone_files = NULL, standalone = FALSE) {
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
  if (is.null(standalone_files)) {
    standalone_files <-
      c(
        "fmesher.cc",
        "x11utils.cc"
      )
  }
  if (standalone) {
    CPPFLAGS <- paste0("-std=c++17",
                       " -I", R.home("include"),
                       " -I", file.path(system.file(package = "Rcpp"), "include"),
                       " -I/usr/local/include",
                       " -Imisc/src_standalone",
                       " -DNDEBUG")
    SOURCE <-
      paste0(c(
        file.path("src", files),
        file.path("misc", "src_standalone", standalone_files)
      ), collapse = " ")
  } else {
    CPPFLAGS <- paste0("-std=c++17",
                       " -I", R.home("include"),
                       " -I", file.path(system.file(package = "Rcpp"), "include"),
                       " -I/usr/local/include",
                       " -DNDEBUG -DFMESHER_WITH_R")
    SOURCE <- paste0(file.path("src", files), collapse = " ")
  }
  cmd <- paste("clang-tidy", SOURCE, "--", CPPFLAGS)
  print(cmd)
  system(cmd, intern = TRUE)
}
