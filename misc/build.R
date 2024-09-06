fmesher_install <- function(repo = ".", debug = FALSE) {
  # Expanded flag support:
  pkgbuild_with_debug <- function(code, flags, debug = TRUE) {
    defaults <- pkgbuild::compiler_flags(debug = debug)
    flags <- unlist(utils::modifyList(as.list(defaults), as.list(flags)))
    pkgbuild:::withr_with_makevars(flags, code)
  }

  flags <- pkgbuild::compiler_flags(debug = debug)
  flag_names <- names(flags)
  if (debug) {
    flags <- paste0(flags, " -DFMESHER_DEBUG=1")
    names(flags) <- flag_names
  }
  flags <- as.list(flags)
  if (identical(repo, ".")) {
    pkgbuild_with_debug(
      # pak::pkg_install(paste0("local::", repo, "?source&reinstall"), ask = FALSE),
      devtools::install(repo),
      flags = c(
        CFLAGS = flags[["CFLAGS"]],
        CXXFLAGS = flags[["CXXFLAGS"]],
        CXX11FLAGS = flags[["CXXFLAGS"]],
        CXX14FLAGS = flags[["CXXFLAGS"]],
        CXX17FLAGS = flags[["CXXFLAGS"]],
        CXX20FLAGS = flags[["CXXFLAGS"]],
        FFLAGS = flags[["FFLAGS"]],
        FCFLAGS = flags[["FCFLAGS"]]
      ),
      debug = debug
    )
  } else {
    pkgbuild_with_debug(
      pak::pkg_install(paste0(repo, "?source&reinstall"), ask = FALSE),
      flags = c(
        CFLAGS = flags[["CFLAGS"]],
        CXXFLAGS = flags[["CXXFLAGS"]],
        CXX11FLAGS = flags[["CXXFLAGS"]],
        CXX14FLAGS = flags[["CXXFLAGS"]],
        CXX17FLAGS = flags[["CXXFLAGS"]],
        CXX20FLAGS = flags[["CXXFLAGS"]],
        FFLAGS = flags[["FFLAGS"]],
        FCFLAGS = flags[["FCFLAGS"]]
      ),
      debug = debug
    )
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
    CPPFLAGS <- paste0(
      "-std=c++17",
      " -I", R.home("include"),
      " -I", file.path(system.file(package = "Rcpp"), "include"),
      " -I/usr/local/include",
      " -Imisc/src_standalone",
      " -DNDEBUG"
    )
    SOURCE <-
      paste0(c(
        file.path("src", files),
        file.path("misc", "src_standalone", standalone_files)
      ), collapse = " ")
  } else {
    CPPFLAGS <- paste0(
      "-std=c++17",
      " -I", R.home("include"),
      " -I", file.path(system.file(package = "Rcpp"), "include"),
      " -I/usr/local/include",
      " -DNDEBUG -DFMESHER_WITH_R"
    )
    SOURCE <- paste0(file.path("src", files), collapse = " ")
  }
  cmd <- paste("clang-tidy", SOURCE, "--", CPPFLAGS)
  print(cmd)
  system(cmd, intern = TRUE)
}
