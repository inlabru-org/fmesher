#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "fmesher_helpers.h"

#include "Rcpp.h"
#include "RcppEigen.h"


//' Command line test
//'
//' @param args_input Input argument list
// [[Rcpp::export]]
Rcpp::List C_cmdline_test(Rcpp::StringVector args_input)
{
  Rcpp::List ret;

  int argc = args_input.size();
  char** argv;
  argv = (char**)calloc(argc + 1, sizeof(char*));
  if (!argv) {
    ret["error"] = "calloc for argv failed";
    return ret;
  }
  argv[argc] = NULL;
  for (size_t i = 0; i < argc; i++) {
    argv[i] = (char*)calloc(args_input(i).size() + 1, sizeof(char));
    if (!argv[i]) {
      ret["error"] = "calloc for some argv[i] failed";
      // TODO: memory cleanup
      return ret;
    }
    strcpy(argv[i], Rcpp::as< std::string >(args_input(i)).c_str());
  }

  gengetopt_args_info args_info;
  struct cmdline_params params;

  FMLOG("checkpoint 1." << std::endl);

  cmdline_init(&args_info);
  cmdline_params_init(&params);

  FMLOG("checkpoint 2." << std::endl);

  /* call the command line parser */
  if (cmdline_ext(argc, argv, &args_info, &params) != 0) {
    cmdline_free(&args_info);
    FMLOG("cmdline failed." << std::endl);
    ret["error"] = "cmdline failed.";
    return ret;
  }

    if (args_info.dump_config_given) {
  FILE* f = fopen("testtest.cfg", "w");
  cmdline_dump(f,&args_info);
  fclose(f);
  }

  ret["narg"] = args_input.size();

//  cmdline_free(&args_info);
//  for (size_t i = 0; i < argc; i++) {
//    free(argv[i]);
//  }
//  free(argv);

  return ret;
}

//' Triangulate
//'
//' @param args_input Input argument list
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_triangulate(Rcpp::List args_input)
{

  return Rcpp::List::create();
}

//' Test the matrix I/O system
//'
//' @param args_input Input argument list
// [[Rcpp::export]]
Rcpp::List C_matrixio_test(Rcpp::List args_input)
{
  //Eigen::SparseMatrix<double> C_qinv(SEXP AA)
  using Eigen::MappedSparseMatrix;
  using Eigen::SparseMatrix;

  Rcpp::List ret;
  return(ret);
}
