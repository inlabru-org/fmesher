// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// C_qinv
Rcpp::List C_qinv(SEXP AA);
RcppExport SEXP _fmesher_C_qinv(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(C_qinv(AA));
    return rcpp_result_gen;
END_RCPP
}
// fmesher_globe_points
Rcpp::NumericMatrix fmesher_globe_points(Rcpp::IntegerVector globe);
RcppExport SEXP _fmesher_fmesher_globe_points(SEXP globeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type globe(globeSEXP);
    rcpp_result_gen = Rcpp::wrap(fmesher_globe_points(globe));
    return rcpp_result_gen;
END_RCPP
}
// fmesher_rcdt
Rcpp::List fmesher_rcdt(Rcpp::List options, Rcpp::NumericMatrix loc, Rcpp::Nullable<Rcpp::IntegerMatrix> tv, Rcpp::Nullable<Rcpp::IntegerMatrix> boundary, Rcpp::Nullable<Rcpp::IntegerMatrix> interior, Rcpp::Nullable<Rcpp::IntegerVector> boundary_grp, Rcpp::Nullable<Rcpp::IntegerVector> interior_grp);
RcppExport SEXP _fmesher_fmesher_rcdt(SEXP optionsSEXP, SEXP locSEXP, SEXP tvSEXP, SEXP boundarySEXP, SEXP interiorSEXP, SEXP boundary_grpSEXP, SEXP interior_grpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type options(optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type loc(locSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerMatrix> >::type tv(tvSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerMatrix> >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerMatrix> >::type interior(interiorSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type boundary_grp(boundary_grpSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type interior_grp(interior_grpSEXP);
    rcpp_result_gen = Rcpp::wrap(fmesher_rcdt(options, loc, tv, boundary, interior, boundary_grp, interior_grp));
    return rcpp_result_gen;
END_RCPP
}
// fmesher_bary
Rcpp::List fmesher_bary(Rcpp::NumericMatrix mesh_loc, Rcpp::IntegerMatrix mesh_tv, Rcpp::NumericMatrix loc, Rcpp::List options);
RcppExport SEXP _fmesher_fmesher_bary(SEXP mesh_locSEXP, SEXP mesh_tvSEXP, SEXP locSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mesh_loc(mesh_locSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type mesh_tv(mesh_tvSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type loc(locSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(fmesher_bary(mesh_loc, mesh_tv, loc, options));
    return rcpp_result_gen;
END_RCPP
}
// fmesher_fem
Rcpp::List fmesher_fem(Rcpp::NumericMatrix mesh_loc, Rcpp::IntegerMatrix mesh_tv, int fem_order_max, Rcpp::List options);
RcppExport SEXP _fmesher_fmesher_fem(SEXP mesh_locSEXP, SEXP mesh_tvSEXP, SEXP fem_order_maxSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mesh_loc(mesh_locSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type mesh_tv(mesh_tvSEXP);
    Rcpp::traits::input_parameter< int >::type fem_order_max(fem_order_maxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(fmesher_fem(mesh_loc, mesh_tv, fem_order_max, options));
    return rcpp_result_gen;
END_RCPP
}
// fmesher_split_lines
Rcpp::List fmesher_split_lines(Rcpp::NumericMatrix mesh_loc, Rcpp::IntegerMatrix mesh_tv, Rcpp::NumericMatrix loc, Rcpp::IntegerMatrix idx, Rcpp::List options);
RcppExport SEXP _fmesher_fmesher_split_lines(SEXP mesh_locSEXP, SEXP mesh_tvSEXP, SEXP locSEXP, SEXP idxSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mesh_loc(mesh_locSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type mesh_tv(mesh_tvSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type loc(locSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type options(optionsSEXP);
    rcpp_result_gen = Rcpp::wrap(fmesher_split_lines(mesh_loc, mesh_tv, loc, idx, options));
    return rcpp_result_gen;
END_RCPP
}
// C_matrixio_test2
Rcpp::List C_matrixio_test2(Rcpp::List args_input);
RcppExport SEXP _fmesher_C_matrixio_test2(SEXP args_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_input(args_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(C_matrixio_test2(args_input));
    return rcpp_result_gen;
END_RCPP
}
// C_matrixio_test
Rcpp::List C_matrixio_test(Rcpp::List args_input);
RcppExport SEXP _fmesher_C_matrixio_test(SEXP args_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args_input(args_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(C_matrixio_test(args_input));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fmesher_C_qinv", (DL_FUNC) &_fmesher_C_qinv, 1},
    {"_fmesher_fmesher_globe_points", (DL_FUNC) &_fmesher_fmesher_globe_points, 1},
    {"_fmesher_fmesher_rcdt", (DL_FUNC) &_fmesher_fmesher_rcdt, 7},
    {"_fmesher_fmesher_bary", (DL_FUNC) &_fmesher_fmesher_bary, 4},
    {"_fmesher_fmesher_fem", (DL_FUNC) &_fmesher_fmesher_fem, 4},
    {"_fmesher_fmesher_split_lines", (DL_FUNC) &_fmesher_fmesher_split_lines, 5},
    {"_fmesher_C_matrixio_test2", (DL_FUNC) &_fmesher_C_matrixio_test2, 1},
    {"_fmesher_C_matrixio_test", (DL_FUNC) &_fmesher_C_matrixio_test, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fmesher(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}