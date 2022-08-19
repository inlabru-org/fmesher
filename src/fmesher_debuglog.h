#ifndef FMESHER_DEBUGLOG_HH
#define FMESHER_DEBUGLOG_HH

#include <Rcpp.h>
#include <iostream>

// Define NDEBUG to disable assert
#include <cassert>

#ifndef WHEREAMI
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#endif

#ifndef FMLOG_
#define FMLOG_(msg) Rcpp::Rcout << WHEREAMI << msg;
#endif

#ifndef FMLOG
#ifdef FMESHER_DEBUG
#define FMLOG(msg) FMLOG_(msg)
#else
#define FMLOG(msg)
#endif
#endif

#endif
