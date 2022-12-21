% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crs.R
\name{fm_CRS_as_list}
\alias{fm_CRS_as_list}
\alias{fm_list_as_CRS}
\alias{fm_CRSargs}
\alias{fm_CRSargs_as_list}
\alias{fm_list_as_CRSargs}
\title{Show expanded CRS arguments}
\usage{
fm_CRS_as_list(x, ...)

fm_list_as_CRS(x, ...)

fm_CRSargs(x, ...)

fm_list_as_CRSargs(x, ...)

fm_CRSargs_as_list(x, ...)
}
\arguments{
\item{x}{An \code{sp::CRS} or \code{inla.CRS} object (for
\code{fm_CRSargs} and \code{fm_CRS_as_list}), a character string (for
\code{fm_CRSargs_as_list}), or a list (for \code{fm_list_as_CRS} and
\code{fm_list_as_CRSargs}).}

\item{\dots}{Additional arguments passed on to other methods.}
}
\value{
For \code{fm_CRSargs} and \code{fm_list_as_CRSargs}, a character
string with PROJ.4 arguments.

For \code{fm_CRS_as_list} and \code{fm_CRSargs_as_list}, a list of
name/value pairs.

For \code{fm_list_as_CRS}, a \code{CRS} or \code{inla.CRS} object.

For \code{fm_list_as_CRSargs()}, a CRS proj4 string for name=value pair list

For \code{fm_CRSargs_as_list()}, a list of name=value pairs from CRS proj4 string
}
\description{
Wrappers for \code{sp::CRS} and \code{inla.CRS} objects to handle the
coordinate reference system argument string.
These methods should no longer be used with PROJ6/rgdal3;
see \code{\link[=fm_wkt]{fm_wkt()}} for a new approach.
}
\examples{

crs0 <- fm_CRS("longlat")
p4s <- fm_CRSargs(crs0)
lst <- fm_CRSargs_as_list(p4s)
crs1 <- fm_list_as_CRS(lst)
lst$a <- 2
crs2 <- fm_CRS(p4s, args = lst)
print(fm_CRSargs(crs0))
print(fm_CRSargs(crs1))
print(fm_CRSargs(crs2))
}
\seealso{
\code{\link[=fm_CRS]{fm_CRS()}}
}
\author{
Finn Lindgren \href{mailto:finn.lindgren@gmail.com}{finn.lindgren@gmail.com}
}
\keyword{internal}