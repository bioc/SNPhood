% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{annotationRegions}
\alias{annotationRegions}
\title{Get the annotation of SNP regions for a \emph{SNPhood} object.}
\usage{
annotationRegions(SNPhood.o, asGRangesObj = FALSE, verbose = FALSE)
}
\arguments{
\item{SNPhood.o}{Object of class \code{\linkS4class{SNPhood}}}

\item{asGRangesObj}{Logical(1). Default FALSE. Should the full annotation be returned (as \code{GRanges} object) or only the annotation of the SNP regions (as character vector)?}

\item{verbose}{Logical(1). Default FALSE. Should the verbose mode (i.e., diagnostic messages during execution of the script) be enabled?}
}
\value{
If \code{asGRangesObj} is set to \code{TRUE}, a \code{GRanges} object is returned. Otherwise, a character vector with the currently stored SNP annotation is returned.
}
\description{
Return the annotation of the SNP regions that are defined in the \code{\linkS4class{SNPhood}} object.
}
\examples{
data(SNPhood.o, package="SNPhood")
IDs.vec = annotationRegions(SNPhood.o, asGRangesObj = FALSE)
IDs.gr  = annotationRegions(SNPhood.o, asGRangesObj = TRUE)
}

