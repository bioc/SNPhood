% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataClasses.R
\name{results}
\alias{results}
\title{Get results of various analyses performed with a \code{SNPhood} object.}
\usage{
results(SNPhood.o, type, elements = NULL)
}
\arguments{
\item{SNPhood.o}{Object of class \code{\linkS4class{SNPhood}}}

\item{type}{Character(1). Name of analyses one wants to retrieve the results for. 
Currently supported are "allelicBias", "clustering", "genotype" and "samplesCorrelation".}

\item{elements}{Character. Default NULL. Which elements of the resulting list structure should be returned? 
If set to NULL, all elements will be returned. Otherwise, if names are provided, only the requested subset elements will be returned.
If type equals "allelicBias", valid values are "pValue", "confIntervalMin", "confIntervalMax", "fractionEstimate", "background", "FDR_results", and "parameters".
If type equals "clustering", valid values are the defined read groups in the object.
If type equals "genotype", valid values are "strongGenotypes", "weakGenotypes", and "invariantGenotypes".
If type equals "samplesCorrelation", valid values are "corTable", and "transl".}
}
\value{
A list with the results of the requested analysis and elements within.
}
\description{
Return the results of a particular analyis that is stored in the \code{SNPhood} object.
}
\examples{
data(SNPhood.o, package="SNPhood")
head(results(SNPhood.o, type="allelicBias", elements = "parameters"))
head(results(SNPhood.o, type="allelicBias"))
}

