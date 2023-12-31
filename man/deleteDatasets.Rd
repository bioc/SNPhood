% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{deleteDatasets}
\alias{deleteDatasets}
\title{Delete a particular set of datasets from a \emph{SNPhood} object.}
\usage{
deleteDatasets(SNPhood.o, datasets = NULL, verbose = TRUE)
}
\arguments{
\item{SNPhood.o}{Object of class \code{\linkS4class{SNPhood}}}

\item{datasets}{Numeric or Character or NULL. Default NULL. Datasets that should be used for plotting, either specified as integer (such as 1, value must be
between 1 and the total number of datasets as defined in the object) or their annotation (name must appear in the dataset names as obtained via the function \code{annotationDatasets}). If set to NULL, all datasets will be considered.}

\item{verbose}{Logical(1). Default TRUE. Should the verbose mode (i.e., diagnostic messages during execution of the script) be enabled?}
}
\value{
an object of class \code{\link{SNPhood}} with the requested datasets removed from all slots.
}
\description{
\code{deleteDatasets} deletes a particular set of datasets from a \code{\link{SNPhood}} object.
Removal is irreversible. It is therefore recommended to save the resulting \code{\link{SNPhood}} object with a new name because the deleted datasets cannot be recovered.
}
\examples{
data(SNPhood.o, package="SNPhood")
SNPhood_mod.o = deleteDatasets(SNPhood.o, c(1,2))
}
\seealso{
\code{\link{deleteRegions}}, \code{\link{deleteReadGroups}}
}

