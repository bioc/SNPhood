% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{plotClusterAverage}
\alias{plotClusterAverage}
\title{Visualize average enrichment per cluster}
\usage{
plotClusterAverage(SNPhood.o, readGroup, dataset, fileToPlot = NULL,
  returnOnlyPlotNotObject = FALSE, verbose = FALSE)
}
\arguments{
\item{SNPhood.o}{Object of class \code{\linkS4class{SNPhood}}}

\item{readGroup}{Character(1). Default NULL. Read group that should be plotted, specified by its name as obtained by the function \code{annotationReadGroups}). If only one read group is defined in the object, this may also be NULL for user conveniance.}

\item{dataset}{Numeric(1) or Character(1). Single dataset that should be used for plotting, either specified as integer (such as 1, value must be
between 1 and the total number of datasets as defined in the object) or its annotation (name must appear in the dataset names as obtained via the function \code{annotationDatasets}).}

\item{fileToPlot}{Character(1) or \code{NULL}. Default \code{NULL}. 
Filename of the PDF file for the output plots.
If set to \code{NULL}, plots will be plotted to the currently active device.}

\item{returnOnlyPlotNotObject}{Logical(1). Default FALSE. If set to TRUE, only the plots are returned but not the actual object. 
Otherwise, for consistancy among the various visualization functions, the \code{\linkS4class{SNPhood}}  object is always returned, while the plots are
either written to a PDF file as specified by the parameter \code{fileToPlot} and/or to the currently active graphics device (i.e., the console usually)}

\item{verbose}{Logical(1). Default FALSE. Should the verbose mode (i.e., diagnostic messages during execution of the script) be enabled?}
}
\value{
the generated \pkg{ggplot2} plot(s) as list for further processing. May contain multiple plots, depending on the function. The plot(s) can then be plotted individually or modified arbitrarily as the user wants. For example, if multiple plots are returned and the plots have been saved in a variable called plots.l, simply type plots.l[[1]] to view the first plot.
}
\description{
\code{plotClusterAverage} visualizes the average reads per cluster. Note that the function \code{plotAndClusterMatrix} has to be executed
before \code{plotClusterAverage} is called for the same read group and dataset
}
\examples{
data(SNPhood.o, package="SNPhood")
plot = plotClusterAverage(SNPhood.o, readGroup = "paternal", dataset = 1)
}
\seealso{
\code{plotAndClusterMatrix}
}

