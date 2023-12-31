% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{getDefaultParameterList}
\alias{getDefaultParameterList}
\title{Helper function to generate a default parameter list as input for the function \emph{analyzeSNPhood}}
\usage{
getDefaultParameterList(path_userRegions = NULL, isPairedEndData = TRUE)
}
\arguments{
\item{path_userRegions}{Character(1). Specify the value of the parameter \code{path_userRegions} 
(absolute path to the user regions file, see the Vignette for details).}

\item{isPairedEndData}{Logical(1). Default TRUE. Are the data paired-end (TRUE) or single-end (FALSE)?}
}
\value{
a named list with default values for the currently supported parameters that can be used as 
input for the function \code{\link{analyzeSNPhood}}:
\itemize{
 \item{readFlag_isPaired: Logical(1), TRUE for paired-end data, NA for single-end}
 \item{readFlag_isProperPair: Logical(1), TRUE}
 \item{readFlag_isUnmappedQuery: Logical(1), FALSE}
 \item{readFlag_hasUnmappedMate: Logical(1), FALSE}
 \item{readFlag_isMinusStrand: Logical(1), NA}
 \item{readFlag_isMateMinusStrand: Logical(1), NA}
 \item{readFlag_isFirstMateRead: Logical(1), NA}
 \item{readFlag_isSecondMateRead: Logical(1), NA}
 \item{readFlag_isNotPrimaryRead: Logical(1), FALSE}
 \item{readFlag_isNotPassingQualityControls: Logical(1), FALSE}
 \item{readFlag_isDuplicate: Logical(1), FALSE}
 \item{readFlag_reverseComplement: Logical(1), FALSE}
 \item{readFlag_simpleCigar: Logical(1), TRUE}
 \item{path_userRegions: Character(1), as given by the function argument path_userRegions}
 \item{zeroBasedCoordinates: Logical(1), FALSE}
 \item{regionSize: Integer(1), 500}
 \item{binSize: Integer(1), 50}
 \item{readGroupSpecific: Logical(1), TRUE}
 \item{strand: Character(1), "both"}
 \item{startOpen: Logical(1), FALSE}
 \item{endOpen: Logical(1), FALSE}
 \item{headerLine: Logical(1), FALSE}
 \item{linesToParse: Integer(1), -1}
 \item{lastBinTreatment: Character(1), "delete"}
 \item{assemblyVersion: Character(1), "hg19"}
 \item{nCores: Integer(1), 1}
 \item{keepAllReadCounts: Logical(1), FALSE}
 \item{normByInput: Logical(1), FALSE}
 \item{normAmongEachOther: Logical(1), TRUE}
 \item{poolDatasets: Logical(1), FALSE}
 
}
For reasons of reduced redundancy, a detailed description of the parameters can be found at the end of the 
main vignette in \code{SNPhood} (\code{browseVignettes("SNPhood")}).
}
\description{
\code{getDefaultParameterList} generates a default parameter list that can be used as input for the function \code{\link{analyzeSNPhood}}. 
The path to the user regions file can optionally be provided as an argument to the function. See the examples for further details.
Before running the function \code{\link{analyzeSNPhood}}, carefully check that the default parameters are suitable for the analysis.
}
\examples{
## Only one parameter can, optionally, be specified when calling the function
par.l = getDefaultParameterList(path_userRegions = "path/to/regions", isPairedEndData = TRUE)
## If the file is not specified, you need to change it 
## before you can execute the function \\code{\\link{analyzeSNPhood}}
par.l = getDefaultParameterList(isPairedEndData = TRUE)
par.l$path_userRegions = "path/to/regions"
}
\seealso{
\code{\link{analyzeSNPhood}}
}

