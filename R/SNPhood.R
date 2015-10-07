#' SNPhood: Investigate, quantify and visualise the epigenomic neighbourhood of SNPs using NGS data
#'
#'For more information and an introduction to the package, see the two vignettes.

#' @section \code{SNPhood} functions:
#' \code{\link{analyzeSNPhood}}
#' \code{\link{annotation}}
#' \code{\link{annotationBins}}
#' \code{\link{annotationBins2}}
#' \code{\link{annotationDatasets}}
#' \code{\link{annotationReadGroups}}
#' \code{\link{annotationRegions}}
#' \code{\link{associateGenotypes}}
#' \code{\link{collectFiles}}
#' \code{\link{convertToAllelicFractions}}
#' \code{\link{counts}}
#' \code{\link{deleteDatasets}}
#' \code{\link{deleteReadGroups}}
#' \code{\link{deleteRegions}}
#' \code{\link{enrichment}}
#' \code{\link{getDefaultParameterList}}
#' \code{\link{mergeReadGroups}}
#' \code{\link{nBins}}
#' \code{\link{nDatasets}}
#' \code{\link{nReadGroups}}
#' \code{\link{nRegions}}
#' \code{\link{parameters}}
#' \code{\link{plotAllelicBiasResults}}
#' \code{\link{plotAllelicBiasResultsOverview}}
#' \code{\link{plotAndCalculateCorrelationDatasets}}
#' \code{\link{plotAndCalculateWeakAndStrongGenotype}}
#' \code{\link{plotAndClusterMatrix}}
#' \code{\link{plotBinCounts}}
#' \code{\link{plotClusterAverage}}
#' \code{\link{plotGenotypesPerCluster}}
#' \code{\link{plotGenotypesPerSNP}}
#' \code{\link{plotRegionCounts}}
#' \code{\link{renameBins}}
#' \code{\link{renameDatasets}}
#' \code{\link{renameReadGroups}}
#' \code{\link{renameRegions}}
#' \code{\link{results}}
#' \code{\link{testForAllelicBiases}}
#' 
#' @section Contact Information: 
#' 
#' We value all the feedback that we receive and will try to reply in a timely manner.
#' Please report any bug that you encounter as well as any feature request that you may have to \email{SNPhood@@gmail.com}.

#' 
#' @return Summary analyses and visualizations for the selected genomic regions with respect to, for example, 
#' their read counts, genotype, and allelic origin
#' @docType package
#' @keywords SNPhood, SNPhood-package
#' @name SNPhood

NULL


#' SNPhood example data
#'
#' This dataset is an example dataset that can be used for exploring the \code{SNPhood} package.
#' For more information, see the workflow vignette of the \code{SNPhood} and \code{SNPhoodData} package, respectively.
#'
#' @docType data
#' @keywords datasets
#' @name SNPhood.o
#' @aliases  SNPhood-data
#' @return an example \code{SNPhood} object from the \code{SNPhoodData} package with read counts for 
#' 174 genomic regions across 2 datasets, three read groups and 100 bins
NULL
