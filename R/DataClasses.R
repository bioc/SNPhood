

#### S4 class definition
#' A class to represent, investigate, quantify and visualise the epigenomic neighbourhood of SNPs using NGS data
#' 
#' The class \code{SNPhood} stores read count-derived information from NGS files for a set of genomic regions of interest as well as associated metadata. 
#' It may additionally contain results of various subsequent analyses and statistical tests. 
#' See the description below or the Vignette for more details.
#' @section Constructors:
#' Currently, a \code{SNPhood} object can only be constructed by executing the main function of the package, \code{\link{analyzeSNPhood}}.
#' @section Accessors:
#' In the following code snippets, \code{SNPhood.o} is a \code{SNPhood} object and
#' \code{readGroupCur} and \code{datasetCur} a particular read group and dataset as defined in \code{SNPhood.o}, respectively.
#' 
#'# Get general annotation of a SNPhood object
#' 
#' \code{annotation(SNPhood.o)}: Get the annotation information, a nested list with multiple components (see names(annotation(SNPhood.o))).
#' 
#'# Get more specific annotation such as number and annotation of regions, datasets, bins, and read groups
#' 
#' \code{nRegions(SNPhood.o)}: Get the number of user regions.
#' 
#' \code{nDatasets(SNPhood.o)}: Get the number of datasets.
#' 
#' \code{nBins(SNPhood.o)}: Get the number of bins.
#' 
#' \code{nReadGroups(SNPhood.o)}: Get the number of read groups.
#' 
#' \code{annotationRegions(SNPhood.o)}: Get the annotation of user regions. 
#' 
#' \code{annotationDatasets(SNPhood.o)}: Get the annotation of datasets. 
#' 
#' \code{annotationBins(SNPhood.o)}: Get the annotation of bins. 
#' 
#' \code{annotationReadGroups(SNPhood.o)}: Get the annotation of read groups. 
#' 
#'# Get the parameters that were used for the analysis
#' 
#' \code{parameters(SNPhood.o)}: Get the parameter information, a nested list with multiple components (see names(parameters(SNPhood.o))).
#' 
#'# Get counts before binning 
#'
#' \code{counts(SNPhood.o, type = "unbinned", readGroup = readGroupCur, dataset = datasetCur)}: Get the counts for each user region before binning.
#' See ?counts for more details.
#' 
#' # If applicable, get counts after binning
#' 
#' \code{counts(SNPhood.o, type = "binned", readGroup = readGroupCur, dataset = datasetCur)}: Get the counts for each user region after binning.
#' See ?counts for more details.
#' 
#' # If applicable, get enrichment after binning
#' 
#' \code{enrichment(SNPhood.o, type = "binned", readGroup = readGroupCur, dataset = datasetCur)}: Get the enrichment for each user region after binning.
#' See ?enrichment for more details.
#' 
#' \strong{In addition, see the workflow vignette (browseVignettes(\"SNPhood\") for a full workflow that uses all accessors.}
#'
#' @slot annotation Named list. Contains various annotation and metadata such as:
#' \itemize{ 
#' \item \code{regions}: An object of class \code{GenomicRanges} that contains the user regions, including annotation and the position 
#' of the original user-provided position before creating regions and bins.  
#' \item \code{genotype}: A list one or two elements, both of which contain genotype-related information, either directly from the sequencing reads
#' or externally derived from a VCF file using the function \code{\link{associateGenotypes}}.
#' \item \code{readGroups}: The names of the read groups that are currently defined.
#' \item \code{files}: Contains a named list with additional information about each processed file, such as type(\code{signal} or \code{input}), \code{files}(a vector of one or multiple filenames), and \code{composite}(\code{TRUE} or \code{FALSE}, indicating if this is a composite file from multiple individual files)
#' }
#' Elements from this slot can be retrieved with the accessor function \code{\link{annotation}}.

# 
#' @slot config Named list. Named list with the parameters as specified in the parameter list and additionally the specific parameters
#'  the function \code{\link{analyzeSNPhood}} was called with (such as \code{onlyPrepareForDatasetCorrelation} and \code{input}). Elements from
#'  this slot can be retrieved with the accessor function \code{\link{parameters}}.

#' @slot readCountsUnbinned Named list (nested). Contains vectors of raw reads counts for each user region (before binning). 
#' The names of the list are the read groups and the filenames of the annotated datasets. Elements from
#'  this slot can be retrieved with the accessor function \code{\link{counts}} using \code{type} = "unbinned".

#' @slot readCountsBinned Named list (nested). Each element contains a matrix of raw reads counts per user region and bin (i.e., after binning). 
#'  The names of the list are the read groups and the filenames of the annotated datasets. 
#' Contains the raw read counts if normalization among all datasets ha sbeen performed (parameter \code{normAmongEachOther} is set to FALSE)
#' and normalized read counts otherwise.\cr \cr 
#' If read counts are recorded allele-specifically (in the following snippet paternal, maternal and ambiguous) for each group, 
#' the structure therefore may look like this:
#' \itemize{
#' \item \code{paternal}:
#'  \itemize{
#'  \item \code{dataset ID 1}: Matrix of read counts for each user region across bins
#'  \item \code{dataset ID 2}: Matrix of read counts for each user region across bins
#'  \item ...
#'  }
#' \item \code{maternal}: See read group paternal, identical structure
#' \item \code{ambiguous}: See read group paternal, identical structure
#' }

#' @slot enrichmentBinned Named list. See the description for the slot \code{readCountsBinned}, 
#' with the only difference that this slot contains the enrichment after normalizing with an input rather than the read counts. 
#' If input normalization is turned off, this slot is empty.

#' @slot additionalResults Named list. Contains additional information from subsequent analyses such as allelic bias tests or results 
#' of the genotype analysis. Initially empty. Different functions write the results in this slot. 
#' Elements from this slot can be retrieved with the accessor function \code{\link{results}}.

#' @keywords SNPhood-class, SNPhood
#' @name SNPhood-class
#' @aliases SNPhood-class
#
setClass("SNPhood", representation = representation(
  annotation        = "list",
  config            = "list", 
  enrichmentBinned        = "list",
  readCountsUnbinned = "list",
  readCountsBinned    = "list",
  additionalResults = "list",
  internal            = "list"
)
)

# SLOT config #

#' Retrieve the parameters of an object.
#' 
#' @param object An object containing parameters with which it was created.
#' @param ... Additional arguments, for use in specific methods.
#' @docType methods
#' @rdname parameters-methods
#' @aliases parameters
#' @export
setGeneric("parameters", function(object, ...) standardGeneric("parameters")) 

#' Retrieve the parameters of a \code{SNPhood} object.
#' 
#' @docType methods
#' @rdname parameters-methods
#' @aliases parameters
#' @return A named list with all parameters and its current values of the \code{SNPhood} object.
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' parameters(SNPhood.o)
#' @export
setMethod("parameters", "SNPhood", function(object, ...) {object@config})


#' Extract count data from a \code{\link{SNPhood}} object.
#' 
#' \code{counts} extracts count data from a \code{\link{SNPhood}} object. The full count data or only a subset can be extracted 
#' by settings the parameters \code{type}, \code{readGroup} and \code{dataset} accordingly. Either the count data
#' for the unbinned or binned SNP regions can be extracted.
#' @template object
#' @param type Character(1). Default "binned". Either "binned" or "unbinned" to extract counts after or before binning the SNP regions, respectively.  
#' @template readGroup
#' @template dataset
#' @param ... not used
#' @return A named nested list with the requested count data, organized after read group and dataset.
#' @aliases counts
#' @docType methods
#' @rdname counts-method
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' str(counts(SNPhood.o))
#' str(counts(SNPhood.o, readGroup = "paternal", dataset = 1))
#' str(counts(SNPhood.o, readGroup = c("maternal", "paternal"), dataset = 1))
#' @seealso \code{\link{SNPhood}}, \code{\link{enrichment}}
#' @export
#' @importFrom BiocGenerics counts
setMethod("counts", "SNPhood", function(object, type = "binned", readGroup = NULL, dataset = NULL, ...) {.getCounts(object, type, readGroup, dataset)})

#' @import checkmate
.getCounts <- function(SNPhood.o, type, readGroup = NULL, dataset = NULL, ...) {
    
    .checkObjectValidity(SNPhood.o)
    
    # TODO: extend and enable datasets and not just one and readGroups and not just one
    
    assertChoice(type, c("binned", "unbinned", "enrichmentBinned"))
    assertSubset(readGroup, SNPhood.o@annotation$readGroups)
    assert(checkNull(dataset), 
           checkIntegerish(dataset, lower = 1, upper = length(SNPhood.o@readCountsUnbinned[[1]]), any.missing = FALSE, min.len = 1, unique = TRUE), 
           checkSubset(dataset, names(SNPhood.o@readCountsUnbinned[[1]])))
    
    
    
    # get the read groups the user does not want to have
    discardReadGroups = NULL
    if (!testNull(readGroup)) {
        discardReadGroups = setdiff(SNPhood.o@annotation$readGroups, readGroup)       
    }
    
    discardDatasets = NULL
    if (!testNull(dataset)) {
        
        if (testSubset(dataset, names(SNPhood.o@readCountsUnbinned[[1]]))) {
            discardDatasets = setdiff(names(SNPhood.o@readCountsUnbinned[[1]]), dataset)       
            
        } else {
            discardDatasets = setdiff(seq_len(length(SNPhood.o@readCountsUnbinned[[1]])), dataset)
        }
    }
    
    if (type == "unbinned") {
        counts.l = SNPhood.o@readCountsUnbinned
    } else if (type == "binned") {
        counts.l = SNPhood.o@readCountsBinned
    } else {
        counts.l = SNPhood.o@enrichmentBinned
    }
    
    
    
    if (!(testNull(readGroup) & testNull(dataset))) {
        
        # Filter read groups
        for (readGroupCur in discardReadGroups) {
            counts.l[[readGroupCur]] = NULL
        }
        
        # Filter datasets
        for (readGroupCur in names(counts.l)) {
            
            # Delete datasets, but presort so the indexes are still correct
            if (!is.null(discardDatasets)) {
                for (datasetCur in sort(discardDatasets, decreasing = TRUE)) {
                    counts.l[[readGroupCur]][[datasetCur]] = NULL
                }
            }
            
            
        } 
        
    }
    
    
    for (readGroupCur in names(counts.l)) {
        
        for (datasetCur in names(counts.l[[readGroupCur]])) {
            
            if (type == "binned" | type == "enrichmentBinned") {
                rownames(counts.l[[readGroupCur]] [[datasetCur]]) = annotationRegions(SNPhood.o)
                colnames(counts.l[[readGroupCur]] [[datasetCur]]) = SNPhood.o@annotation$bins
            } else {
                names(counts.l[[readGroupCur]] [[datasetCur]]) = annotationRegions(SNPhood.o)
                
            }
            
        }
        
    } 
    
    # Reduce the complexity of the reuslting count object
    # First for the read group
    if (length(counts.l) == 1) {
        counts.l = counts.l[[1]]
        
        # Repeat it for the dataset
        if (length(counts.l) == 1) {
            counts.l = counts.l[[1]]
        }
        
    } else { # If multiple read groups are defined, but only one dataset
       
        for (i in seq_len(length(counts.l))) {
            if (length(counts.l[[i]]) == 1) {
                counts.l[[i]] = counts.l[[i]][[1]]
            }
        }
        
    }
    
    if (length(counts.l) == 0) {
        
        if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
            warning("Returning an empty list, could not find the requested data in the object. The parameter \"onlyPrepareForDatasetCorrelation\" has been set to TRUE. Rrun the function analyzeSNPhood again and set the parameter to FALSE. See also the help pages.")
            
        } else {
            warning("Returning an empty list, could not find the requested data in the object. Did you ask for the correct type of data? See also the help pages.")
            
        }
    }
    
    return(counts.l)
    
}


#' Extract enrichment data from an object.
#' 
#' @param object An object containing enrichment information.
#' @param ... Additional arguments, for use in specific methods.
#' @return Enrichment of the object or the objects components.
#' @export
#' @docType methods
#' @rdname enrichment-methods
#' @aliases enrichment
setGeneric("enrichment", function(object, ...) standardGeneric("enrichment")) 

#' Extract enrichment data from a \code{\link{SNPhood}} object.
#' 
#' \code{enrichment} extracts enrichment data from a \code{\link{SNPhood}} object. The full count data or only a subset can be extracted 
#' by settings the parameters \code{type}, \code{readGroup} and \code{dataset} accordingly. Principally, either the count data
#' for the unbinned or binned SNP regions can be extracted.
#' @template readGroup
#' @template dataset
#' @return Named list with the requested enrichment matrices from the \code{\link{SNPhood}} object, organized by read group and dataset
#' @seealso \code{\link{counts}}
#' 
#' @export
#' @docType methods
#' @rdname enrichment-methods
#' @aliases enrichment
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' str(enrichment(SNPhood.o), list.len=5)
setMethod("enrichment", "SNPhood", function(object, readGroup = NULL, dataset = NULL, ...) {.getEnrichment(object, readGroup, dataset)})

.getEnrichment <- function(SNPhood.o, readGroup = NULL, dataset = NULL) {
    
    .getCounts(SNPhood.o, type = "enrichmentBinned", readGroup, dataset)
}



# SLOT ANNOTATION #

#' Retrieve the annotation of a \code{SNPhood} object.
#' 
#' Specific elements within the annotation slot may also be extracted by using the \code{elements} parameter.
#' @template object
#' @param elements Character. The name of the element(s) in the annotation slot to be extracted. 
#' If set to \code{NULL}, the full annotation slot is returned.
#' @param ... not supported
#' @export
#' @return If only a single value for \code{elements} is provided, the element is returned directly. 
#' If multiple values are provided, a named list with the requested elements is returned.
#' @docType methods
#' @rdname annotation-methods
#' @aliases annotation
#' @importFrom BiocGenerics annotation
#' @examples 
#' data(SNPhood.o, package="SNPhood")
#' annotation(SNPhood.o)
#' annotation(SNPhood.o, elements = "regions")
#' annotation(SNPhood.o, elements = c("regions", "bins"))
setMethod("annotation", "SNPhood", function(object, elements = NULL, ...) {.getAnnotation(object, elements)})

.getAnnotation <- function(SNPhood.o, elements) {
    
    .checkObjectValidity(SNPhood.o)
    
    result = NULL
    
    assertSubset(elements, names(SNPhood.o@annotation))
    
    if (!testNull(elements)) {
        
        result = SNPhood.o@annotation[elements]
        
        if (length(elements) == 1) {
            result = result[[elements]]
        } 

    } else {
        
        result = SNPhood.o@annotation
    }
    
    
    result
    
}



#' Get results of various analyses performed with a \code{SNPhood} object.
#' 
#' Return the results of a particular analyis that is stored in the \code{SNPhood} object.
#'
#' @template SNPhood
#' @param type Character(1). Name of analyses one wants to retrieve the results for. 
#' Currently supported are "allelicBias", "clustering", "genotype" and "samplesCorrelation".
#' @param elements Character. Default NULL. Which elements of the resulting list structure should be returned? 
#' If set to NULL, all elements will be returned. Otherwise, if names are provided, only the requested subset elements will be returned.
#' If type equals "allelicBias", valid values are "pValue", "confIntervalMin", "confIntervalMax", "fractionEstimate", "background", "FDR_results", and "parameters".
#' If type equals "clustering", valid values are the defined read groups in the object.
#' If type equals "genotype", valid values are "strongGenotypes", "weakGenotypes", and "invariantGenotypes".
#' If type equals "samplesCorrelation", valid values are "corTable", and "transl".

#' @return A list with the results of the requested analysis and elements within.
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' head(results(SNPhood.o, type="allelicBias", elements = "parameters"))
#' head(results(SNPhood.o, type="allelicBias"))
#' @export
#' @import checkmate
#' @importFrom utils head
results <- function(SNPhood.o, type, elements = NULL) {
    
    # Check types and validity of arguments 
    .checkObjectValidity(SNPhood.o)
    assertChoice(type, SNPhood.o@internal$addResultsElementsAdded)
    
    assert(checkNull(elements), 
           checkIntegerish(elements, any.missing = FALSE, lower = 1, upper = length(SNPhood.o@additionalResults[[type]])), 
           checkSubset(elements, names(SNPhood.o@additionalResults[[type]]))
    )
    
    if (type == "allelicBias") {
        
        if (testNull(SNPhood.o@additionalResults$allelicBias) | testNull(names(SNPhood.o@additionalResults$allelicBias))) {
            stop("Could not find results of allelic bias test, did you call testForAllelicBiases before?")
        }

        results = SNPhood.o@additionalResults$allelicBias
        
        # Add row and coumn names for each matrix. Only do it now to save space
        elemsToAnnotate = c("pValue", "confIntervalMin", "confIntervalMax", "fractionEstimate")

        
        for (elemCur in elemsToAnnotate) {
            for (datasetCur in names(results[[elemCur]])) {
                dimnames(results[[elemCur]][[datasetCur]]) = list(annotationRegions(SNPhood.o), SNPhood.o@annotation$bins)
            }
           
        }

        if (results$parameters$calculateBackgroundDistribution) {
            for (datasetCur in names(results[["background"]])) {
                dimnames(results[[elemCur]][[datasetCur]]) = list(annotationRegions(SNPhood.o), SNPhood.o@annotation$bins)
            }
        }
        

        if (!testNull(elements)) {
            for (elementCur in names(results)) {
                if (!elementCur %in% elements) results[[elementCur]] = NULL
            }
            
            results = results[[elements]]
        }

        
        return(results)
        
    } else if (type == "genotype") {
       
        if (testNull(SNPhood.o@additionalResults[[type]]) | testNull(names(SNPhood.o@additionalResults[[type]]))) {
            stop("Could not find results of genotype analysis. Did you run the function calculateWeakAndStrongGenotype before?")
        }

        
    } else if (type == "samplesCorrelation") {
        
        if (testNull(SNPhood.o@additionalResults[[type]]) | testNull(names(SNPhood.o@additionalResults[[type]]))) {
            stop("Could not find results of samplesCorrelation analysis. Did you run the function calcCorrelationSamples before?")
        }
        

 
    } else if (type == "clustering") {
        
        if (testNull(SNPhood.o@additionalResults[[type]]) | testNull(names(SNPhood.o@additionalResults[[type]]))) {
            stop("Could not find results of samplesCorrelation analysis. Did you run the function clusterCountMatrix before?")
        }
        
    }
    
    if (length(elements) == 0) {
        results = SNPhood.o@additionalResults[[type]]
    } else if (length(elements) == 1) {
        results = SNPhood.o@additionalResults[[type]][[elements]]
    } else {
        results = SNPhood.o@additionalResults[[type]][elements]
    }
    
    return(results)
 
}



# OTHER FUNCTIONS #


# setMethod("show",
#             "SNPhood",
#             function(object) {
#                 
#               SNPhood.o = object    
#               .checkObjectValidity(SNPhood.o)
#                 
#               cat("OBJECT OF CLASS \"", class(SNPhood.o),"\"\n", sep = "")
#               
#               nDatasets = nDatasets(SNPhood.o)
#               nAlleles  = nReadGroups(SNPhood.o)
#               nBins     = nBins(SNPhood.o)
#               nRegions  = nRegions(SNPhood.o)
#               
#               cat("\n || Object Summary ||\n")
#               
#               addedText = ""
#               if (SNPhood.o@internal$calcEnrichment) {
#                   addedText = "and enrichment data "
#               }
#               cat("  Stores read counts ",addedText, "for ", nRegions," SNP regions (each of which is divided into ", nBins,
#                   " bins) and \n  metadata for a total of ", nDatasets," datasets (individuals) across ", nAlleles," read group(s)\n", sep = "")
#               
#               cat("\n || SLOTS ||\n")
#               cat("\n")
#               
#               cat(" | Slot \"annotation\" |\n  Stores annotation for the SNP regions and files. See the accessor function *annotation* for details.\n", sep = "")
#                             
#               cat("\n")
#       
#               
#               cat(" | Slot \"config\" |\n", sep = "")
#               cat("  Stores parameters and the general configuration. See the accessor function *parameters* for details.\n", sep = "")
# 
#               
#               cat("\n")
#               
#               cat(" | Slot \"readCountsBinned\" |\n  Stores read counts for each SNP region and bin. See the accessor function *counts* for details.\n", sep = "")
#               
#               cat("\n")
#                             
#               cat(" | Slot \"readCountsUnbinned\" |\n  Stores the total read counts for each SNP region. See the accessor function *counts* for details.\n", sep = "")
#               
#               cat("\n")
#               
#               if (SNPhood.o@internal$calcEnrichment) {
#                 cat(" | Slot \"enrichmentBinned\" |\n  Stores enrichment values for each SNP region and bin if applicable. See the accessor function *enrichment* for details.\n", sep = "")
#               }
#               cat("\n")
# 
#               
#               cat(" | Slot \"additionalResults\" |\n  Stores additional results such as those from allelic bias tests or genotype analyses. See the accessor function results for details.\n", sep = "")
#               
#               
#             }
#             )


setMethod("show",
          "SNPhood",
          function(object) {
              
              SNPhood.o = object    
              .checkObjectValidity(SNPhood.o)
              
              SNPhood.o@internal$disableObjectIntegrityChecking = TRUE
              
              cat("Class:", class(SNPhood.o),"\n")
              
              nDatasetsCur = nDatasets(SNPhood.o)
              nAllelesCur  = nReadGroups(SNPhood.o)
              nBinsCur     = nBins(SNPhood.o)
              nRegionsCur  = nRegions(SNPhood.o)
              
              if (SNPhood.o@config$onlyPrepareForDatasetCorrelation)
                  cat("!NOTE: Object incomplete because the parameter \"onlyPrepareForDatasetCorrelation\" in the function \"analyzeSNPhood\" has been set to TRUE\n")
              

              if (SNPhood.o@internal$countType == "enrichment") {
                  type = "Enrichments over input"
              } else if(SNPhood.o@internal$countType == "readCountsRaw") {
                  type = "Raw read counts"
              } else {
                  type = "Normalized read counts"
              }
              cat("Summary:", type, "for", nRegionsCur, "user regions across",nDatasetsCur,"datasets,",nAllelesCur, "read groups and", nBinsCur,"bins\n")
              
              
              cat("Metadata:\n")
              
              if (nRegionsCur > 4) {
                  string = paste0(paste0(head(annotationRegions(SNPhood.o),2), collapse = " "), " ... ", paste0(tail(annotationRegions(SNPhood.o),2), collapse = " "))
              } else {
                  string = paste0(annotationRegions(SNPhood.o), collapse = " ")
              }
              cat(" User regions (",nRegionsCur,"): ", string, "\n", sep = "")
              
              if (nDatasetsCur > 4) {
                  string = paste0(paste0(head(annotationDatasets(SNPhood.o),2), collapse = " "), " ... ", paste0(tail(annotationDatasets(SNPhood.o),2), collapse = " "))
              } else {
                  string = paste0(annotationDatasets(SNPhood.o), collapse = " ")
              }
              cat(" Datasets (",nDatasetsCur,"): ", string, "\n", sep = "")
              
              
              if (nAllelesCur > 4) {
                  string = paste0(paste0(head(annotationReadGroups(SNPhood.o),2), collapse = " "), " ... ", paste0(tail(annotationReadGroups(SNPhood.o),2), collapse = " "))
              } else {
                  string = paste0(annotationReadGroups(SNPhood.o), collapse = " ")
              } 
              cat(" Read groups (",nAllelesCur,"): ", string, "\n", sep = "")
              
              if (nBinsCur > 4) {
                  string = paste0(paste0(head(annotationBins(SNPhood.o),2), collapse = " "), " ... ", paste0(tail(annotationBins(SNPhood.o),2), collapse = " "))
              } else {
                  string = paste0(annotationBins(SNPhood.o), collapse = " ")
              } 
              cat(" Bins (",nBinsCur,"): ", string, "\n", sep = "")
              
              if (!is.null(SNPhood.o@annotation$genotype$external)) {
                  
                  colnamesVCF  = names(SNPhood.o@annotation$genotype$external)[-(1:3)]
                  vcf.unique.l = strsplit(colnamesVCF, split = ":", fixed = TRUE)
                  vcf.unique = unique(sapply(vcf.unique.l, "[[", 1))
                  
                  if (length(vcf.unique) > 4) {
                      string = paste0(paste0(head(vcf.unique, 2), collapse = " "), " ... ", paste0(tail(vcf.unique, 2), collapse = " "))
                  } else {
                      string = paste0(vcf.unique, collapse = " ")
                  } 
                  
                  cat(" Genotypes (2): reads-derived + VCF-derived (", string, ")", sep = "")
              } else {
                  cat(" Genotypes (1): reads-derived", sep = "")
              }
              cat("\n")
              
              if (length(SNPhood.o@internal$addResultsElementsAdded) > 0) {
                  cat("Stores analyses results for:", paste0(SNPhood.o@internal$addResultsElementsAdded,collapse = ", "),"\n")
              }
              
              SNPhood.o@internal$disableObjectIntegrityChecking = FALSE
              
  
          }
)


# SNPhood OBJECT VALIDITY #


#' @import checkmate
#' @importFrom methods validObject
.checkObjectValidity <- function(object, verbose = FALSE) {

    assertClass(object, "SNPhood")
    
    if (!is.null(object@internal$disableObjectIntegrityChecking)) {
        if (object@internal$disableObjectIntegrityChecking == FALSE) {
            if (verbose) message("Check object integrity and validity. For large objects, this may take some time. Use the function changeObjectIntegrityChecking to disable this check for the object.")
            
        }
    }
    
    res = validObject(object, test = TRUE, complete = TRUE)
    if (testFlag(res)) {
        if (!res) {
            stop("The SNPhood object is not valid for the following reason(s): \n\n", res)
        }
    } else {
        stop("The SNPhood object is not valid for the following reason(s): \n\n", res)
    }   
    
}

#' @import checkmate
#' @import GenomicRanges
# @importFrom GenomicRanges mcols
.validSNPhoodObj <- function(object) {
    
    valid = TRUE
    msg   = NULL
    
    # Skip validity check if object has not yet been fully constructed
    if (object@internal$disableObjectIntegrityChecking) return(TRUE) 
    
    
    if (object@config$onlyPrepareForDatasetCorrelation) return(TRUE) 

    
    # SLOT ANNOTATION #
    
    
    nRows = length(object@annotation$regions)
    nCols = length(object@annotation$bins)
    nDatasets = length(object@annotation$files)
    
    if (testNull(object@annotation) | testNull(names(object@annotation))) {
        
        valid = FALSE
        msg = c(msg, "Slot \"annotation\" or its names must not be NULL.\n")
        
    } else {
        
        validNames = c("regions", "genotype", "files", "readGroups", "bins")
        if (!testSubset(validNames, names(object@annotation), empty.ok = FALSE)) {
            valid = FALSE
            msg = c(msg, "Slot \"annotation\" must contain the elements ",paste0(validNames, collapse = ","),".\n")
        }
        
        ## Read groups ##
        if (!testCharacter(object@annotation$readGroups, min.chars = 1, any.missing = FALSE, min.len = 1)) {
            valid = FALSE
            msg = c(msg, "Element \"readGroups\" in slot \"annotation\" must be a character vector with at least one element.\n")
        }  
        
        ## Files ##
        if (!testList(object@annotation$files, any.missing = FALSE, min.len = 1, types = c("list"))) {
            valid = FALSE
            msg = c(msg, "Element \"files\" in slot \"annotation\" must be a named list of lists with at least one element.\n")
        } 
        
        for (i in seq_len(length(object@annotation$files))) {
            if (!testList(object@annotation$files[[i]], any.missing = FALSE, min.len = 4, types = c("character", "logical", "list"))) {
                valid = FALSE
                msg = c(msg, "Each element in \"files\" in slot \"annotation\" must be a list with four elements of type \"character\", \"logical\" or \"list\".\n")
            } 
        }
        
        ## Bins ##
        if (!object@config$onlyPrepareForDatasetCorrelation) {
            
           
            if (!object@config$normByInput) {
                lenMatrix = ncol(object@readCountsBinned[[1]][[1]])
            } else {
                lenMatrix = ncol(object@enrichmentBinned[[1]][[1]])
            }
            
            if (!testCharacter(object@annotation$bins, min.chars = 1, any.missing = FALSE, len = lenMatrix)) {
                valid = FALSE
                msg = c(msg, "The elements bins in the slot \"annotation\" must be of type character with a minimal length of 1 and no missing elements.\n")
            }  
           
        }
       
        
        ## Genotype ##
        validNames = c("readsDerived", "external", "heterozygosity", "mismatches")
        if (!testSubset(names(object@annotation$genotype), validNames)) {
            
            valid = FALSE
            msg = c(msg, "Element \"genotype\" in slot annotation must must only contain the elements ",paste0(validNames, collapse = ","),".\n")
            
        }
            
        # 1. readsDerived
        for (i in seq_len(length(object@annotation$readGroups))) {
            
            if (!testSubset(names(object@annotation$genotype$readsDerived)[i], object@annotation$readGroups[i])) {
                valid = FALSE
                msg = c(msg, "The names of the read groups are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
            }
            
            for (j in seq_len(length(object@annotation$files))) {
                
                if (!testSubset(names(object@annotation$genotype$readsDerived[[i]])[j], names(object@annotation$files)[j])) {
                    valid = FALSE
                    msg = c(msg, "The names of the datasets are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
                }
                
                if (!testMatrix(object@annotation$genotype$readsDerived[[i]][[j]], any.missing = FALSE, nrows = 4, ncols = length(object@annotation$regions))) {
                    valid = FALSE
                    msg = c(msg, "The genotype distribution element is invalid in the slot \"annotation$genotype$readsDerived\".\n")
                }
                
                if (!testSubset(rownames(object@annotation$genotype$readsDerived[[i]][[j]]), c("A","C","G","T"))) {
                    valid = FALSE
                    msg = c(msg, "The names of the datasets are incorrect for the element readsDerived in the slot \"annotation$genotype\".\n")
                }
                
                if (!testIntegerish(object@annotation$genotype$readsDerived[[i]][[j]], lower = 0, any.missing = FALSE)) {
                    valid = FALSE
                    msg = c(msg, "The genotype distribution element is invalid in the slot \"annotation$genotype$readsDerived\".\n")
                }
                
            }
        }
        
        # 2. External
        
        if (!testNull(object@annotation$genotype$external)) {
        
            if (!testDataFrame(object@annotation$genotype$external, min.cols = 4, all.missing = FALSE, types = c("character", "logical"))) {
                valid = FALSE
                msg = c(msg, "The element genotype$external in the slot \"annotation\" must be data frame.\n")
            }
            
            validNames = c("alleleRef", "alleleAlt", "genotypeMismatch")
            if (!testSubset(colnames(object@annotation$genotype$external)[1:3],  validNames)) {
                valid = FALSE
                msg = c(msg, "The element \"external\" in the slot \"annotation\" must contain only contain the following elements: ",paste0(validNames, collapse = ","),".\n")
            }
           
            validNames = c()
            for (i in seq_len(length(object@annotation$files))) {
                validNames = c(validNames, paste0(object@annotation$files[[i]]$genotypeFile,":",object@annotation$files[[i]]$genotypeFileSampleName))
            }
            
            nColsGenotype = ncol(object@annotation$genotype$external)
            if (nColsGenotype > 3) {
                testSubset(colnames(object@annotation$genotype$external)[4:nColsGenotype],  validNames)
            }
            
            if (!testLogical(object@annotation$genotype$external$genotypeMismatch)) {
                valid = FALSE
                msg = c(msg, "The column genotypeMismatch in genotype$external in the slot \"annotation\" must be of type logical.\n")
            }
            
            validValues = c("A","C","G","T", NA)
            if (!testSubset(object@annotation$genotype$external$alleleRef, validValues)) {
                valid = FALSE
                msg = c(msg, "The column alleleRef in genotype$external in the slot \"annotation\" must only contain the values ",paste0(validValues, collapse = ","),".\n")
            }
        
        }
        
        #3. heterozygosity
        if (!testMatrix(object@annotation$genotype$heterozygosity, mode = "logical", nrows = nRows, ncols = nDatasets + 1, any.missing = TRUE)) {
            valid = FALSE
            msg = c(msg, "The element \"heterozygosity\" in genotype in the slot \"annotation\" must be a matrix containing only logical or NA values.\n")
        }

        
        
        ## Regions ##  
        if (!testClass(object@annotation$regions, "GRanges")) {
            valid = FALSE
            msg = c(msg, "Element \"regions\" in slot \"annotation\" must be an object of class \"GRanges\".\n")
        } 
        
        
        if (!assert(checkNull(object@annotation$regions), checkSubset(c("annotation", "SNPPos"), names(mcols(object@annotation$regions)), empty.ok = FALSE))) {
            valid = FALSE
            msg = c(msg, "Element \"regions\" in slot \"annotation\" must contain at least the metadata \"annotation\" and \"SNPPos\"\n")
        }
        
        if (length(object@annotation$regions) < 1) {
            valid = FALSE
            msg = c(msg, "Element \"regions\" in slot \"annotation\" must contain at least one SNP region\n")
        }  
    }
    
    
    # Slot config #
    
    if (testNull(object@config)) {
        valid = FALSE
        msg = c(msg, "Slot \"config\" must not be NULL\n")
    }
    
    
    if (!testList(object@config, min.len = 1)) {
        valid = FALSE
        msg = c(msg, "Slot \"config\" must be a list\n")
        
    } else {
        
        if (!testSubset(c("onlyPrepareForDatasetCorrelation", "input", names(.getListOfSupportedParameters())), 
                        names(object@config))) {
            valid = FALSE
            msg = c(msg, "Slot \"config\" must contain all required parameters as elements, but at least one is missing.\n")
        }
        
        if (!testDataFrame(object@config$input, types = c("logical","character"), min.cols = 4, min.rows = 1, col.names = "named") |
            !testSubset(c("signal", "input", "individual", "genotype"), colnames(object@config$input)) ) {
            
            valid = FALSE
            msg = c(msg, "Element \"input\" in \"config\" must contain the columns signal, input, genotype and individual.\n")
        }
        
        
    }

    # Slot INTERNAL #
    
    if (testNull(object@internal)) {
        valid = FALSE
        msg = c(msg, "Slot \"internal\" must not be NULL\n")
        
    } else {
        
        requiredElems = c("isAllelicRatio", 
                          "sizeFactors", 
                          "mergedReadGroups", 
                          "plot_origBinSNPPosition",
                          "plot_labelBins",
                          "countType",
                          "readWidth",
                          "readStartPos",
                          "globalBackground",
                          "disableObjectIntegrityChecking",
                          "addResultsElementsAdded"
                          )
        
        
        if (!testSubset(requiredElems, names(object@internal))) {
            valid = FALSE
            msg = c(msg, "Slot \"internal\" must contain the elements ",paste0(requiredElems,collapse = ","),".\n")
        }
        
        if (!testFlag(object@internal$isAllelicRatio)) {
            valid = FALSE
            msg = c(msg, "The element isAllelicRatio in the slot \"internal\" must be of type logical.\n")
        }
        
        if (!testList(object@internal$sizeFactors, any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
            !testSetEqual(object@annotation$readGroups, names(object@internal$sizeFactors))) {
            valid = FALSE
            msg = c(msg, "Elements sizeFactors in the slot \"internal\" must be a list of lists with the number of elements identical to the element readGroups in slot annotation.\n")
        }  
        
        if (!testFlag(object@internal$mergedReadGroups)) {
            valid = FALSE
            msg = c(msg, "The element mergedReadGroups in the slot \"internal\" must be of type logical.\n")
        }

        
        if (!testFlag(object@internal$disableObjectIntegrityChecking)) {
            valid = FALSE
            msg = c(msg, "The element disableObjectIntegrityChecking in the slot \"internal\" must be of type logical.\n")
        }
        
        if (!testCharacter(object@internal$addResultsElementsAdded, min.chars = 1, any.missing = FALSE)) {
            valid = FALSE
            msg = c(msg, "The element addResultsElementsAdded in the slot \"internal\" must be of type character.\n")
        }
        
        
        if (!testNumeric(object@internal$plot_labelBins, any.missing = FALSE, len = length(object@annotation$bins), unique = TRUE)) {
            valid = FALSE
            msg = c(msg, "The element plot_labelBins in the slot \"internal\" must be of type character.\n")
        }
        
        if (!testSubset(names(object@internal$plot_labelBins), seq_len(length(object@annotation$bins)))) {
            valid = FALSE
            msg = c(msg, "The names of the element plot_labelBins in the slot \"internal\" must be from 1:", length(object@annotation$bins),".\n")
        }
        
        if (!testNumber(object@internal$plot_origBinSNPPosition, lower = 1)) {
            valid = FALSE
            msg = c(msg, "The element plot_origBinSNPPosition in the slot \"internal\" must be a single number.\n")
        }
        
        validValues = c("readCountsRaw", "readCountsNormalized", "enrichment")
        if (!testSubset(object@internal$countType, validValues)) {
            valid = FALSE
            msg = c(msg, "The element countType in the slot \"internal\" must be one of ",paste0(validValues,collapse = ","), ".\n")
        }
        
        if (object@internal$countType == "enrichment" & !parameters(object)$normByInput) {
            valid = FALSE
            msg = c(msg, "The element countType in the slot \"internal\" must match with the value of the parameter normByInput.\n")
        }


        if (length(object@annotation$readGroups) > 1) {
            
            for (i in seq_len(length(object@annotation$readGroups))) {
                
                if (!testSubset(names(object@internal$readStartPos)[i], object@annotation$readGroups[i])) {
                    valid = FALSE
                    msg = c(msg, "The names of the read groups are incorrect for the element readStartPos in the slot \"internal\":",paste0(names(object@internal$readStartPos),collapse = ",")," but expected ",paste0(object@annotation$readGroups,collapse = ","),".\n")
                }
                
                for (j in seq_len(length(object@annotation$files))) {
                    
                    if (!testSubset(names(object@internal$readStartPos[[i]])[j], names(object@annotation$files)[j])) {
                        valid = FALSE
                        msg = c(msg, "The names of the datasets are incorrect for the element readStartPos in the slot \"internal\".\n")
                    }
                    
                    # Check length, should equal nRegions
                    if (length(object@internal$readStartPos[[i]][[j]]) != length(object@annotation$regions)) {
                        valid = FALSE
                        msg = c(msg, "The dimension of the element readStartPos are incorrect in the slot \"internal\".\n")
                    }
    
                
                    for (l in seq_len(length(object@annotation$regions))) {
                        if (!testInteger(object@internal$readStartPos[[i]][[j]][[l]], lower = 1, any.missing = FALSE)) {
                            valid = FALSE
                            msg = c(msg, "At least one element in readStartPos in the slot \"internal\" is invalid and does not contain a vector of type Integer.\n")
                        } 
                    }
                   
                }
            }
        }
    }
    
    
    # Slot additionalResults #
    
    #assertSubset(elemsToAnnotate, names(results), empty.ok = FALSE)
    #assertLogical(results$parameters$calculateBackgroundDistribution, len = 1, any.missing = FALSE)
    
    #TODO
    
    
    
    #####################################################################
    # Slot readCountsUnbinned and readCountsBinned and enrichmentBinned #
    #####################################################################
    
    if (testNull(object@readCountsUnbinned) | testNull(object@readCountsBinned) | testNull(object@enrichmentBinned)) {
        valid = FALSE
        msg = c(msg, "Slots \"readCountsUnbinned\", \"readCountsBinned\" , and \"enrichmentBinned\"  must not be NULL\n")
        
    } else {
        
        # Test if names are correct
        for (readGroup in object@annotation$readGroups) {
            
            testSetEqual(names(object@annotation$files), names(object@readCountsUnbinned[[readGroup]]))
            
            if (!object@config$onlyPrepareForDatasetCorrelation)  {
               
                 testSetEqual(names(object@annotation$files), names(object@readCountsBinned[[readGroup]])) 
                
                # Test enrichment slot also
                if (length(object@enrichmentBinned[[1]]) > 0) {
                    testSetEqual(names(object@annotation$files), names(object@enrichmentBinned[[readGroup]])) 
                }
            } 
                
            
           
        }
        
        if (!object@config$onlyPrepareForDatasetCorrelation)  {
        
            if (!testList(object@readCountsUnbinned, any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
                !testList(object@readCountsBinned,   any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
                !testList(object@enrichmentBinned,       any.missing = FALSE, len = length(object@annotation$readGroups), types = "list") |
                !testSetEqual(object@annotation$readGroups, names(object@readCountsUnbinned)) | 
                !testSetEqual(object@annotation$readGroups, names(object@readCountsBinned))    |
                !testSetEqual(object@annotation$readGroups, names(object@enrichmentBinned))     
            ) {
                valid = FALSE
                msg = c(msg, "Slots \"readCountsUnbinned\", \"readCountsBinned\" and \"enrichmentBinned\" must be named lists of lists with the number of elements identical to the element readGroups in slot annotation\n")
            }
            
            # test validity and equality of names within counts of regions and bins
            if (length(object@readCountsBinned[[1]]) > 0) {
                
                for (i in seq_len(length(object@annotation$readGroups))) {
                    if (!testSetEqual(names(object@readCountsUnbinned[[i]]), names(object@readCountsBinned[[i]]))) {
                        valid = FALSE
                        msg = c(msg, "Slots \"readCountsUnbinned\" and \"readCountsBinned\" must contain identical names.\n")
                        
                    }
                }
                
            }  
            
            # Integrity checks for the enrichmentBinned slot
            if (length(object@enrichmentBinned[[1]]) > 0) {
                
                for (i in seq_len(length(object@annotation$readGroups))) {
                    for (j in names(object@annotation$files)) {
                        
                        # Skip input files, enrichment has only been calculated for the signal files
                        if (object@annotation$files[[j]]$type == "input") {
                            next
                        }
                        
                        if (!j %in% names(object@enrichmentBinned[[i]])) {
                            valid = FALSE
                            msg = c(msg, "Could not find element ", j, " in slot \"enrichmentBinned\".\n")
                            
                        }
                        
                        
                    } # end all files
  
                } # end all read groups
                
            }  # end if enrichment slot not empty
        }

        
        
        anyMatrixValueMissingAllowed = FALSE
        if (object@internal$isAllelicRatio | object@config$normByInput) {
            anyMatrixValueMissingAllowed = TRUE
        } 
        
        # Test dimensions and matrix of all elements in readCountsUnbinned and readCountsBinned
        for (i in seq_len(length(object@annotation$readGroups))) {
            
            for (j in names(object@annotation$files)) {

                # Skip check under certain conditions due to memory saving
                if (!object@config$keepAllReadCounts & object@config$normByInput) {
                    
                    next
                }
                
                # Test slot readCountsUnbinned
                if (!testVector(object@readCountsUnbinned[[i]][[j]], strict = TRUE, any.missing = anyMatrixValueMissingAllowed, len = nRows)) {
                    valid = FALSE
                    msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )              
                }
                
                if (!object@internal$isAllelicRatio) {
                    if (!testNumeric(object@readCountsUnbinned[[i]][[j]], lower = 0, any.missing = FALSE, len = nRows)) {
                        valid = FALSE
                        msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))               
                    }
                } else {
                    if (!testNumeric(object@readCountsUnbinned[[i]][[j]], lower = 0, upper = 1, any.missing = TRUE, len = nRows)) {
                        valid = FALSE
                        msg = c(msg, paste0("Slot \"readCountsUnbinned\" must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))            
                    }
                }

                if (!object@config$onlyPrepareForDatasetCorrelation)  {
                    # Test slot readCountsBinned         
                    if (nCols > 1) {
                        if (!testMatrix(object@readCountsBinned[[i]][[j]], any.missing = anyMatrixValueMissingAllowed, nrows = nRows, ncols = nCols)) {
                            valid = FALSE
                            msg = c(msg, paste0("Slot \"readCountsBinned\" 1 must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )             

                        }
                    } else {
                        if (!testNumeric(object@readCountsBinned[[i]][[j]], any.missing = anyMatrixValueMissingAllowed, len = nRows)) {
                            valid = FALSE
                            msg = c(msg, paste0("Slot \"readCountsBinned\" 2 must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n") )             
                        }
                    }
                    
                   
                    
                    if (!object@internal$isAllelicRatio) {
                        
                        if (object@internal$countType == "readCountsRaw") {
                            if (!testIntegerish(object@readCountsBinned[[i]][[j]], lower = 0, any.missing = FALSE)) {
                                valid = FALSE
                                msg = c(msg, paste0("Slot \"readCountsBinned\" 3 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))              
                            }
                        } else {
                            if (!testNumeric(object@readCountsBinned[[i]][[j]], lower = 0, any.missing = TRUE)) {
                                valid = FALSE
                                msg = c(msg, paste0("Slot \"readCountsBinned\" 4 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
                            }
                        }
                        
                    } else {
                        
                        if (!testNumeric(object@readCountsBinned[[i]][[j]], lower = 0, upper = 1, any.missing = TRUE)) {
                            valid = FALSE
                            msg = c(msg, paste0("Slot \"readCountsBinned\" 5 must contain only numerical vectors of read counts for each element. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
                        }
                    }
 
                    
                } # end if prepare only for correlation
                
            } # end all files
            

            # Test slot enrichmentBinned
            if (object@config$normByInput) {
                
                for (j in names(object@annotation$files)) {
                
                    # Skip input files, enrichment has only been calculated for the signal files
                    if (object@annotation$files[[j]]$type == "input") {
                        next
                    }
                    
                    if (!testMatrix(object@enrichmentBinned[[i]][[j]], any.missing =  anyMatrixValueMissingAllowed, nrows = nRows, ncols = nCols)) {
                        valid = FALSE
                        msg = c(msg, paste0("Slot \"enrichmentBinned\" must contain only matrices of read counts for each SNP and bin. However, at position [[", i,"]] [[", j,"]], a violation was found.\n"))             
                    }
                    
                }
            }
            
            
            
        } # end all read groups
        
    }
    
    
    if (valid) TRUE else msg
    
}

setValidity("SNPhood", .validSNPhoodObj)        
