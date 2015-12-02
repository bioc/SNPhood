
#' Calculate and plot correlation of region read counts among pairs of input files.
#'
#' \code{plotAndCalculateCorrelationDatasets} calculates and plots the pairwise correlation of all pairs of input files with among each other. 
#' The main purpose is to identify artefacts with particular files that should subsequently be excluded.
#' The correlation is based on the raw region read counts(i.e., before binning). 
#' The results of the correlation analysis are stored in the \code{\linkS4class{SNPhood}} object.
#' If the \code{corrplot} package is available, it will be used to produce a nice visualization of the correlation matrix.
#'
#' @template SNPhood
#' @template fileToPlot
#' @param corMeasure Character(1). Default "pearson". The correlation measure that should be used to compare between pairs of samples. 
#' Either \code{pearson}, \code{spearman}, or \code{kendall}. 
#' @param ... Additional arguments for the \code{corrplot.mixed} function from the \code{corrplot} package (if available).
#' @return An object of type \code{SNPhood}, with the results of the correlation analysis stored in the slot "additionalResults".
#' They can be retrieved via the helper function \code{results} for further investigation.
#' The results consist of a named list with two elements: A correlation matrix of the region read counts across all input files and a 
#' translation table to correlate the input files with the abbreviations from the correlation plot.
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' # Plot directly, using Pearson correlation
#' SNPhood.o = plotAndCalculateCorrelationDatasets(SNPhood.o)
#' # Plot to a PDF file
#' SNPhood.o = plotAndCalculateCorrelationDatasets(SNPhood.o, fileToPlot = "res.pdf")
#' # Using Spearman correlation instead of Pearson
#' SNPhood.o = plotAndCalculateCorrelationDatasets(SNPhood.o, corMeasure = "spearman")
#' @export
#' @import checkmate S4Vectors
#' @importFrom grDevices pdf dev.off
plotAndCalculateCorrelationDatasets <- function(SNPhood.o, fileToPlot = NULL, corMeasure = "pearson", ...) {
    
    # Check types and validity of arguments    
    .checkObjectValidity(SNPhood.o)
    assert(checkNull(fileToPlot), 
           checkCharacter(fileToPlot, min.chars = 1, any.missing = FALSE, len = 1))
    assertChoice(corMeasure, choices = c("pearson", "spearman", "kendall"))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
    }
    
    if (!is.null(fileToPlot)) pdf(fileToPlot)
    
    nRowsMatrix = nDatasets(SNPhood.o)
    nColsMatrix = nRegions(SNPhood.o)
    
    for (alleleCur in annotationReadGroups(SNPhood.o)) {
        
        # Create a matrix out of all region counts
        readCounts.m = matrix(0, nrow = nRowsMatrix, ncol = nColsMatrix) 
        rownames(readCounts.m) = names(SNPhood.o@readCountsUnbinned[[1]])
        
        for (i in seq_len(nRowsMatrix)) {
            
            readCounts.m[i,] = readCounts.m[i,]  + SNPhood.o@readCountsUnbinned[[alleleCur]][[i]]
        }
        
    }
    
    # Remove 0 rows
    rowSums = rowSums(readCounts.m)
    
    indexFiles_zeroReadCounts = which(rowSums == 0)
    
    if (length(indexFiles_zeroReadCounts) > 0) {
        warning("At least one file or individual contains only zero read counts, removing from correlation analysis(", rownames(readCounts.m)[indexFiles_zeroReadCounts],")")
        readCounts.m = readCounts.m[-which(rowSums == 0),]
    }
    
    # Remove 0 cols
    colSums = colSums(readCounts.m)
    
    nSNSPs_zeroCounts = length(which(colSums == 0))
    if (nSNSPs_zeroCounts > 0) {
        warning(nSNSPs_zeroCounts, " SNP region(s) with no read counts across all files, removing from correlation analysis...")
        readCounts.m = readCounts.m[,-which(colSums == 0)]
    }
    
    # Transpose the matrix as this is the required format
    readCounts.m = t(readCounts.m)
    
    
    M <- cor(readCounts.m, method = corMeasure)  
    M2 = M
    
    index_signal = 0
    index_input  = 0
    
    colnamesNew.vec = c()
    
    
    for (i in seq_len(ncol(M2))) {
        
        type = SNPhood.o@annotation$files[[colnames(M2)[i]]]$type
        stopifnot(!is.null(type))
        
        if (type == "signal") {
            index_signal = index_signal + 1
            colnamesNew.vec = c(colnamesNew.vec, paste0("S", index_signal))
            
        } else if (type == "input")  {
            index_input = index_input + 1
            colnamesNew.vec = c(colnamesNew.vec, paste0("I", index_input))
        } else {
            stop("Unknown type ", type," in slot annotation$files[[", i,"]], object inconsistent")
        }
    } 
    
    
    
    colnames(M2) = colnamesNew.vec
    rownames(M2) = colnamesNew.vec
    
    transl.df = data.frame(orig = rownames(M), plot = rownames(M2))
    
    # Finally, visualize
    
    if (requireNamespace("corrplot", quietly = TRUE)) {
        corrplot::corrplot.mixed(M2, ...)
    } else {
        ## When the corrplot package is not available
        warning("Default visualization not possible because the corrplot package is not available. You may consider installing it manually for improved visualization.")
        print(levelplot(M2), scales = list(y = c(-1,1)), xlab = "Datasets", ylab = "Datasets")
        
    }
    
    
    
    if (!is.null(fileToPlot)) dev.off()
    
    SNPhood.o@additionalResults$samplesCorrelation = list(corTable = M, transl = transl.df)
    
    SNPhood.o@internal$addResultsElementsAdded = c(SNPhood.o@internal$addResultsElementsAdded, "samplesCorrelation")
    
    
    SNPhood.o
    
}


#' Graphically summarize the results of the allelic bias analysis for a specific dataset and region.
#'
#' \code{plotAllelicBiasResults} graphically summarizes the results of the allelic bias analysis for a specific dataset and region.
#' Three plots are generated, each of which focuses on a different aspect of the allelic bias analysis across the selected user region.
#' 
#' The first plot shows the estimates of the allelic fraction, along with confidence intervals for the estimate
#' according to the parameters chosen when the function \code{testForAllelicBias} was called. Fraction estimates for which the corresponding p-values
#' are deemed significant according to the value of the parameter \code{signThreshold} are highlighted (see also the legend). At 0.5, the estimated
#' allelic fraction if there was no allelic bias, a horizontal line is drawn. 
#' 
#' The second plot shows the p-values (-log 10 transformed, so that smaller p-values have higher transformed values).
#' In analogy to the estimates of the allelic fraction, significant p-values are highlighted. The -log10 transformed significance threshold 
#' (according to the parameter \code{signThreshold}) appears as a horizontal line.
#' 
#' Finally, the third plot shows the distribution of the read counts across all read groups. In addition, the genotype distribution for each
#' read group is given (see the Vignette for details). This helps to identify allelic biases based on genotype differences among read groups.
#'
#' @template SNPhood
#' @template dataset
#' @template region
#' @template signThreshold
#' @template readGroupColors
#' @template fileToPlot
#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' SNPhood.o = testForAllelicBiases(SNPhood.o, readGroups = c("maternal", "paternal"))
#' # Leave all parameters with their standard values
#' plots = plotAllelicBiasResults(SNPhood.o)
#' 
#' # Change the colors
#' plots = plotAllelicBiasResults(SNPhood.o, readGroupColors = c("blue", "red", "gray"))
#' 
#' # Alter the significance threshold
#' plots = plotAllelicBiasResults(SNPhood.o, signThreshold = 0.01)
#' 
#' @export
#' @import checkmate
#' @import grid
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggplotGrob ggplot aes geom_point geom_line xlab ylab coord_cartesian ggtitle geom_hline theme element_text element_blank scale_colour_manual scale_fill_discrete scale_shape_manual labs 
plotAllelicBiasResults <- function(SNPhood.o, dataset = 1, region = 1, signThreshold = 0.05, readGroupColors = NULL, fileToPlot = NULL) {
    
    # Check types and validity of arguments    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    if (length(region) > 1) {
        stop("It is only possible to plot one region. Change the parameter region accordingly.")
    }
    
    if (!exists("allelicBias", where = SNPhood.o@additionalResults)) {
        stop("Element allelicBias not found in slot \"additionalResults\". Execute the function testForAllelicBiases first.")
    }
    
    region = .checkAndConvertRegionArgument(SNPhood.o, region, nullAllowed = FALSE, maxLength = 1) 
    dataset = .checkAndConvertDatasetArgument(SNPhood.o, dataset, nullAllowed = FALSE, maxLength = 1)
    
    assertPercentage(signThreshold)
    
    assert(checkNull(readGroupColors),
           checkCharacter(readGroupColors, min.chars = 1, len = nReadGroups(SNPhood.o))
    )
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
        pdf(fileToPlot)
    }
    
    #########################################
    # First plot with p value distribution #
    ######################################### 
    
    data.df = data.frame(bin = 1:nBins(SNPhood.o))
    data.df$value     = SNPhood.o@additionalResults$allelicBias$fractionEstimate[[dataset]][region,]
    data.df$confLower = SNPhood.o@additionalResults$allelicBias$confIntervalMin[[dataset]][region,]
    data.df$confUpper = SNPhood.o@additionalResults$allelicBias$confIntervalMax[[dataset]][region,]
    
    data.df$pvalue    = SNPhood.o@additionalResults$allelicBias$pValue[[dataset]][region,]
    data.df$sign      = data.df$pvalue <= signThreshold
    
    data2.df = data.df[which(data.df$sign),]
    
    data2.df = data.df
    
    mainLabel = paste0("Allelic bias results for the region around the SNP\n", .produceTitleForPlot(SNPhood.o, region))
    
    legendLabelConfInterval = paste0("Lower and upper\n",round(SNPhood.o@additionalResults$allelicBias$parameters$confLevel * 100, 0),"% CI")
    
    sizePoints = 3
    
    p1 <- ggplot(data.df, aes(x = bin, y = value, shape = as.factor(sign)))  + geom_point(colour = "black", size = sizePoints) 
    p1 <- p1 + .getVerticalLineForGGPlot(SNPhood.o@internal$plot_origBinSNPPosition)
    p1 <- p1 + scale_shape_manual(name = "Allelic fraction\nsignificant", values = c(1,19))
    p1 <- p1 + geom_line(aes(x = bin, y = confLower, colour = "lowerConf", shape = NULL))
    p1 <- p1 + geom_line(aes(x = bin, y = confUpper, colour = "lowerConf", shape = NULL))
    p1 <- p1 + ylab("Allelic fraction\n")
    p1 <- p1 + coord_cartesian(ylim = c(-0.1, 1.1))
    p1 <- p1 + ggtitle(mainLabel)
    p1 <- p1 + .getBinAxisLabelsForGGPlot(SNPhood.o@internal$plot_labelBins)
    
    p1 <- p1 + geom_hline(yintercept = 0.5, linetype = "dotted")
    p1 <- p1 + .getThemeForGGPlot()
    p1 <- p1 + theme(axis.title.x = element_blank(), plot.title = element_text(size = 12))
    p1 <- p1 + scale_colour_manual("Confidence\nintervals (CI)", 
                                   breaks = c("lowerConf"),
                                   labels = c(legendLabelConfInterval),
                                   values = c("gray"))
    
    
    
    #########################################
    # Second plot with p value distribution #
    #########################################
    
    data.df = data.frame(bin = 1:nBins(SNPhood.o))
    data.df$value       = SNPhood.o@additionalResults$allelicBias$pValue[[dataset]][region,]
    data.df$valueTransf = -log(data.df$value ,10)
    data.df$sign        = data.df$value <= signThreshold
    
    pValueSigThresholdTransformed = -log(signThreshold, 10)
    
    ylim = c(-0.1, max(max(data.df$valueTransf), pValueSigThresholdTransformed))
    ylim[2] = ylim[2] + 0.1 * ylim[2]
    
    p2 <- ggplot(data.df, aes(x=bin, y=valueTransf, shape=as.factor(sign)))  + .getVerticalLineForGGPlot(SNPhood.o@internal$plot_origBinSNPPosition) + ylab("-log10 p-value\n") + geom_point(colour = "red", size = 3)
    p2 <- p2 + scale_shape_manual(values = c(1,19))
    p2 <- p2 + coord_cartesian(ylim = ylim)
    p2 <- p2 + .getBinAxisLabelsForGGPlot(SNPhood.o@internal$plot_labelBins)
    
    p2 <- p2 + .getThemeForGGPlot()
    p2 <- p2 + theme(axis.title.x=element_blank())
    p2 <- p2 + geom_hline(yintercept = pValueSigThresholdTransformed, linetype = "dotted")
    p2 <- p2 + labs(shape = "p-value\nsignificant") + scale_fill_discrete(labels = c("no","yes"))
    
    ##################################################################
    # Third plot with a region summary for the particular individual #
    #################################### #############################
    
    p3 = plotBinCounts(SNPhood.o, region = region, readGroups = NULL, datasets = dataset, readGroupColors = readGroupColors, 
                       ylim = NULL, addGenotype = TRUE, addTitle = FALSE, plotGraph = FALSE)
    
    # Now force the width of all graphs to be identical so they align nicely and plot it using a grid
    gA <- ggplotGrob(p1)
    gB <- ggplotGrob(p2)
    gC <- ggplotGrob(p3)
    
    # Take width of the last plot as reference width
    widthPlot = gC$widths
    gA$widths <- widthPlot
    gB$widths <- widthPlot
    
    # Produce grid
    grid.arrange(gA, gB, gC, nrow = 3, newpage = TRUE)
    
    if (!testNull(fileToPlot)) dev.off()  
    
    list(plot1=p1, plot2=p2, plot3=p3)
    
}

#' Visualize counts or enrichment for a particular region across bins, datasets, and read groups.
#'
#' \code{plotBinCounts} visualizes counts or enrichment for a particular region across bins, datasets, and read groups. 
#' Many graphical parameters can be adjusted to suit the needs of the user, see below. 
#' 
#'
#' @template SNPhood
#' @template region
#' @template readGroups
#' @template datasets
#' @template readGroupColors
#' @template ylim
#' @param addGenotype Logical(1). Default TRUE. Should the genotype distribution for each read group at the original user position be displayed in the legend in addition?
#' See the Vignette for more details how this distribution is determined.
#' @param plotGenotypeRatio  Logical(1). Default FALSE. Should the ratio of the genotypes be plotted instead of the count or enrichment values?
#' Only applicable if the number of read groups to be plotted is 2. Setting this parameter to TRUE may result in ratios across bins that are interrupted due to
#' zero counts (and a resulting division by zero, which can therefore not be displayed). Also, ratios cannot be plotted if the genotype for the selected
#' regions could not be determined due to the lack of reads overlapping with the particular region (see the Vignette for details).
#' @param addTitle Logical(1). Default TRUE. Should the plot contain a title that summarizes the genomic region that is visualized?
#' @template colorPalette
#' @template plotGraph
#' @template fileToPlot
#' @template maxWidthLabels
#' 
#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' 
#' # Plot the first region, all parameters with their default values
#' plot = plotBinCounts(SNPhood.o)
#' 
#' # Plot the second region for the first dataset, using specific colors for the read groups.
#' plot = plotBinCounts(SNPhood.o, region = 2, dataset = 1, readGroupColors = c("red","blue","gray"))
#' 
#' # Plot the first region for the first dataset and the genotype ratio. Save the plot in a variable
#' plot = plotBinCounts(SNPhood.o, region = 1, readGroups = c("maternal", "paternal"), dataset = 1, plotGenotypeRatio = TRUE)
#' @export
#' @importFrom ggplot2 ggplot aes ggtitle xlab ylab geom_line coord_cartesian labs scale_colour_manual
plotBinCounts <- function(SNPhood.o, region = 1, readGroups = NULL, datasets = NULL, readGroupColors = NULL, ylim = NULL,
                          addGenotype = TRUE, plotGenotypeRatio = FALSE, addTitle = TRUE, colorPalette = "Set1", plotGraph = TRUE,  fileToPlot = NULL, maxWidthLabels = 25) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    if (testNull(readGroups)) {
        readGroups = annotationReadGroups(SNPhood.o)
    }
    
    if (testNull(datasets)) {
        datasets = annotationDatasets(SNPhood.o)
    }
    
    region = .checkAndConvertRegionArgument(SNPhood.o, region, nullAllowed = FALSE, maxLength = 1)
    
    assertSubset(readGroups, annotationReadGroups(SNPhood.o))
    
    assert(checkNull(readGroupColors),
           checkCharacter(readGroupColors, min.chars = 1, len = length(readGroups))
    )
    
    
    datasets = .checkAndConvertDatasetArgument(SNPhood.o, datasets, nullAllowed = FALSE, maxLength = NULL, returnNames = FALSE) 
    
    
    if (length(datasets) > 1 & !testNull(readGroupColors)) {
        warning("Set parameter readGroupColors to NULL because multiple datasets are going to be plotted. This parameter is only incorporated if a single dataset is plotted")
    }
    
    
    assertIntegerish(maxWidthLabels, lower = 1, any.missing = FALSE, len = 1)
    
    assert(checkNull(ylim),
           checkIntegerish(ylim, any.missing = FALSE, len = 2)
    )
    
    allowedPalettes = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
    assertChoice(colorPalette, allowedPalettes)
    
    assertFlag(addGenotype)
    assertFlag(plotGenotypeRatio)
    assertFlag(addTitle)
    assertFlag(plotGraph)
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
        pdf(fileToPlot)
    }
    
    
    if (plotGenotypeRatio & length(readGroups) != 2) {
        warning("Cannot plot genotype ratio as the number of read groups to plot is ",length(readGroups)," and not 2.")
        plotGenotypeRatio = FALSE
    }
    
    
    # Produce the data
    plot.df = data.frame(label = c(), dataset = c(), readGroup = c(), bin = c(), value = c())
    
    # Collect the most common genotype for each individual and read group
    mostCommonGenotype.df = data.frame(dataset = c(), readGroup = c(), base = c(), freq = c(), stringsAsFactors = FALSE)
    
    for (i in 1:length(readGroups)) {
        
        readGroupCur = readGroups[i]
        
        for (j in 1:length(datasets)) {
            
            datasetCur = datasets[j]
            
            # Construct label name
            label = paste0(readGroupCur)
            if (nReadGroups(SNPhood.o) == 1) label = "all" 
            
            if (length(datasets) > 1) {
                datasetName = strtrim(paste0(annotationDatasets(SNPhood.o)[datasetCur],": ",label), maxWidthLabels)
            } else {
                datasetName = strtrim(paste0(label), maxWidthLabels)
            }
            
            
            # Add genotype
            if (addGenotype) {
                genotype = SNPhood.o@annotation$genotype$readsDerived[[readGroupCur]][[datasetCur]][,region]
                genotype = sort(genotype,decreasing = TRUE)
                nonZeroReads = length(which(genotype > 0))
                
                genotypeStr = "NA"
                
                if (plotGenotypeRatio) {
                    mostCommonGenotype.df = rbind(mostCommonGenotype.df, 
                                                  data.frame(
                                                      dataset  = datasetCur, 
                                                      readGroup = readGroupCur,
                                                      base = ifelse(sum(genotype) > 0, names(genotype[1]), NA), 
                                                      freq = ifelse(sum(genotype) > 0, genotype[1], 0),
                                                      stringsAsFactors = FALSE
                                                  ))
                }
                
                if (nonZeroReads > 0) {
                    genotype_nonZero = genotype[which(genotype > 0)]
                    genotypeStr = paste0(names(genotype_nonZero),":",genotype_nonZero,collapse = ",")
                    
                }
                
                datasetName = paste0(datasetName, " (",genotypeStr,")")
            }
            
            # Add to df rowwise
            for (k in 1:nBins(SNPhood.o)) {
                plot.df = rbind(plot.df, data.frame(label = datasetName, dataset = datasetCur, readGroup = readGroupCur, bin = k, value = SNPhood.o@readCountsBinned[[readGroupCur]][[datasetCur]][region,k]))
            }
            
        }
        
    }
    
    if (plotGenotypeRatio) {
        
        excludeDatasets = c()
        
        # Get the most common genotype for each individual and determine which genotype should be the nominator
        summary.df = data.frame(dataset = datasets, refReadGroup = NA, mostFrequentBase = NA, altReadGroup = NA, leastFrequentBase = NA, stringsAsFactors = FALSE)
        
        for (datasetCur in datasets) {
            
            maxCur = which.max(mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$freq)
            minCur = which.min(mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$freq)
            summary.df$refReadGroup[which(summary.df$dataset == datasetCur)]      = mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$readGroup[maxCur]
            summary.df$mostFrequentBase[which(summary.df$dataset == datasetCur)]  = mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$base[maxCur]   
            summary.df$altReadGroup[which(summary.df$dataset == datasetCur)]      = mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$readGroup[minCur]
            summary.df$leastFrequentBase[which(summary.df$dataset == datasetCur)] = mostCommonGenotype.df[mostCommonGenotype.df$dataset == datasetCur,]$base[minCur]   
            
        }
        
        # Determine the two alleles for which the ratio should be plotted
        refBases = table(summary.df$mostFrequentBase)
        altBases = table(summary.df$leastFrequentBase)
        
        allBases = unique(c(summary.df$mostFrequentBase, summary.df$leastFrequentBase))
        allBases = allBases[!is.na(allBases)]
        
        indexRef = which(!names(refBases) %in% names(altBases))
        indexAlt = which(!names(altBases) %in% names(refBases))
        
        if (length(indexRef) == 0 & length(indexAlt) == 0) {
            
            # Select the base that is predominant in mostFrequentBase as refAllele
            refBase = names(refBases[1])
            altBase = allBases[which(allBases != refBase)][1]
            
            # if the SNP is homozygous for the given individual, altbase will be NA
            if (is.na(altBase)) {
                excludeDatasets = summary.df$dataset[summary.df$mostFrequentBase == summary.df$leastFrequentBase]
            }
            
        } else if (length(indexRef) != 0 & length(indexAlt) == 0) {
            
            refBase = names(refBases)[indexRef[1]]
            altBase = allBases[which(allBases != refBase)][1]
            
            
            
        } else if (length(indexRef) == 0 & length(indexAlt) != 0) {
            
            altBase = names(altBases)[indexAlt[1]]
            refBase = allBases[which(allBases != altBase)][1]
            
        } else {
            
            refBase = names(refBases)[indexRef[1]]
            altBase = names(altBases)[indexAlt[1]]
        }
        
        # Check if any datasets are excluded at this point
        if (length(excludeDatasets) > 0) {
            warning("Excluding dataset ", paste0(excludeDatasets, collapse = ","), " from plot because the SNP is homozygous.")
            summary.df = summary.df[-which(summary.df$dataset %in% excludeDatasets),]
        }
        
        # Delete datasets if they have missing values
        
        excludeDatasets = summary.df$dataset[unique(which(is.na(summary.df$mostFrequentBase) | is.na(summary.df$leastFrequentBase)))]
        
        if (length(excludeDatasets) > 0) {
            warning("Excluding dataset ", paste0(excludeDatasets, collapse = ","), " from plot because of missing genotypes.")
            summary.df = summary.df[-which(summary.df$dataset %in% excludeDatasets),]
        } 
        
        # Reset
        excludeDatasets = c()
        
        # Check the remaining datasets for further criteria. 
        # Consistency checks: Does the refbase appear in all the datasets at least in one read group? if not, cannot plot
        for (i in seq_len(nrow(summary.df))) {
            
            if ((summary.df$mostFrequentBase[i]  != refBase & summary.df$mostFrequentBase[i]  != refBase) |
                (summary.df$leastFrequentBase[i] != refBase & summary.df$leastFrequentBase[i] != altBase)) {
                excludeDatasets = c(excludeDatasets, summary.df$dataset[i])                
            }
            
        }
        
        if (length(excludeDatasets) > 0) {
            warning("Excluding dataset ", paste0(excludeDatasets, collapse=","), " from plot because of genotypes that differ from two selected genotypes for which the ratio is produced.")
            summary.df = summary.df[-which(summary.df$dataset %in% excludeDatasets),]
        } 
        
        
        if (nrow(summary.df) == 0) {
            
            plotGenotypeRatio = FALSE
            warning("No datasets left to plot the genotype ratio. Plotting normal read counts instead.")
            
        } else {
            
            # Construct the ratio, how to handle cases when the denominator is 0? Draw points not lines
            plot2.df = data.frame(label = c(), dataset = c(), bin = c(), value = c(), stringsAsFactors = FALSE)
            
            for (datasetCur in summary.df$dataset) {
                
                for (binCur in 1:nBins(SNPhood.o)) {
                    plot.subset.df = plot.df[plot.df$bin == binCur & plot.df$dataset == datasetCur,]
                    refReadGroup = summary.df$refReadGroup[which(summary.df$dataset == datasetCur)]
                    altReadGroup = summary.df$altReadGroup[which(summary.df$dataset == datasetCur)]
                    
                    ratio = plot.subset.df$value[which(plot.subset.df$readGroup == refReadGroup)] / plot.subset.df$value[which(plot.subset.df$readGroup == altReadGroup)] 
                    
                    # TODO: Ratio of one genotype with respect to the other
                    #                     if (plot) {
                    #                         ratio = 
                    #                     }
                    
                    if (!is.finite(ratio)) {
                        ratio = NA
                    }
                    
                    labelCur = annotationDatasets(SNPhood.o)[datasetCur]
                    
                    plot2.df = rbind(plot2.df, data.frame(label = labelCur, dataset = datasetCur, bin = binCur, value = ratio, stringsAsFactors = FALSE))
                    
                }
            }
            
            # Replace the original data frame by the new one
            plot.df = plot2.df
        }
        
        
        
    } # end if (plotGenotypeRatio)
    
    
    customYLimits = FALSE
    if (!testNull(ylim)) {
        customYLimits = TRUE
        
        if (ylim[1] > ylim[2]) {
            stop("Value for parameter ylim not correct: The second value must be larger than the first value.")
        }
        
    } else {
        ylim = c(-0.5, max(plot.df$value, na.rm = TRUE) + 1)
    }
    
    if (addTitle) {
        mainLabel = paste0(.getBinLabelYAxis(SNPhood.o), " for the region around the SNP\n", .produceTitleForPlot(SNPhood.o, region))
    }
    if (plotGenotypeRatio) {
        mainLabel = paste0("Genotype ratio of ",.getBinLabelYAxis(SNPhood.o)," for the region around the SNP\n", .produceTitleForPlot(SNPhood.o, region))
    }
    
    legendTitle = "Dataset: read group"
    
    if (length(datasets) == 1) {
        legendTitle = paste0("Dataset ",annotationDatasets(SNPhood.o)[datasets],":\nread group")
        
        if (addGenotype) {
            legendTitle = paste0(legendTitle, " (genotype)")
        }
        
    } else {
        if (addGenotype) {
            legendTitle = paste0(legendTitle, "\n(genotype)")
        }
    }
    
    if (plotGenotypeRatio) {
        
        legendTitle = paste0("Dataset")
        if (length(datasets) == 1) {
            legendTitle = annotationDatasets(SNPhood.o)[datasets]
        } 
        
    }
    
    ylabLabel = paste0(.getBinLabelYAxis(SNPhood.o),"\n")
    if (plotGenotypeRatio) {
        ylabLabel = paste0("Ratio of ",.getBinLabelYAxis(SNPhood.o)," with genotype ", refBase," over ", altBase, "\n(read groups ", paste0(readGroups, collapse = ","),")\n")
    }
    
    p <- ggplot(plot.df, aes(x = bin, y = value, colour = label)) + xlab(.getBinLabelXAxis(SNPhood.o)) + ylab(ylabLabel) + geom_line()
    p <- p + coord_cartesian(ylim = ylim)
    if (addTitle) p <- p + ggtitle(mainLabel)
    p <- p + .getBinAxisLabelsForGGPlot(SNPhood.o@internal$plot_labelBins)
    p <- p + .getVerticalLineForGGPlot(SNPhood.o@internal$plot_origBinSNPPosition)
    if (plotGenotypeRatio) p <- p + geom_hline(yintercept = 1, linetype = "dashed")
    p <- p + .getThemeForGGPlot()
    p <- p + labs(colour = legendTitle)
    
    if (length(datasets) == 1) {
        if (!testNull(readGroupColors)) p <- p + scale_colour_manual(values = readGroupColors)
    } else {
        p <- p + scale_colour_manual(values = .generateColorsForReadGroupsAndDatasets(length(datasets), length(readGroups), colorPalette, saturationMin = 0.3, deleteFirst = FALSE, deleteLast = TRUE))
    } 
    
    
    
    if (plotGraph) {
        
        if (!testNull(fileToPlot)) {
            print(p)
            dev.off()  
        } 
        
        return(print(p))
    } else {
        return(p)
    }
    
    
}

#' Visualize the results of the allelic bias analysis across regions or a user-defined genomic range
#'
#' \code{plotBinCounts} visualizes the results of the allelic bias analysis across regions or a user-defined genomic range.
#' Note that only the results of a particular chromosome can be visualized. It is therefore only possible if the regions to be visualized 
#' are located on one particular chromosome; otherwise, an error is thrown.
#'
#' @template SNPhood
#' @template regions
#' @template datasets 
#' @template plotChr
#' @template plotStartPos
#' @template plotEndPos
#' @template ylim
#' @template plotRegionBoundaries
#' @template plotRegionLabels
#' @template signThreshold
#' @param pValueSummary Character(1). Default "min". Either "min" or "median". If set to "min", for each region, the minimum p-value across all bins
#' is displayed as a representative result for the region. This is in analogy to how the background caclulation for the FDR calculation works,
#' see the vignette for details. If set to "median", the median p-value is calculated for each region and plotted. This may facilitate to identify
#' regions for which a lot of bins have low p-values.
#' @template maxWidthLabels
#' @template colorPalette
#' @template sizePoints
#' @template plotGraph
#' @template fileToPlot
#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' 
#' # Plot the allelic bias results for the first region using default values for all parameters
#' plots = plotAllelicBiasResultsOverview(SNPhood.o)
#' 
#' # Plot the allelic bias results for the full chr21
#' plots = plotAllelicBiasResultsOverview(SNPhood.o, regions = NULL, plotChr = "chr21")
#' 
#' @export
plotAllelicBiasResultsOverview <- function(SNPhood.o, regions = 1, datasets = NULL, plotChr = NULL, plotStartPos = NULL, plotEndPos = NULL, ylim = NULL, plotRegionBoundaries = FALSE, plotRegionLabels = FALSE,      
                                           signThreshold = 0.05, pValueSummary = "min", maxWidthLabels = 25, colorPalette = "Set1", sizePoints = 4, plotGraph = TRUE, fileToPlot = NULL) {
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    .plotRegionFeatures(SNPhood.o, regions = regions, readGroups = NULL, datasets = datasets, mergeReadGroupCounts = FALSE, plotChr = plotChr, plotStartPos = plotStartPos, plotEndPos = plotEndPos, ylim = ylim, plotRegionBoundaries = plotRegionBoundaries, plotRegionLabels = plotRegionLabels, plotAllelicBiasResults = TRUE,         
                        signThreshold = signThreshold, pValueSummary = pValueSummary, maxWidthLabels = maxWidthLabels, colorPalette = colorPalette,  sizePoints = sizePoints, plotGraph = plotGraph, fileToPlot = fileToPlot)
    
}  


#'  Visualize the raw read counts across regions or a user-defined genomic range
#'
#' \code{plotBinCounts} visualizes the raw read counts (i.e., before binning user regions) across regions or a user-defined genomic range.
#' Note that only the results of a particular chromosome can be visualized. It is therefore only possible if the regions to be visualized 
#' are located on one particular chromosome; otherwise, an error is thrown.
#' @template SNPhood
#' @template regions
#' @template datasets 
#' @template readGroups
#' @param mergeReadGroupCounts Logical(1). Default FALSE. Should the read groups be merged for visualization purposes?
#' @template plotChr
#' @template plotStartPos
#' @template plotEndPos
#' @template ylim
#' @template plotRegionBoundaries
#' @template plotRegionLabels
#' @template maxWidthLabels
#' @template colorPalette
#' @template sizePoints
#' @template plotGraph
#' @template fileToPlot
#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' 
#' # Plot the read counts for the first ten regions
#' plot = plotRegionCounts(SNPhood.o, regions = 1:10)
#' 
#' # Plot the read counts for the full chr21
#' plot = plotRegionCounts(SNPhood.o, plotChr = "chr21")
#' 
#' # Plot the read counts for the full chr21, merge read group counts and decrease the point size
#' plot = plotRegionCounts(SNPhood.o, plotChr = "chr21", sizePoints = 2, mergeReadGroupCounts = TRUE)
#' @export
#' @import checkmate 
plotRegionCounts <- function(SNPhood.o, regions = NULL, datasets = NULL, readGroups = NULL, mergeReadGroupCounts = FALSE, plotChr = NULL, plotStartPos = NULL, plotEndPos = NULL, ylim = NULL, plotRegionBoundaries = FALSE, plotRegionLabels = FALSE,           
                             maxWidthLabels = 25, colorPalette = "Set1",  sizePoints = 4, plotGraph = TRUE, fileToPlot = NULL) {
    
    assertFlag(mergeReadGroupCounts) 
    
    .plotRegionFeatures(SNPhood.o, regions = regions, readGroups = readGroups, datasets = datasets, mergeReadGroupCounts = mergeReadGroupCounts, plotChr = plotChr, plotStartPos = plotStartPos, plotEndPos = plotEndPos, ylim = ylim, plotRegionBoundaries =  plotRegionBoundaries, plotRegionLabels = plotRegionLabels, plotAllelicBiasResults = FALSE,         
                        maxWidthLabels = maxWidthLabels, colorPalette = colorPalette, sizePoints = sizePoints, plotGraph = plotGraph, fileToPlot = NULL)
}

#' @import checkmate S4Vectors
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes xlab ylab geom_point coord_cartesian ggtitle geom_vline labs scale_colour_manual scale_x_continuous theme
.plotRegionFeatures <- function(SNPhood.o, regions = NULL, readGroups = NULL, datasets = NULL, mergeReadGroupCounts = FALSE, plotChr = NULL, plotStartPos = NULL, plotEndPos = NULL, ylim = NULL, plotRegionBoundaries = FALSE, plotRegionLabels = FALSE, plotAllelicBiasResults = FALSE,         
                                signThreshold = 0.05, pValueSummary = "min", maxWidthLabels = 25, colorPalette = "Set1", sizePoints = 4, plotGraph = TRUE, fileToPlot = NULL) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (testNull(datasets)) {
        datasets = annotationDatasets(SNPhood.o)
    }
    
    readGroupsOrigLabel = NULL
    
    if (testNull(readGroups)) {
        readGroupsOrigLabel = "all"
        readGroups = annotationReadGroups(SNPhood.o)
    }
    
    assertSubset(readGroups, annotationReadGroups(SNPhood.o))
    
    regions = .checkAndConvertRegionArgument(SNPhood.o, regions, nullAllowed = TRUE, maxLength = NULL)
    datasets = .checkAndConvertDatasetArgument(SNPhood.o, datasets, nullAllowed = FALSE, maxLength = NULL, returnNames = FALSE) 
    
    
    genomeAssembly.df = .getGenomeData(SNPhood.o@config$assemblyVersion)
    
    assert(checkNull(plotChr),
           checkCharacter(plotChr, min.chars = 1, len = 1, any.missing = FALSE))
    
    if (!testNull(plotChr)) assertChoice(plotChr, as.character(genomeAssembly.df$chr))
    
    assert(checkNull(plotStartPos),
           checkIntegerish(plotStartPos, lower = 1, upper = genomeAssembly.df$size[which(genomeAssembly.df$chr == plotChr)], len = 1))
    
    assert(checkNull(plotEndPos),
           checkIntegerish(plotEndPos  , lower = 1, upper = genomeAssembly.df$size[which(genomeAssembly.df$chr == plotChr)], len = 1))
    
    assert(checkNull(ylim),
           checkIntegerish(ylim, any.missing = FALSE, len = 2)
    )
    
    allowedPalettes = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
    assertChoice(colorPalette, allowedPalettes)
    
    assertIntegerish(maxWidthLabels, lower = 1, any.missing = FALSE, len = 1)
    
    assertFlag(plotRegionBoundaries)
    assertFlag(plotAllelicBiasResults)
    assertChoice(pValueSummary, c("min", "median"))
    assertFlag(plotGraph)
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
        pdf(fileToPlot)
    }
    
    assertPercentage(signThreshold)
    
    
    if (!testNull(plotChr)) {
        
        if (testNull(plotStartPos)) plotStartPos = 1
        if (testNull(plotEndPos))   plotEndPos = genomeAssembly.df$size[which(genomeAssembly.df$chr == plotChr)]
        
    }
    
    if (!testNull(regions) & !testNull(plotChr)) {
        warning("Both regions and plotChr have been provided. The values for plotChr, plotStartPos and plotEndPos will be set to NULL and the requested SNP regions will be plotted.")
        plotChr = NULL
        plotStartPos = NULL
        plotEndPos = NULL
    }
    
    if (testNull(regions) & testNull(plotChr)) {
        stop("Neither values for regions nor plotChr have been provided.")
    }
    
    
    
    
    
    # Plot a specific user region, this is a user-friendly option to automatically set the location correctly
    if (!is.null(regions)) {
        
        # Check if located on different chromosomes. If yes, abort
        chr = as.character(seqnames(SNPhood.o@annotation$regions[regions]))
        
        if (length(unique(chr)) > 1) {
            stop("The requested SNP regions lie on different chromosomes (",paste0(unique(chr), collapse =","),"). Only one particular chromosome can be plotted at once.")
        }
        
        startPos = as.numeric(start(SNPhood.o@annotation$regions[regions]))
        endPos   = as.numeric(end  (SNPhood.o@annotation$regions[regions]))
        SNPPos   = as.numeric(mcols(SNPhood.o@annotation$regions[regions])$SNPPos)
        
        plotChr = unique(chr)
        plotStartPos = min(startPos)
        plotEndPos = max(endPos)
        
    }
    
    # Determine which SNPs fall into the specified range
    SNPIndex = which(start(annotation(SNPhood.o)$regions) >= plotStartPos & 
                         end(annotation(SNPhood.o)$regions) <= plotEndPos &
                         as.character(seqnames(annotation(SNPhood.o)$regions)) == plotChr)
    
    regions = SNPIndex
    
    # TODO: Avoid redundancy
    startPos = as.numeric(start(SNPhood.o@annotation$regions[SNPIndex]))
    endPos   = as.numeric(end  (SNPhood.o@annotation$regions[SNPIndex]))
    SNPPos   = as.numeric(mcols(SNPhood.o@annotation$regions[SNPIndex])$SNPPos)
    SNPName  = names(SNPhood.o@annotation$regions[SNPIndex])
    
    if (length(SNPIndex) == 0) {
        stop("No SNPs and counts to display between ", plotChr,":", plotStartPos,"-", plotEndPos,"(", plotEndPos - plotStartPos," base pairs). Try to change the location.", sep = "")
        
    }
    
    # Produce the data
    plot.df = data.frame(datasetAndReadGroup = c(), SNPName = c(), SNPPos = c(), value = c(), stringsAsFactors = FALSE)
    
    if (!plotAllelicBiasResults) {
        
        
        for (i in 1:length(readGroups)) {
            
            readGroupCur = readGroups[i]
            
            for (j in 1:length(datasets)) {
                
                datasetCur = datasets[j]
                
                # Construct label name
                if (mergeReadGroupCounts) {
                    
                    label = paste0(readGroups,collapse = ",")
                    if (!testNull(readGroupsOrigLabel)) label = "all" 
                    
                    datasetName = paste0(j,": ",strtrim(label, maxWidthLabels))
                    
                } else {
                    datasetName = paste0(j,": ",strtrim(readGroups[i], maxWidthLabels))
                    
                }
                
                for (k in 1:length(regions)) {
                    
                    regionCur = regions[k]
                    
                    count = SNPhood.o@readCountsUnbinned[[readGroupCur]][[datasetCur]][regionCur]
                    
                    # Add to data frame only if read group counts are not merged
                    addRow = TRUE
                    if (mergeReadGroupCounts) {
                        
                        index = which(plot.df$datasetAndReadGroup == datasetName & 
                                          plot.df$SNPName == SNPName[k] &  
                                          plot.df$SNPPos == SNPPos[k])
                        
                        stopifnot(length(index) <= 1)
                        
                        if (length(index == 1)) {
                            addRow = FALSE
                            
                            plot.df$value[index] = plot.df$value[index] + count 
                        }
                        
                    } 
                    
                    if (addRow) {
                        plot.df = rbind(plot.df, data.frame(datasetAndReadGroup = datasetName, 
                                                            SNPName = SNPName[k],
                                                            SNPPos = SNPPos[k],
                                                            value = count, stringsAsFactors = FALSE))
                    }
                    
                }
                
            }
            
        }
        
    } else  {#  plotAllelicBiasResults == TRUE
        
        for (j in seq_len(length(datasets))) {
            
            datasetName = paste0(datasets[j])
            
            for (k in 1:length(regions)) {
                regionsCur = regions[k]
                
                if (pValueSummary == "min") {
                    val = min(SNPhood.o@additionalResults$allelicBias$pValue[[j]][regionsCur,])
                } else {
                    val = median(SNPhood.o@additionalResults$allelicBias$pValue[[j]][regionsCur,])
                }
                
                
                plot.df = rbind(plot.df, data.frame(datasetAndReadGroup = datasetName, 
                                                    SNPName = SNPName[k],
                                                    SNPPos = SNPPos[k],
                                                    value = val, 
                                                    stringsAsFactors = FALSE))
            }
        }
        
        # Subset with only the significant ones
        plot.df$sign  = plot.df$value <= signThreshold
        nSign = length(which(plot.df$sign == TRUE))
        
        plot.df$value = -log(plot.df$value ,10)  
        pValueSigThresholdTransformed = -log(signThreshold, 10)
        
    }
    
    customYLimits = FALSE
    if (!testNull(ylim)) {
        customYLimits = TRUE
        
        if (ylim[1] > ylim[2]) {
            stop("Value for parameter ylim not correct: The second value must be larger than the first value.")
        }
        
    } else {
        ylim = c(-5, max(plot.df$value) + 5)
        if (plotAllelicBiasResults) {
            ylim = c(-0.5, max(max(plot.df$value), pValueSigThresholdTransformed))
        }
        ylim[2] = ylim[2] + 0.1 * ylim[2]
    }
    
    yLabLabel = paste0(.getBinLabelYAxis(SNPhood.o), "\n")
    if (plotAllelicBiasResults) {
        
        if (pValueSummary == "min") {
            yLabLabel = "Most significant bin per SNP region\n (-log 10 p-value)\n"
        } else {
            yLabLabel = "Median p-value per SNP region\n (-log 10 transformed)\n"
        }
        
    }
    
    
    xLabLabel = paste0("\nPosition on chromosome ", plotChr,"(in bp)")
    
    if (plotRegionLabels) {
        xLabLabel = paste0("\nPosition on chromosome ", plotChr, ":", plotStartPos,"-", plotEndPos)
        
        axis.df = plot.df[,c("SNPName", "SNPPos")]
        axis.df = axis.df[!duplicated(axis.df),]
        axis.df = axis.df[order(axis.df$SNPPos, decreasing = FALSE),]
    }
    
    mainLabel = paste0(.getBinLabelYAxis(SNPhood.o), " for the region\n",plotChr, ":", plotStartPos,"-", plotEndPos," (covering ",length(regions)," SNPs)")
    if (plotAllelicBiasResults) {
        
        labelDataset = ifelse(length(datasets) > 1, "datasets", "dataset")
        mainLabel = paste0("Allelic bias overview (", 
                           SNPhood.o@additionalResults$allelicBias$parameters$readGroupsTested[1], " vs. ",
                           SNPhood.o@additionalResults$allelicBias$parameters$readGroupsTested[2], 
                           ") for the region\n",plotChr, ":", plotStartPos,"-", plotEndPos," across ",
                           length(datasets), " ",labelDataset, " (covering ",length(regions)," SNPs). \n",
                           "Significant results: ",nSign, " out of ", nrow(plot.df))
        
    }
    
    legendTitle = "Dataset: read group"
    if (plotAllelicBiasResults) {
        legendTitle = "Dataset"
    }
    
    if (plotAllelicBiasResults) {
        p <- ggplot(plot.df, aes(x = SNPPos, y = value, colour = datasetAndReadGroup, shape = as.factor(sign))) + xlab(xLabLabel) + ylab(yLabLabel) + geom_point(size = sizePoints)
        
        # handle special cases when all values are significant or non significant
        shapeValues = c(1,19)
        if (all(plot.df$sign))   shapeValues = 19
        if (all(!plot.df$sign))  shapeValues = 1
        p <- p + scale_shape_manual(values = shapeValues)
    } else {
        p <- ggplot(plot.df, aes(x = SNPPos, y = value, colour = datasetAndReadGroup)) + xlab(xLabLabel) + ylab(yLabLabel) + geom_point(size = sizePoints, shape = 19)
    }
    
    
    p <- p + coord_cartesian(ylim = ylim)
    p <- p + ggtitle(mainLabel)
    if (plotRegionBoundaries) p <- p + geom_vline(xintercept = unique(c(startPos, endPos)), size = 0.2, linetype = "dashed")
    if (plotAllelicBiasResults)  p <- p + geom_hline(yintercept = pValueSigThresholdTransformed, linetype = "dotted")
    p <- p + .getThemeForGGPlot(verticalLines = FALSE)
    p <- p + labs(colour = legendTitle)
    p <- p + scale_colour_manual(values = .generateColorsForReadGroupsAndDatasets(nDatasets(SNPhood.o), nReadGroups(SNPhood.o), colorPalette, saturationMin = 0.3, deleteFirst = FALSE, deleteLast = TRUE))
    
    if (plotAllelicBiasResults)  p <- p + labs(shape = "p-value\nsignificant") + scale_fill_discrete(labels = c("no","yes"))
    
    if (plotRegionLabels) {
        p <- p + scale_x_continuous(breaks = axis.df$SNPPos, labels = axis.df$SNPName)
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    
    if (plotGraph) {
        
        if (!testNull(fileToPlot)) {
            print(p)
            dev.off()  
        } 
        return(print(p))
    } else {
        return(p)
    }
    
}



#' Clustering of read counts or enrichmens across bins for a specific dataset and read group
#' 
#' \code{plotAndClusterMatrix} can be used to cluster regions such as SNPs based on their local neighbourhood. 
#' The underlying clustering is done using partitioning around medoids (PAM). For more details, see the vignette.
#' 
#' @template SNPhood
#' @template readGroup
#' @template dataset
#' @template nClusters
#' @template normalize
#' @param clustersToPlot Integer. Default NULL. Vector of clusters that should be plotted. If set to NULL, all clusters from the clustering result
#' will be plotted. Otherwise, only the clusters as specified by the user are plotted, omitting regions belonging to other clusters. This
#' is useful to, for example, only display regions that show a bin-dependent pattern and are not invariant across the whole user region.
#' @template fileToPlot
#' @template verbose
#' @param ... Additional graphical parameters that can be used to modify the output of the function levelplot (panel.levelplot). 
#' See ?levelplot for details.
#' @return The clustering reports the cluster in which each SNP falls, the average silhouette for pam clustering, plots for the clustered reads and a summary plot of average reads per cluster across the region being analyzed. 
#' @export
#' @examples 
#' data(SNPhood.o, package="SNPhood")
#' SNPhood.o = plotAndClusterMatrix(SNPhood.o, readGroup = "paternal", dataset = 1, nClustersVec = c(3:6))
#' SNPhood.o = plotAndClusterMatrix(SNPhood.o, readGroup = "paternal", dataset = 1, normalize = FALSE)
#' @import checkmate
#' @importFrom lattice levelplot

plotAndClusterMatrix <- function(SNPhood.o, readGroup, dataset, nClustersVec = 3, normalize = TRUE, clustersToPlot = NULL, fileToPlot = NULL, verbose = TRUE, ...) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    assertIntegerish(nClustersVec, any.missing = FALSE, min.len = 1, unique = TRUE, lower = 2)
    assertChoice(readGroup, annotationReadGroups(SNPhood.o))   
    
    dataset = .checkAndConvertDatasetArgument(SNPhood.o, dataset, nullAllowed = FALSE, maxLength = 1) 
    assertFlag(normalize)
    
    
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
    }
    
    
    if (SNPhood.o@internal$calcEnrichment & !parameters(SNPhood.o)$normByInput) {
        stop("Cannot do clustering on an enrichment matrix because input normalization has not been performed.")
    }
    
    type = ifelse(SNPhood.o@internal$calcEnrichment, "enrichmentBinned", "binned")
    
    target.m = counts(SNPhood.o, type = type, readGroup = readGroup, dataset = dataset)
    
    if (normalize) {
        target.m = .normalizeMatrixForClustering(target.m, verbose)
    }
    
    dimnames(target.m) = list(annotationRegions(SNPhood.o), annotationBins(SNPhood.o))
    
    
    # Init if never done before
    if (testNull(SNPhood.o@additionalResults$clustering)) {
        SNPhood.o@additionalResults$clustering = list()
        
        for (readGroupCur in annotationReadGroups(SNPhood.o)) {
            SNPhood.o@additionalResults$clustering[[readGroupCur]] = list()
        }
        SNPhood.o@internal$addResultsElementsAdded = c(SNPhood.o@internal$addResultsElementsAdded, "clustering")
    }
    
    if (testNull(SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]])) {
        SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]] = list()
    }
    
    if (!is.null(fileToPlot)) pdf(fileToPlot)
    
    for (nClusters in nClustersVec) {
        
        resultsMissing = ifelse(testNull(SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]][[paste0("nClusters",nClusters)]]), TRUE, FALSE)
        if (resultsMissing) {
            if (verbose) message("Performing clustering for ", nClusters, " clusters.")
            SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]][[paste0("nClusters",nClusters)]] = .pamClustering(SNPhood.o, target.m, nClusters, ...)
            
        }
        
        # Plot, but only selected clusters
        clusterMatrix = SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]][[paste0("nClusters",nClusters)]]$clusteringMatrix[[1]]
        assert(checkNull(clustersToPlot), checkSubset(clustersToPlot, unique(clusterMatrix$cluster)))
        p = .generateClusterPlot(clusterMatrix, SNPhood.o, clustersToPlot = clustersToPlot)
        print(p)
        
    }
    
    if (!is.null(fileToPlot)) dev.off()
    

    SNPhood.o
    
}


#' Visualize average enrichment per cluster
#' 
#' \code{plotClusterAverage} visualizes the average reads per cluster. Note that the function \code{plotAndClusterMatrix} has to be executed
#' before \code{plotClusterAverage} is called for the same read group and dataset
#' @template SNPhood
#' @template readGroup
#' @template dataset
#' @template fileToPlot
#' @template verbose
#' @seealso \code{plotAndClusterMatrix}
#' @export 
#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' plot = plotClusterAverage(SNPhood.o, readGroup = "paternal", dataset = 1)

#' @import checkmate
plotClusterAverage <- function(SNPhood.o, readGroup, dataset, fileToPlot = NULL, verbose = TRUE) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    assertChoice(readGroup, annotationReadGroups(SNPhood.o))  
    
    dataset = .checkAndConvertDatasetArgument(SNPhood.o, dataset, nullAllowed = FALSE, maxLength = 1) 
    
    # Retrieve the clustering results
    if (testNull(SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]])) {
        stop("Could not find clustering results for read group ", readGroup, " and dataset ", dataset,". Has the function plotAndClusterMatrix be called before?")
    } 
    
    clusteringResults = SNPhood.o@additionalResults$clustering[[readGroup]][[dataset]]
    
    nClusteringResults = length(clusteringResults)
    
    if (nClusteringResults > 1 & testNull(fileToPlot)) {
        warning("Multiple clustering results found, summarizing all of them. Multiple plots will be produced. Specify a filename for the parameter fileToPlot to see them all.")
    }    
    
    plots.l = list()
    
    if (!testNull(fileToPlot)) pdf(fileToPlot)
    
    for (i in seq_len(nClusteringResults)) {
        plots.l[[i]] = .plotClusterAverage(SNPhood.o, clusteringResults[[i]], fileToPlot = fileToPlot)
    }
    
    if (!testNull(fileToPlot)) dev.off()  
    
    
    plots.l
}

#' @import checkmate
#' @importFrom reshape2 melt
#' @importFrom grDevices pdf dev.off
.plotClusterAverage <- function(SNPhood.o, clusteringResults, fileToPlot = NULL) {
    
    .checkObjectValidity(SNPhood.o)
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
    }
    
    PamObj = clusteringResults
    assertList(PamObj, any.missing = FALSE, min.len = 1)
    assertSubset(c("clusteringMatrix","plots"), names(PamObj))
    
    pam_cluster <- PamObj$clusteringMatrix[[1]]
    stopifnot(ncol(pam_cluster)>2)
    splitting_pam <- split(pam_cluster, pam_cluster$cluster)
    colsTake <- grep("bin", colnames(pam_cluster))
    split_pam <- lapply(splitting_pam, function(x) x[colsTake])
    split_pam <- lapply(split_pam, function(x) colMeans(x))
    split_pam_melt  <- reshape2::melt(split_pam)
    
    allBins<-SNPhood.o@internal$plot_labelBins
    
    takethis <- rep(allBins, length(unique(split_pam_melt$L1)))
    bins <- rep(1:length(colsTake), length(unique(split_pam_melt$L1)))
    split_pam_melt <- cbind(split_pam_melt, bins)
    split_pam_melt <- cbind(split_pam_melt, takethis)
    colnames(split_pam_melt) <- c("value","Cluster","bins","Coord")
    
    ## tabulating frequencies :
    table_frequencies <- reshape2::melt(table(pam_cluster$cluster))
    
    freqLabel <- c()
    for (i in 1:nrow(table_frequencies)) {
        temp <- paste0(table_frequencies[i,1]," - ", paste(table_frequencies[i,2])," regions")
        freqLabel <- c(freqLabel, temp)
    }
    split_pam_melt$Cluster <- factor(split_pam_melt$Cluster)
    levels(split_pam_melt$Cluster) <- freqLabel
    
    
    posVerticalLinewhich <- which(allBins == SNPhood.o@internal$plot_origBinSNPPosition)
    
    
    p <- ggplot(split_pam_melt, aes(x = bins, y = value)) + geom_line(aes(col = Cluster, group = Cluster))
    p <- p + xlab(.getBinLabelXAxis(SNPhood.o)) + ylab("Average enrichment per cluster") + theme_bw()
    p <- p + .getBinAxisLabelsForGGPlot(SNPhood.o@internal$plot_labelBins)
    p <- p + .getVerticalLineForGGPlot(SNPhood.o@internal$plot_origBinSNPPosition)
    p <- p + .getThemeForGGPlot()
    
    
    print(p)
    
    
    
    return(p)
    
}


#' Visualize average counts/enrichment based on strong and weak genotypes.
#' 
#' The function \code{plotGenotypesPerCluster} plots average clusters per genotype based on the clustering results of
#' the strong an weak genotype analysis (see \code{\link{plotAndCalculateWeakAndStrongGenotype}}), which has to be executed before. 
#' @template SNPhood 
#' @param printBinLabels Logical(1). Default TRUE. Should the bin labels be printed? 
#' If multiple clusters are plotted simultaenously, bin labels might overlap, in which case \code{printBinLabels} can be set to FALSE.
#' @template fileToPlot
#' @param printPlot Logical(1). Default TRUE. Should the plots be printed? Only relevant if \code{fileToPlot} is set to NULL; otherwise, the plots
#' are always printed to the output file.
#' @template ggplotReturn
#' @export
#' @seealso \code{\link{plotAndCalculateWeakAndStrongGenotype}}
#' @examples 
#' data(SNPhood.o, package="SNPhood")
#' SNPhood_merged.o = mergeReadGroups(SNPhood.o)
#' SNPhood_merged.o = plotAndCalculateWeakAndStrongGenotype(SNPhood_merged.o)
#' plot = plotGenotypesPerCluster(SNPhood_merged.o, printPlot = FALSE)

#' @importFrom ggplot2 ggplot scale_x_discrete aes geom_line facet_grid theme_bw geom_vline element_blank ggtitle
#' @importFrom grDevices pdf dev.off
#' @importFrom reshape2 melt
#' @import checkmate
plotGenotypesPerCluster = function(SNPhood.o, printBinLabels = TRUE, fileToPlot = NULL, printPlot = TRUE) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    if (nReadGroups(SNPhood.o) > 1) {
        stop(.getErrorMessageReadGroupSpecificty())
    }
    
    if (testNull(SNPhood.o@additionalResults$genotype)) {
        stop("Could not find genotype and cluster results. Has the function plotAndCalculateWeakAndStrongGenotype ba called before?")
    }
    
#     allGenotypes = c("weak", "strong")
#     assert(checkNull(type), checkChoice(type, allGenotypes))
#     
#     if (testNull(type)) {
#         type = allGenotypes
#     }
    
    assertSubset(names(SNPhood.o@additionalResults$genotype),choices = c("strongGenotypes","weakGenotypes","invariantGenotypes","clustering"))
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
    }
    
    plotsGenotypes = vector("list",length(type))
    
    plotsGenotypes = lapply(plotsGenotypes, 
                                  function(x) vector("list",length(SNPhood.o@additionalResults$genotype$clustering$strongGenotypes[[1]]$clusteringMatrix)))
 
    
    header = c("strong","weak")
    

    nPlots = 1
    
    for (j in 1:length(plotsGenotypes)) {
        
        for (i in 1:length(SNPhood.o@additionalResults$genotype$clustering[[j]])) {
 
            nPlots = nPlots + 1
            clusters = SNPhood.o@additionalResults$genotype$clustering[[j]][[i]]$clusteringMatrix[[1]]$cluster
            GenotypeClusters = cbind(SNPhood.o@additionalResults$genotype[[j]], clusters)
            GenotypeClusters_byCluster = split(GenotypeClusters,GenotypeClusters$clusters)
            GenotypeClusters_byClusterGenotype = lapply(GenotypeClusters_byCluster, function(x) split(x,x$genotype))
            averagePerclusterPerGenotype =  lapply(GenotypeClusters_byClusterGenotype, function(x) lapply(x,function(y) colMeans(y[,1:nBins(SNPhood.o)])))
            melt_Genotypes = reshape2::melt(averagePerclusterPerGenotype)
            index = rep(1:nBins(SNPhood.o),length(GenotypeClusters_byClusterGenotype) * 3)
            melt_Genotypes = cbind(melt_Genotypes,index)  
            linePos = SNPhood.o@internal$plot_origBinSNPPosition
            colnames(melt_Genotypes) = c("value","Genotype","L1", "index")
            melt_Genotypes$L1 = paste("Cluster",melt_Genotypes$L1)
            
            p = ggplot(melt_Genotypes,aes(x=index,y=value)) + geom_line(aes(col=Genotype, group=Genotype))
            p <- p + facet_grid(.~L1) + theme_bw() + geom_vline(xintercept=linePos, linetype = "longdash")
           
            p <- p + ylab("Average Reads") 
            if (printBinLabels) {
                p <- p + scale_x_discrete(labels = SNPhood.o@internal$plot_labelBins)
            } else {
                p <- p + scale_x_discrete(labels = element_blank())
            }
            p <- p + xlab(.getBinLabelXAxis(SNPhood.o = SNPhood.o))
            p <- p + ggtitle(paste("Clusters of", header[j], "genotypes ( Number of clusters: ",length(unique(clusters)),")"))
            
            plotsGenotypes[[j]][[i]] = p
        }
        
        
    }
    
    if (!testNull(fileToPlot)) pdf(fileToPlot)
    
    if (testNull(fileToPlot) & nPlots > 1 & printPlot) {
        warning("Multiple plots have been produced but no output file has been given. Only the last plot may be visible.")
    }
    
    for (j in 1:length(plotsGenotypes)) {
        
        for (i in 1:length(SNPhood.o@additionalResults$genotype$clustering[[j]])) {
            if (printPlot | !testNull(fileToPlot)) print(plotsGenotypes[[j]][[i]])
        }
    }
    
    if (!testNull(fileToPlot)) dev.off()
    

    plotsGenotypes
}



#' Plot genotype frequencies of regions across datasets. 
#' 
#' Creates bar plots for the distribution of genotype frequencies of regions across individuals.

#' @template SNPhood
#' @template regions
#' @template fileToPlot

#' @template ggplotReturn
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' plot = plotGenotypesPerSNP(SNPhood.o, regions=1:20)
#' @seealso \code{\link{plotAndClusterMatrix}}
#' @export
#' @import checkmate 
#' @importFrom ggplot2 ggplot aes geom_bar ylab scale_y_discrete theme element_text
#' @importFrom reshape2 melt
#' @importFrom grDevices pdf dev.off
plotGenotypesPerSNP <- function(SNPhood.o, regions = NULL, fileToPlot = NULL) {
    
    .checkObjectValidity(SNPhood.o)
    
    if (SNPhood.o@config$onlyPrepareForDatasetCorrelation) {
        stop(.getErrorForOnlyPrepareSamplesCorrelation())
    }
    
    if (testNull(regions)) {
        regions = annotationRegions(SNPhood.o)
    }
    
    regions = .checkAndConvertRegionArgument(SNPhood.o, regions, nullAllowed = TRUE, maxLength = NULL) 
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    
    if (testNull(annotation(SNPhood.o)$genotype$external)) {
        stop("Could not find genotype information in the object. Has the genotype been integrated?")
    }
    
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
        pdf(fileToPlot)
    }
    
    
    genotype = annotation(SNPhood.o)$genotype$external[regions, ]
    genotype<-genotype[,4:ncol(genotype)]
    rownames(genotype) <- annotationRegions(SNPhood.o)[regions]
    
    rowsToDelete = which(rowSums(is.na(genotype)) == ncol(genotype))
    if (length(rowsToDelete) > 0) {
        genotypes_individuals_refined <- genotype[-rowsToDelete, 
                                                  ]
        message("Cannot display ", length(rowsToDelete), " SNP genotypes because all genotypes were NA. Remaining number of SNPs to display: ", 
                nrow(genotypes_individuals_refined))
    } else {
        genotypes_individuals_refined <- genotype
    }
    
    
    individualsNA <- apply(genotypes_individuals_refined, 2, function(x) !all(is.na(x)))
    Nas <- which(!individualsNA)
    if (length(Nas) > 0) {
        genotypes_individuals_refined <- genotypes_individuals_refined[,-Nas]
        message(paste("Individual", Nas," does not have any genotype information and will be excluded \n"))
    }
    
    for (i in 1:ncol(genotypes_individuals_refined)){
        splitGenotypes <-strsplit(x = as.character(genotypes_individuals_refined[,i]),"[^0-9]+")
        splitGenotypes <-lapply(splitGenotypes,as.numeric)
        genotypes_individuals_refined[,i]<-unlist(lapply(splitGenotypes,sum))
        
        
    }
    
    ones <- c()
    twos <- c()
    zeros <- c()
    for (i in 1:nrow(genotypes_individuals_refined)) {
        ones.temp <- length(which(genotypes_individuals_refined[i, 
                                                                ] == 1))
        ones <- c(ones, ones.temp)
        twos.temp <- length(which(genotypes_individuals_refined[i, 
                                                                ] == 2))
        twos <- c(twos, twos.temp)
        zeros.temp <- length(which(genotypes_individuals_refined[i, 
                                                                 ] == 0))
        zeros <- c(zeros, zeros.temp)
    }
    
    
    sumEachGenotype <- cbind(zeros, ones, twos)
    rownames(sumEachGenotype) <- rownames(genotypes_individuals_refined)
    sumEachGenotype_melt <- reshape2::melt(t(sumEachGenotype))
    colnames(sumEachGenotype_melt) <- c("Genotype", "SNP", "value")
    sumEachGenotype_melt$SNP = as.character(sumEachGenotype_melt$SNP)
    sumEachGenotype_melt$Genotype <- as.factor(sumEachGenotype_melt$Genotype)
    levels(sumEachGenotype_melt$Genotype) <- c("1", "2", "0")
    
    
    characCount <- nchar(sumEachGenotype_melt$SNP)
    longChars <- which(characCount > 9)
    sumEachGenotype_melt$SNP[longChars] <- paste(substr(sumEachGenotype_melt$SNP[longChars], 1, 9), ".", sep = "")
    
    sumEachGenotype_melt$Genotype <- factor(sumEachGenotype_melt$Genotype, levels = c(0,1,2))
    
    p <- ggplot(sumEachGenotype_melt, aes(x = SNP, y = value, fill = Genotype)) + geom_bar(position = "dodge", stat = "identity") 
    p <- p + ylab("Number of datasets per genotype")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    
    if (!testNull(fileToPlot)) {
        
        print(p)
        dev.off()
        
    } else {
        print(p)
    }
    
    return(p)
}


#' Visualizes and calculates strong and weak genotypes.
#' 
#' The function \code{plotAndCalculateWeakAndStrongGenotype} finds the strongest and weakest genotypes based on reads extracted around each region. Strong and weak genotypes are found using the reads extracted from SNPhood and their corresponding genotypes as found by the function \code{associateGenotypes}
#' Note the reads have to be merged using the function \code{mergeReadGroups} before running this function. 
#' @template SNPhood
#' @template normalize
#' @template nClusters
#' @template fileToPlot
#' @template verbose
#' @return Modified \code{\linkS4class{SNPhood}} object with the results of the analysis stored in the object. 
#' Specifically, a matrix for average reads per SNP for datasets which have strong and weak genotypes, respectively, are stored in the
#' slot additionalResults$genotype. 
#' The SNPs which have invariant genotypes across all the samples being analyzed are also saved. 
#' In addition, clustering on the strong and weak genotype read mateices are reportd as in the function \code{\link{plotAndClusterMatrix}}.
#' @export
#' @examples
#' data(SNPhood.o, package="SNPhood")
#' SNPhood_merged.o = mergeReadGroups(SNPhood.o)
#' SNPhood_merged.o = plotAndCalculateWeakAndStrongGenotype(SNPhood_merged.o, nClustersVec = 6)
#' SNPhood_merged.o = plotAndCalculateWeakAndStrongGenotype(SNPhood_merged.o, nClustersVec = 2:6, verbose = FALSE)
#' @import checkmate
plotAndCalculateWeakAndStrongGenotype <- function(SNPhood.o, normalize = TRUE, nClustersVec = 3, fileToPlot = NULL, verbose = FALSE) {
    
    .checkObjectValidity(SNPhood.o)
    assertFlag(normalize)
    assertIntegerish(nClustersVec, any.missing = FALSE, min.len = 1, unique = TRUE, lower = 2)
    assertFlag(verbose)
    
    assert(checkNull(fileToPlot), checkCharacter(fileToPlot, min.chars = 1, len = 1))
    if (!testNull(fileToPlot)) {
        assertDirectory(dirname(fileToPlot), access = "r")
    }
    
    if (nReadGroups(SNPhood.o) > 1) {
        stop(.getErrorMessageReadGroupSpecificty())
    }
    
    
    countsObj <- counts(SNPhood.o, type = "binned")$allReadGroups
    genotypeAssociated <- annotation(SNPhood.o)$genotype$external
    rownames(genotypeAssociated) <- rownames(annotation(SNPhood.o)$genotype$external)
    finalCounts <- lapply(countsObj, function(x) x[match(rownames(genotypeAssociated), 
                                                         annotationRegions(SNPhood.o)), ])
    
    if (ncol(genotypeAssociated) < 4) {
        stop("Cannot summarize genotypes because no genotypoes have been added to slot annotation$genotype$external. Has the function associateGenotypes be executed?")
    }
    
    
    onlyGenotypes <- genotypeAssociated[,4:length(genotypeAssociated)] # Getting the columns with only genotypes
    onlyGenotypes <- as.matrix(onlyGenotypes) 
    
    for (i in 1:ncol(onlyGenotypes)) {
        if (all(is.na(onlyGenotypes[,i]))) { #identifying if genotype information for inividual is NA
            onlyGenotypes[,i] <- NA
        } else {# Split and designate 0,1,2
            splitGenotypes <- strsplit(x = onlyGenotypes[,i],"[^0-9]+")
            splitGenotypes <- lapply(splitGenotypes,as.numeric)
            onlyGenotypes[,i] <- unlist(lapply(splitGenotypes,sum))
            
        }
    }
    
    filterdGenotypes <- as.matrix(onlyGenotypes[, apply(onlyGenotypes, 2, function(x) !all(is.na(x)))]) # getting columns (individuals) for whom some genotype info is available
    filterdGenotypes <- as.matrix(filterdGenotypes)
    if (ncol(filterdGenotypes) == 1) {
        stop("The comparison requrires genotype information atleast from 2 individuals. Only one of the individuals has genotype information for the regions being analysed.")
    }else{
        if (verbose) message("Identifying strong and weak genotypes.")
    }
    
    
    invariantGenotypes = which(apply(filterdGenotypes,1,function(x) length(unique(x))) == 1)
    invariantGenotypes = apply(filterdGenotypes[invariantGenotypes,],1,unique)
    naGenotypes = apply(filterdGenotypes,1,function(x) all(!is.na(x)))
    filterdGenotypes = filterdGenotypes[naGenotypes,]
    
    SNPhood.o@readCountsBinned[[annotationReadGroups(SNPhood.o)]] <- lapply(SNPhood.o@readCountsBinned[[annotationReadGroups(SNPhood.o)]], 
                                                                            function(x){ row.names(x)<-rownames(onlyGenotypes); x})
    
    refinedIndividualFiles <- lapply(SNPhood.o@readCountsBinned[[annotationReadGroups(SNPhood.o)]], 
                                     function(x) x[match(row.names(filterdGenotypes), row.names(x)), ])
    refinedIndividualFiles <- refinedIndividualFiles[apply(onlyGenotypes,2,function(x) !all(is.na(x)))]
    
    if (normalize) {
        for (i in 1:length(refinedIndividualFiles)) {
            refinedIndividualFiles[[i]] = .normalizeMatrixForClustering(refinedIndividualFiles[[i]], verbose)
        }
    }
    
    reads = lapply(refinedIndividualFiles,function(x) rowSums(x))
    reads = do.call(cbind,reads)
    genotype_strong <- c()
    genotype_temp <- vector("list", 3)
    chooseGenotype <- c()
    zeros <- c()
    ones <- c()
    twos <- c()
    strong_genotypes <- matrix(data = ,nrow = nrow(filterdGenotypes),ncol = nBins(SNPhood.o))
    
    weak_genotypes <- matrix(data = ,nrow = nrow(filterdGenotypes),ncol = nBins(SNPhood.o))
    
    genotype_weak <- c()
    for (i in 1:nrow(filterdGenotypes)) {
        genotype_temp[[1]] <- which(filterdGenotypes[i, ] == 0)
        genotype_temp[[2]] <- which(filterdGenotypes[i, ] == 1)
        genotype_temp[[3]] <- which(filterdGenotypes[i, ] == 2)
        
        zeros <- sum(reads[i, which(filterdGenotypes[i, ] == 0)])/length(genotype_temp[[1]])
        ones <- sum(reads[i, which(filterdGenotypes[i, ] == 1)])/length(genotype_temp[[2]])
        twos <- sum(reads[i, which(filterdGenotypes[i, ] == 2)])/length(genotype_temp[[3]])
        compare_genotype <- c(zeros, ones, twos)
        max_genotype <- which.max(compare_genotype)
        min_genotype <- which.min(compare_genotype)
        genotype_strong <- c(genotype_strong, max_genotype)
        genotype_weak <- c(genotype_weak, min_genotype)
        chooseGenotype <- genotype_temp[[max_genotype]]
        chooseGenotypeWeak <- genotype_temp[[min_genotype]]
        strong_genotypes[i,] = colMeans(do.call(rbind,lapply(refinedIndividualFiles[chooseGenotype], function(x) x[i,])))
        weak_genotypes[i,] = colMeans(do.call(rbind,lapply(refinedIndividualFiles[chooseGenotypeWeak], function(x) x[i,])))
    }
    rownames(strong_genotypes) <- row.names(refinedIndividualFiles[[1]])
    rownames(weak_genotypes) <- row.names(refinedIndividualFiles[[1]])
    strong_genotypes <- cbind(strong_genotypes, genotype_strong)
    weak_genotypes <- cbind(weak_genotypes, genotype_weak)
    weak_genotypes <- as.data.frame(weak_genotypes, stringsAsFactors = FALSE)
    strong_genotypes <- as.data.frame(strong_genotypes, stringsAsFactors = FALSE)
    strong_genotypes$genotype_strong <- factor(strong_genotypes$genotype_strong)
    levels(strong_genotypes$genotype_strong) <- c(0, 1, 2)
    weak_genotypes$genotype_weak <- factor(weak_genotypes$genotype_weak)
    levels(weak_genotypes$genotype_weak) <- c(0, 1, 2)
    binNames <- paste0("bin_",1:nBins(SNPhood.o))
    colnamesGenotypes = c(binNames,"genotype")
    colnames(strong_genotypes) = colnamesGenotypes
    colnames(weak_genotypes) = colnamesGenotypes
    
    SNPhood.o@additionalResults$genotype$strongGenotypes <- strong_genotypes
    SNPhood.o@additionalResults$genotype$weakGenotypes   <- weak_genotypes
    SNPhood.o@additionalResults$genotype$invariantGenotypes <- invariantGenotypes
    
    if (!is.null(fileToPlot)) pdf(fileToPlot)
    
    for (nClusters in nClustersVec) {
        if (verbose) message("Performing clustering for ", nClusters, " clusters.")
        SNPhood.o@additionalResults$genotype$clustering$strongGenotypes[[paste0("nClusters",nClusters)]] = .pamClustering(SNPhood.o, as.matrix(strong_genotypes[,1:nBins(SNPhood.o)]), nClusters)
        SNPhood.o@additionalResults$genotype$clustering$weakGenotypes[[paste0("nClusters",nClusters)]]   = .pamClustering(SNPhood.o, as.matrix(weak_genotypes  [,1:nBins(SNPhood.o)]), nClusters)
    }
    
    if (!is.null(fileToPlot)) dev.off()
    
    SNPhood.o
}



#' @importFrom ggplot2 theme_bw theme element_blank
.getThemeForGGPlot <- function(verticalLines = TRUE) {
    
    if (verticalLines) {
        theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank()) + theme_bw()
    } else {
        theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(),
              panel.grid.minor.x=element_blank(),panel.grid.major.x=element_blank()) + theme_bw()
    }
    
}

#' @importFrom ggplot2 geom_vline
.getVerticalLineForGGPlot <- function(xPos) {
    
    geom_vline(xintercept = xPos, linetype = "longdash", size = 0.5)
}

#' @importFrom ggplot2 scale_x_discrete
.getBinAxisLabelsForGGPlot <- function(binPos) {
    
    scale_x_discrete(labels = binPos)
}

#' @import checkmate
.produceTitleForPlot <- function(SNPhood.o, region) {
    
    .checkObjectValidity(SNPhood.o)
    assertIntegerish(region, lower = 1, upper = nRegions(SNPhood.o), len = 1)
    
    plotChr          = as.character(seqnames(SNPhood.o@annotation$regions)  [region]) 
    plotStartPos     = as.character(start  (SNPhood.o@annotation$regions)   [region]) 
    plotEndPos       = as.character(end    (SNPhood.o@annotation$regions)   [region])  
    mainLabel = paste0(annotationRegions(SNPhood.o)[region],
                       " (", plotChr, ":",  plotStartPos, "-", plotEndPos, ")"
    )
    
    mainLabel
    
}

#' @import checkmate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette 
.generateColorsForReadGroupsAndDatasets <- function(nDatasets, nReadGroups, paletteName, saturationMin = 0.3, deleteFirst = FALSE, deleteLast = TRUE) {
    
    allowedPalettesMaxColors = c(8, 8, 12, 9, 8, 9, 8, 12)
    names(allowedPalettesMaxColors) = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
    
    assertIntegerish(nDatasets, any.missing = FALSE, len = 1, lower = 1)
    assertIntegerish(nReadGroups, any.missing = FALSE, len = 1, lower = 1)
    assertChoice(paletteName, names(allowedPalettesMaxColors))
    assertPercentage(saturationMin)
    assertFlag(deleteFirst)
    assertFlag(deleteLast)
    
    if (nDatasets <= allowedPalettesMaxColors[paletteName]) {
        
        mainColors = brewer.pal(allowedPalettesMaxColors[paletteName],paletteName)[1:nDatasets]
        
    } else {
        # Delete the last color, as this is sometimes gray and does not work with the desaturation
        nColorsIncrease = as.integer(deleteFirst) + as.integer(deleteLast)
        mainColors = colorRampPalette(brewer.pal(allowedPalettesMaxColors[paletteName],paletteName))(nDatasets + nColorsIncrease)
        if (deleteFirst)  mainColors = mainColors[-1]
        if (deleteLast)   mainColors = mainColors[-length(mainColors)]
        
        
        
    }
    
    
    
    
    saturationMax = 1
    saturationDiffStep = (saturationMax - saturationMin) / (nReadGroups - 1)
    saturation.vec = seq(saturationMax, saturationMin, -saturationDiffStep)
    
    col.vec = c()
    for (i in 1:nReadGroups) {
        for (j in 1:nDatasets) {
            col.vec = c(col.vec, .desat(cols = mainColors[j], sat = saturation.vec[i]))
        }
    }
    
    col.vec
}

#' @importFrom grDevices rgb2hsv col2rgb hsv
.desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

.getBinLabelXAxis <- function(SNPhood.o) {
    
    paste0("\nDistance from SNP in ",SNPhood.o@config$binSize," bp bins")
}

.getBinLabelYAxis <- function(SNPhood.o) {
    
    label = "Read counts"
    
    if (SNPhood.o@internal$countType == "readCountsNormalized") {
        label = "Read counts (normalized)"
    } else if (SNPhood.o@internal$countType == "enrichment") {
        label = "Enrichment"
    } 
    
    label
    
}

.getErrorMessageReadGroupSpecificty <- function() {
    
    paste0("Read counts are stored specifically for each read group in the SNPhood object. This function requires non-allele-specific reads, however. Run the function mergeReadGroups first and try again.")
    
}

