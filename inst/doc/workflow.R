## ----<knitr, echo=FALSE, message=FALSE, results="hide"-------------------
library("knitr")
opts_chunk$set(
  tidy = TRUE,
  dev = "png",
  fig.width = 11,
  cache = TRUE,
  message = FALSE)

## ----loadLibaries, echo=TRUE, message=FALSE, results="hide"--------------
library(SNPhoodData)
library(SNPhood)

## ----checkFiles, echo=TRUE-----------------------------------------------
(files = list.files(pattern = "*",system.file("extdata", package = "SNPhoodData"),full.names = TRUE))
fileUserRegions = files[grep(".txt", files)]
fileGenotypes   = files[grep("genotypes", files)]

## ----getParameters, echo=TRUE--------------------------------------------
(par.l = getDefaultParameterList(path_userRegions = fileUserRegions))

## ----changeParameters, echo=TRUE-----------------------------------------
# Verify that you do not have zero-based coordinates
par.l$zeroBasedCoordinates
par.l$readGroupSpecific
par.l$poolDatasets = TRUE
par.l$binSize = 25

## ----createFileList, echo=TRUE-------------------------------------------
patternBAMFiles = paste0(dirname(files[3]), "/*.bam")
(files.df = collectFiles(patternFiles = patternBAMFiles, verbose = TRUE))

## ----assignIndividualID, echo=TRUE---------------------------------------
files.df$individual = c("GM10847", "GM10847", "GM12890", "GM12890")
files.df

## ----qualityTestPrep, echo=TRUE, results="hide"--------------------------
par.l$poolDatasets = FALSE
SNPhood.o = analyzeSNPhood(par.l, files.df, TRUE, verbose = TRUE)

## ----qualityTest, echo=TRUE----------------------------------------------
SNPhood.o = plotAndCalculateCorrelationDatasets(SNPhood.o, fileToPlot = NULL)
corrResults = results(SNPhood.o, type = "samplesCorrelation")
mean(corrResults$corTable[lower.tri(corrResults$corTable)])

## ----runAnalysis, echo=TRUE, results="hide"------------------------------
par.l$poolDatasets = TRUE
SNPhood.o = analyzeSNPhood(par.l, files.df, onlyPrepareForDatasetCorrelation = FALSE, verbose = TRUE)

## ----SNPhoodObj, echo=TRUE, results="hide"-------------------------------
SNPhood.o

## Retrieve number of regions, datasets, bins, and read groups
nRegions(SNPhood.o)
nDatasets(SNPhood.o)
nBins(SNPhood.o)
nReadGroups(SNPhood.o)

## Retrieve general annotation of SNPhood object with all its different elements
names(annotation(SNPhood.o))
annotation(SNPhood.o)$regions
annotation(SNPhood.o)$files

## Retrieve the parameters that were used for the analysis
head(parameters(SNPhood.o))
names(parameters(SNPhood.o))

## Retrieve annotation of regions
head(annotationRegions(SNPhood.o))

## Retrieve annotation of bins
annotationBins(SNPhood.o)
head(annotationBins2(SNPhood.o,fullAnnotation = TRUE))
SNP_names = c("rs7275860","rs76473124")
head(annotationBins2(SNPhood.o, regions = 1:10, fullAnnotation = FALSE))
annotationBins2(SNPhood.o, regions = SNP_names, fullAnnotation = TRUE)

## Retrieve annotation of datasets
annotationDatasets(SNPhood.o)
annotationReadGroups(SNPhood.o)

## Extract counts after the binning

# Extract one count matrix from the paternal read group from the first dataset
head(counts(SNPhood.o, type = "binned", readGroup = "paternal", dataset = 1))

# Extract count matrices from all read groups from the first dataset
str(counts(SNPhood.o, type = "binned", readGroup = NULL, dataset = 1))

# Extract count matrices from all read groups from the first dataset (using its name)
DataSetName <- annotationDatasets(SNPhood.o)[1]
str(counts(SNPhood.o, type = "binned", readGroup = NULL, dataset = DataSetName))

# Extract count matrices from all read groups from the all dataset
str(counts(SNPhood.o, type = "binned", dataset = NULL))

## Similarly, you can also extract counts before the binning
head(counts(SNPhood.o, type = "unbinned", readGroup = "paternal", dataset = 1))

## If you had enrichments instead of counts, you would employ the enrichments method in analogy to counts
enrichment(SNPhood.o, readGroup = "paternal")


## ----SNPhoodObj2, echo=TRUE, results="hide"------------------------------
SNPhood.o

## Rename regions, datasets, bins, and read groups
mapping = as.list(paste0(annotationRegions(SNPhood.o),".newName"))
names(mapping) = annotationRegions(SNPhood.o)
SNPhood_mod.o = renameRegions(SNPhood.o, mapping)

mapping = list("Individual1", "Individual2")
names(mapping) = annotationDatasets(SNPhood.o)
SNPhood_mod.o = renameDatasets(SNPhood.o, mapping)

mapping = list("Bin1_NEW")
names(mapping) = annotationBins(SNPhood.o)[1]
SNPhood_mod.o = renameBins(SNPhood.o, mapping)

mapping = list("a", "b", "c")
names(mapping) = annotationReadGroups(SNPhood.o)
SNPhood_mod.o = renameReadGroups(SNPhood.o, mapping)

## Delete regions, datasets, and read groups (deleting bins is still in development)
SNPhood_mod.o = deleteRegions(SNPhood.o, regions = 1:5)
SNPhood_mod.o = deleteRegions(SNPhood.o, regions = c("rs9984805", "rs59121565"))

SNPhood_mod.o = deleteDatasets(SNPhood.o, datasets = 1)
SNPhood_mod.o = deleteDatasets(SNPhood.o, datasets = "GM12890")

# For read groups, we currently support only a name referral
SNPhood_mod.o = deleteReadGroups(SNPhood.o, readGroups = "paternal")

## Merge read groups
SNPhood_merged.o = mergeReadGroups(SNPhood.o)
nReadGroups(SNPhood_merged.o)
annotationReadGroups(SNPhood_merged.o)

## ----visualizeCounts, echo=TRUE, results="hide"--------------------------

plotBinCounts(SNPhood.o, region = 2)
plotBinCounts(SNPhood.o, region = 2, plotGenotypeRatio = TRUE, readGroups = c("paternal","maternal"))

plotRegionCounts(SNPhood.o, regions = 1:5, plotRegionBoundaries = TRUE, sizePoints = 2, plotRegionLabels = TRUE, mergeReadGroupCounts = TRUE)
plotRegionCounts(SNPhood.o, regions = NULL, plotChr = "chr21", sizePoints = 2)


## ----allelicBias, echo=TRUE----------------------------------------------

# Run the analysis, perform no time-consuming background calculation for now
SNPhood.o = testForAllelicBiases(SNPhood.o, readGroups = c("paternal", "maternal"), verbose = FALSE)

# Extract the results of the analysis, again using the results function
names(results(SNPhood.o, type = "allelicBias"))
head(results(SNPhood.o, type = "allelicBias", elements = "pValue")[[1]], 4)


## ----allelicBias2, echo=TRUE---------------------------------------------
plotAllelicBiasResultsOverview(SNPhood.o, regions = NULL, plotChr = "chr21", signThreshold = 0.05, pValueSummary = "min")
plotAllelicBiasResultsOverview(SNPhood.o, regions = 3:5, plotRegionBoundaries = TRUE, plotRegionLabels = TRUE, signThreshold = 0.01, pValueSummary = "min")

## ----allelicBias3, echo=TRUE, fig.height=7-------------------------------
plots = plotAllelicBiasResults(SNPhood.o, region = 2, signThreshold = 0.05, readGroupColors = c("blue", "red", "gray"))
plots = plotAllelicBiasResults(SNPhood.o, region = 7, signThreshold = 0.05, readGroupColors = c("blue", "red", "gray"))

## ----clusterCountMatrix, echo=TRUE---------------------------------------
SNPhood.o = plotAndClusterMatrix(SNPhood.o, readGroup = "paternal", nClustersVec = 2, dataset = 1, verbose = FALSE)
SNPhood.o = plotAndClusterMatrix(SNPhood.o, readGroup = "paternal", nClustersVec = 5, dataset = 1, verbose = FALSE)
str(results(SNPhood.o, type = "clustering", elements = "paternal"), list.len = 3)

## ----clusterCountMatrix2, echo=TRUE--------------------------------------
SNPhood.o = plotAndClusterMatrix(SNPhood.o, readGroup = "paternal", nClustersVec = 5, dataset = 1, clustersToPlot = 2:5, verbose = FALSE)

## ----summarizeClusters, echo=TRUE----------------------------------------
p = plotClusterAverage(SNPhood.o, readGroup = "paternal", dataset = 1)

## ----associategenotype, echo=TRUE----------------------------------------
(mapping = genotypeMapping = data.frame(samples = annotationDatasets(SNPhood.o),
                                       genotypeFile = rep(fileGenotypes, 2),
                                       sampleName = c("NA10847", "NA12890")
                             ))
SNPhood.o = associateGenotypes(SNPhood.o, mapping)

## ----plotGenotypesPerSNP, echo=TRUE--------------------------------------
p = plotGenotypesPerSNP(SNPhood.o, regions = 1:20)

## ----strongAndWeakGenotype2, echo=TRUE-----------------------------------
SNPhood_merged.o = mergeReadGroups(SNPhood.o)
SNPhood_merged.o = plotAndCalculateWeakAndStrongGenotype(SNPhood_merged.o, normalize = TRUE, nClustersVec = 3)


## ----plotGenotypePerCluster, echo=TRUE-----------------------------------
p = plotGenotypesPerCluster(SNPhood_merged.o, printPlot = FALSE, printBinLabels = TRUE)
p[[1]]

