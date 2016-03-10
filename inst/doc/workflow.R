## ----<knitr, echo=FALSE, message=FALSE, results="hide"-------------------
library("knitr")
opts_chunk$set(
  tidy = TRUE,
  dev = "png",
  fig.width = 11,
  cache = TRUE,
  message = FALSE)

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

