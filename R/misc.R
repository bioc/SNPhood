#' @import checkmate
#' @importFrom data.table data.table :=
#' @importFrom stats binom.test
.calcBinomTestVector <- function(successes, failures, returnType="p.value", indexResult=1, conf.level=0.95) {
    
    # Check types and validity of arguments       
    assertIntegerish(successes, lower = 0)  
    assertIntegerish(failures, lower = 0) 
    assertChoice(returnType, c("p.value","conf.int", "statistic", "parameter", "estimate", "null.value", "alternative", "method", "data.name"))
    assertInt(indexResult, lower = 1) 
    assertPercentage(conf.level)
    
    dt  =  data.table(successes, failures)
    dt[, pp:=binom.test(c(successes, failures), n = NULL, conf.level = conf.level)[returnType][[1]][indexResult], by = list(successes, failures)]$pp
}

.appendList <- function(x, val) 
{
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    
    hasNames = TRUE
    if (is.null( xnames)) {
        hasNames = FALSE
    }
    
    for (v in seq_len(length(val))) {
        
        if (hasNames) {
            x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
                .appendList(x[[v]], val[[v]])
            else c(x[[v]], val[[v]])
        } else {
            x[[v]] <- if (is.list(x[[v]]) && is.list(val[[v]])) 
                .appendList(x[[v]], val[[v]])
            else c(x[[v]], val[[v]])
        }
       
    }
    x
}


#' @import checkmate
.checkAndCreateIndexFile <- function(BAMfile) {
    
    # Check types and validity of arguments     
    assertFile(BAMfile, access = "r")  
    
    indexFile = paste0(BAMfile,".bai")
    
    if (!testFile(indexFile, access = "r")) {  
        stop("Could not find index file ", indexFile,"for BAM file", BAMfile,". Try to create it with samtools index, for example", sep = "")
    }
}

#' @import checkmate Rsamtools
.constructScanBamFlagsGen <- function(isPaired = NA, 
                                       isProperPair = NA,
                                       isUnmappedQuery = NA,
                                       hasUnmappedMate = NA, 
                                       isMinusStrand = NA,
                                       isMateMinusStrand = NA,
                                       isFirstMateRead = NA, 
                                       isSecondMateRead = NA,         
                                       isSecondaryAlignment = NA,   
                                       isNotPassingQualityControls = NA,
                                       isDuplicate = NA ) {
    
    
    assertFlag(isPaired, na.ok = TRUE)
    assertFlag(isProperPair, na.ok = TRUE)
    assertFlag(isUnmappedQuery, na.ok = TRUE)
    assertFlag(hasUnmappedMate, na.ok = TRUE)
    assertFlag(isMinusStrand, na.ok = TRUE)
    assertFlag(isMateMinusStrand, na.ok = TRUE)
    assertFlag(isFirstMateRead, na.ok = TRUE)
    assertFlag(isSecondMateRead, na.ok = TRUE)
    assertFlag(isSecondaryAlignment, na.ok = TRUE)
    assertFlag(isNotPassingQualityControls, na.ok = TRUE)
    assertFlag(isDuplicate, na.ok = TRUE)    
    
    flags = scanBamFlag(isPaired                   = isPaired, 
                        isProperPair                = isProperPair, 
                        isUnmappedQuery             = isUnmappedQuery, 
                        hasUnmappedMate             = hasUnmappedMate, 
                        isMinusStrand               = isMinusStrand, 
                        isMateMinusStrand           = isMateMinusStrand,
                        isFirstMateRead             = isFirstMateRead, 
                        isSecondMateRead            = isSecondMateRead, 
                        isNotPassingQualityControls = isNotPassingQualityControls, 
                        isDuplicate                 = isDuplicate
    )
    
    flags
    
}

# .readBigFilesChunkwise <- function(fileName, buf.size=5e7, idsToKeep, colWithIDs ) {
#     s = file.info(fileName)$size
#     in.file = file(fileName, "r" )
#     buf=""
#     res = list()
#     i=1
#     while(s > 0 ) {
#         n = min(c(buf.size, s ) )
#         buf = paste(buf, readChar(in.file, n ), sep = "" )
#         
#         s = s - n
#         r = strsplit(buf, "\n", fixed=T, useBytes=T)[[1]]
#         n=nchar(buf)
#         if(substr(buf, n, n)=="\n" ) {
#             res[[i]] = r
#             buf = ""
#         } else {
#             res[[i]] = head(r,-1)
#             buf = tail(r,1)
#         }
#         i=i+1
#     }
#     close(in.file )
#     c(unlist(res), buf )
# }

#' @import checkmate
.collectFilesGen <- function(patternFiles, recursive = FALSE, ignoreCase = TRUE, verbose) {
    
    # Check types and validity of arguments   
    assertCharacter(patternFiles, any.missing = FALSE, min.chars = 1, min.len = 1)
    assertFlag(recursive)
    assertFlag(ignoreCase)
    assertFlag(verbose)  
    
    patternGiven = ifelse(length(patternFiles) == 1, TRUE, FALSE)

    if (patternGiven) {
        directory = dirname(patternFiles)
        patternFilesRed = basename(patternFiles)
        
        ####################################################
        # Identify the BAM files that need to be processed #
        ####################################################
        
        if (verbose) message("Search for files with the pattern '", patternFilesRed,"' in directory ", directory," (recursive: ", recursive,", case-sensitive:",!ignoreCase,")")
        filesToProcess.vec = list.files(path = directory, pattern = patternFilesRed, full.names = TRUE, recursive = recursive, ignore.case = ignoreCase)
        
    } else {
        
        filesToProcess.vec = unique(patternFiles)

        # Check if readable and existent
        for (i in seq_len(length(filesToProcess.vec))) {
            assertFile(filesToProcess.vec[i], access = "r")
        }

        # Check if they are all BAM files
        if (!all(grepl(filesToProcess.vec, pattern = "*.bam$", ignore.case = TRUE))) {
            stop("All specified files (",paste0(filesToProcess.vec, collapse = ","), ") must be BAM files, but at least one file did not end with \".bam\"")
        }
       

    }
   
    # Make sure only files ending with the correct ending are processed
    index = which(grepl(paste0(".","bai","$"), filesToProcess.vec, perl = TRUE))
    if (length(index) > 0) {
        filesToProcess.vec = filesToProcess.vec[-index]
    }
    
    
    if (verbose) message("Found the following files:\n", paste(filesToProcess.vec, collapse = "\n "))
    
    
    if (patternGiven) {
        nFilesToProcess = length(filesToProcess.vec)
        
        if (nFilesToProcess == 0) {
            warning(paste("No files to process in folder ", directory," that fulfill the desired criteria.", sep = ""))
        } else {
            
            if (!is.null(directory)) {
                
            } else {
                if (verbose) message("The following ", nFilesToProcess," file will be processed:\n ", paste(filesToProcess.vec, collapse = "\n "))
                
            }
        }
    }
    filesToProcess.vec
    
}


#' @import checkmate Rsamtools
.extractFromBAMGen <- function(regions.gr, file, fieldsToRetrieve = scanBamWhat(), flags, readGroupSpecificity = TRUE, simpleCigar = FALSE, reverseComplement = FALSE, verbose = TRUE) {
    
    assertClass(regions.gr, "GRanges") 
    assertFile(file, access = "r")
    assertSubset(fieldsToRetrieve, scanBamWhat(), empty.ok = FALSE)  
    assertInteger(flags, min.len = 1, any.missing = FALSE)
    assertFlag(readGroupSpecificity)
    assertFlag(simpleCigar)
    assertFlag(reverseComplement)
    assertFlag(verbose)   
    
    if (readGroupSpecificity) {
        
        param  =  ScanBamParam(which = regions.gr, what = fieldsToRetrieve, flag = flags, reverseComplement = reverseComplement, simpleCigar = simpleCigar, tag = "RG" )
        
    } else {
        
        param  =  ScanBamParam(which = regions.gr, what = fieldsToRetrieve, flag = flags, reverseComplement = reverseComplement, simpleCigar = simpleCigar)
        
    }
    
    
    return(scanBam(file, param = param))

}


#' @import checkmate
.filterReads <- function(bamObj, strand.vec, alleleCur) {
        
    assertList(bamObj, any.missing = FALSE, min.len = 1)
    assertSubset(strand.vec, c("+","-","*"))
    assertCharacter(alleleCur, min.chars = 1, len = 1)  
    
    nReadsDeleted = 0
    
    
    for (i in seq_len(length(bamObj))) {
        
        indexToDelete = c()
        
        # Filter strand
        if (!is.na(strand.vec[i]) & strand.vec[i] != "*") {
            # Filter only those with the desired and specific strand. This cannot be done before because in ScanBamParam, data types are coerced to IRangesList, which does not include strand information. 
            # The flags are also not suitable because the strand information flag is set globally and cannot be individual for each entry.
            indexToDelete = which(strand.vec[i] != bamObj[[i]]$strand)
        }
        
        # Filter read group      
        if (!is.na(alleleCur) & alleleCur != "allReadGroups") {
            
            # Check if read groups are actually defined in the file. If not, there is nothing to do. Important: the where parameter here
            if (exists("tag", where = bamObj[[i]])) { 
                if (exists("RG", where = bamObj[[i]]$tag)) {  
                    
                    
                    index = which(bamObj[[i]]$tag$RG != alleleCur)
                    if (length(index) > 0) indexToDelete = c(indexToDelete, index)
                } else {
                    warning("Read groups cannot be extracted from BAM file (no tag RG found) for region ",i,".")
                    
                }
                
            } else {
                
                warning("Read groups cannot be extracted from BAM file (no tag RG found) for region ",i,".")
            }
            
        } 
        
        if (length(indexToDelete) > 0) {
            # Delete the reads completely in strand, pos, and qwidth 
            # Potential improvement: Unlist bam[[1]] and create indexToDeleteSignal to delete all information in one step
            bamObj[[i]]$strand = bamObj[[i]]$strand[-indexToDelete]
            bamObj[[i]]$pos    = bamObj[[i]]$pos   [-indexToDelete]
            
            if (exists("qwidth", where = bamObj[[i]])) bamObj[[i]]$qwidth = bamObj[[i]]$qwidth[-indexToDelete]
            if (exists("seq"   , where = bamObj[[i]])) bamObj[[i]]$seq    = bamObj[[i]]$seq   [-indexToDelete]
            
            if (exists("tag", where = bamObj[[i]])) { 
                if (exists("RG", where = bamObj[[i]]$tag)) {  
                    bamObj[[i]]$tag$RG = bamObj[[i]]$tag$RG[-indexToDelete]      
                }
            }
            #bamObj[[i]]$rname  = bamObj[[i]]$rname [- indexToDelete]
            
            nReadsDeleted =  nReadsDeleted + length(indexToDelete)
        }
        
    }
    
    list(filteredReads = nReadsDeleted,
          bamObj        = bamObj)
    
}

#' @import checkmate
# TODO. Delete? not needed
generateDefaultReadFlags <- function(pairedEndReads = TRUE) {
    
    assertFlag(pairedEndReads)
    
    par.l = list(  
        
        "readFlag_isPaired" = TRUE, 
        "readFlag_isProperPair" = TRUE ,
        "readFlag_isUnmappedQuery" = FALSE, 
        "readFlag_hasUnmappedMate" = FALSE, 
        "readFlag_isMinusStrand" = NA, 
        "readFlag_isMateMinusStrand" = NA,
        "readFlag_isFirstMateRead" = NA, 
        "readFlag_isSecondMateRead" = NA, 
        "readFlag_isNotPrimaryRead" = FALSE,
        "readFlag_isNotPassingQualityControls" = FALSE, 
        "readFlag_isDuplicate" =  FALSE
    )
    
    if (!pairedEndReads) {
        
        par.l$readFlag_isPaired = NA 
        par.l$readFlag_isProperPair = NA 
        par.l$readFlag_isUnmappedQuery = NA 
        par.l$readFlag_hasUnmappedMate = NA 
        par.l$readFlag_isMateMinusStrand = NA 
        par.l$readFlag_isFirstMateRead = NA 
        par.l$readFlag_isSecondMateRead = NA 
        par.l$readFlag_isNotPrimaryRead = NA 
        
        
    }
    
    .constructScanBamFlags(par.l)
}

#' @import checkmate GenomeInfoDb
.getGenomeData <- function(genome, includeChrM = FALSE) {
    
    
    # See ?fetchExtendedChromInfoFromUCSC for a list of supported chromosome sizes
    assertChoice(genome, c("hg38", "hg19", "hg18", "mm10", "mm9", "dm3", "sacCer3", "sacCer2"))
    assertFlag(includeChrM)

    # library("BSgenome") may also be used together with seqlengths(NFKB)=seqlengths(Hsapiens), 
    # but this requires the download of hundreds of data for a particular genome
    
    # use a try-catch construct in case the default approach to determine chromosome sizes fails
    result = tryCatch( {
        chromInfo.df  =  suppressWarnings(fetchExtendedChromInfoFromUCSC(genome))
        chromInfo.df = chromInfo.df[which(chromInfo.df$SequenceRole == "assembled-molecule"),]
        
        rowChrM = which(chromInfo.df$UCSC_seqlevel == "chrM")
        if (!includeChrM & length(rowChrM) == 1) {  
            chromInfo.df =  chromInfo.df[-rowChrM,]
        }
        # Quickly transform to a regular data frame
        chrSizes.df = data.frame(chr = chromInfo.df$UCSC_seqlevel, size = as.numeric(chromInfo.df$UCSC_seqlength))
        
        rownames(chrSizes.df) = chrSizes.df$chr
        
        return(chrSizes.df)
        
    }, error = function(e) {
        
        warning("Could not obtain chromosome sizes using GenomeInfoDb, try fallback implementation...\n")
        assertChoice(genome, c("hg19","hg38"))

        if (genome == "hg19") {
            chrSizes = list(
                "chr1" = 249250621,
                "chr2" = 243199373,
                "chr3" = 198022430,
                "chr4" = 191154276,
                "chr5" = 180915260,
                "chr6" = 171115067,
                "chr7" = 159138663,
                "chrX" = 155270560,
                "chr8" = 146364022,
                "chr9" = 141213431,
                "chr10" = 135534747,
                "chr11" = 135006516,
                "chr12" = 133851895,
                "chr13" = 115169878,
                "chr14" = 107349540,
                "chr15" = 102531392,
                "chr16" = 90354753,
                "chr17" = 81195210,
                "chr18" = 78077248,
                "chr20" = 63025520,
                "chrY" = 59373566,
                "chr19" = 59128983,
                "chr22" = 51304566,
                "chr21" = 48129895,
                "chrM" = 16571
            )
        } else if (genome == "hg38") {
            chrSizes = list(
                "chr1" = 248956422,
                "chr2" =   242193529,
                "chr3" =     198295559,
                "chr4" =     190214555,
                "chr5" =     181538259,
                "chr6" =     170805979,
                "chr7" =     159345973,
                "chrX" =     156040895,
                "chr8" =     145138636,
                "chr9" =     138394717,
                "chr11" =     135086622,
                "chr10" =     133797422,
                "chr12" =     133275309,
                "chr13" =     114364328,
                "chr14" =     107043718,
                "chr15" =     101991189,
                "chr16" =     90338345,
                "chr17" =     83257441,
                "chr18" =     80373285,
                "chr20" =     64444167,
                "chr19" =     58617616,
                "chrY" =     57227415,
                "chr22" =     50818468,
                "chr21" =     46709983,
                "chrM" =     16569
            )
        }
        
        
        if (!includeChrM) chrSizes[["chrM"]] = NULL
        # Quickly transform to a regular data frame
        chrSizes.df = data.frame(chr = names(chrSizes), size = as.numeric(as.data.frame(chrSizes)[1,]))
    
        rownames(chrSizes.df) = chrSizes.df$chr
        chrSizes.df
    }
    
    )
    
    return(result)

    
}


#' @import checkmate
.getParameters <- function(parsValid.l, parName= NULL) {
    
    assertList(parsValid.l, min.len = 1, all.missing = FALSE)
    assert(checkNull(parName), checkCharacter(parName, any.missing = FALSE, min.chars = 1, len = 1))
    
    parName.df = data.frame(parName = names(parsValid.l), parType = as.character(unlist(parsValid.l)), stringsAsFactors = FALSE)
    
    rownames(parName.df) = names(parsValid.l)
    
    if (is.null(parName)) {
        
        return(parName.df)
        
    } else {
        
        rowIndex = which(parName.df$parName == parName)
        if (length(rowIndex) != 1) {
            stop("Could not find value of parameter ", parName," in data frame.")
        }
        
        return(parName.df[rowIndex,]$parType)
        
    }
    
}



# .initParallelisation  =   function(nCores) {
#     
#     assertInt(nCores, lower = 1)
#     
#     if (nCores == 1) return()
#     
#     nCoresAvail = BiocParallel::multicoreWorkers()
#     
#     assertInt(nCoresAvail, lower = 2)
#     
#     # Init the correct number of processors
#     if (nCoresAvail >= nCores) {
#         
#         p = BiocParallel::register(BiocParallel::MulticoreParam(nCores))
#         
#     } else {
#         
#         stop("The machine provides only ", nCoresAvail," cores, but ", nCores, " were requested. Adjust the parameter value for \"nCores\" accordingly.")
#     }
#     
# }


#' @import checkmate
#' @importFrom utils read.table
.parseBed6File <- function(path_bedFile, headerLine = FALSE, linesToParse = -1, verbose= FALSE) {
    
    assertFile(path_bedFile, access = "r")
    assertFlag(headerLine)
    assert(checkInt(linesToParse, lower = -1, upper = -1), checkInt(linesToParse, lower = 1))  
    assertFlag(verbose)
    
    userRegions.df = read.table(file = path_bedFile, header = headerLine, nrows = linesToParse, stringsAsFactors = FALSE)

    if (nrow(userRegions.df) == 0) {
        stop("No entries found in file ", path_bedFile,"\n", sep = "")
    }
    
    nCols = ncol(userRegions.df)
    nColsOrig = nCols
    
    if (nCols < 2 | nCols > 6) {
        stop("Wrong number of columns in file ", path_bedFile,": Expected 2-6, found only", nCols,". At least two columns are needed (chromosome, start\n", sep = "")      
    } 
    
    # Check column 1, always present
    assert(checkInteger(userRegions.df[,1], lower = 1, any.missing = FALSE, min.len = 1), checkCharacter(userRegions.df[,1], min.chars = 4, pattern = "chr", any.missing = FALSE, min.len = 1))
    
    # Column 2: (Start) position, always present
    assertInteger(userRegions.df[,2], lower = 0, any.missing = FALSE)
    
    # Set the column names
    if (nCols == 2)  {
        colnames(userRegions.df) = c("chr","start")
    }
    if (nCols == 3)  {
        colnames(userRegions.df) = c("chr","start","annotation")
        assertCharacter(as.character(userRegions.df[,3]), min.chars = 1, any.missing = FALSE, min.len = 1)
    }
    if (nCols == 4)  {
        colnames(userRegions.df) = c("chr","start","end","annotation")
        assertInteger(userRegions.df[,3], lower = 0, any.missing = FALSE)        
        assertCharacter(userRegions.df[,4], min.chars = 1, any.missing = FALSE, min.len = 1)
        
    }
    if (nCols == 5)  {
        colnames(userRegions.df) = c("chr","start","end","annotation","strand")
        assertInteger(userRegions.df[,3], lower = 0, any.missing = FALSE)        
        assertCharacter(userRegions.df[,4], min.chars = 1, any.missing = FALSE, min.len = 1)
        assertSubset(userRegions.df[,5], c("+","-","*","."))
    }
    if (nCols == 6)  {
        colnames(userRegions.df) = c("chr","start","end","annotation","score","strand")
        assertInteger(userRegions.df[,3], lower = 0, any.missing = FALSE)        
        assertCharacter(userRegions.df[,4], min.chars = 1, any.missing = FALSE, min.len = 1)
        assertInteger(userRegions.df[,5], any.missing = FALSE)  
        assertSubset(userRegions.df[,6], c("+","-","*","."))
        
    }
    
    if (nCols < 3) {
        
        userRegions.df$annotation = paste0("seq", seq_len(nrow(userRegions.df)),":", userRegions.df[,1],"_", userRegions.df[,2])
        
    }
    
    if (nCols < 6) {
        userRegions.df$score = rep(".", nrow(userRegions.df))
    }
    
    if (nCols < 5) {
        userRegions.df$strand = rep("*", nrow(userRegions.df))
        if (verbose) message("   Add strand information for all regions. Assume strand is irrelevant (*)...")
    } 
    
    if (nCols < 4) {
        
        userRegions.df$end = as.numeric(userRegions.df[,2])
    } 
    
    userRegions.df$annotation = as.character(userRegions.df$annotation)
    
    
    # Do a few more tests to see if the content is valid  
    userRegions.df$strand[which(userRegions.df$strand == ".")] = "*"
    
    # Check if the chromosome names conform to the "standard" chr{Nr}
    if (length(which(grepl("^chr", userRegions.df$chr, perl = TRUE))) != nrow(userRegions.df)) {
        if (verbose) message("   Modify chromosome names because they do not start with \"chr\"")
        userRegions.df$chr = paste0("chr", userRegions.df$chr)
    }
    
    # Test if end is an integer
    if (!all(floor(userRegions.df$end) == userRegions.df$end) | !all(floor(userRegions.df$start) == userRegions.df$start)) {
        stop("Either start or end positions do not contain valid numbers. Please check if you, for example, accidentally include the header line. In that case, adjust the parameter headerLine accordingly.")
    }
    
    
    if (max(userRegions.df$end) > 536870912) {
        stop("At least one end position is larger than 536870912 (", max(userRegions.df$end),"), which is not supported by the IRange package\n", sep = "")
    }
    
    userRegions.df = userRegions.df[, c("chr","start","end","annotation","score","strand")]
    
    list(userRegions.df, nColsOrig)
    
}

#' @import checkmate S4Vectors IRanges GenomeInfoDb GenomicRanges
.parseAndProcessUserRegions <- function(par.l, chrSizes.df, verbose = TRUE)  {
    
    assertList(par.l, min.len = 1, all.missing = FALSE)  
    assertDataFrame(chrSizes.df, min.rows = 1, min.cols = 1)
    assertFlag(verbose)
    
    listFunction = .parseBed6File(par.l$path_userRegions, par.l$headerLine, par.l$linesToParse, verbose = verbose)
    if (verbose) message(" Parse the user regions file")
    
    userRegions.df = listFunction[[1]]
    nColOrig       = listFunction[[2]]
    
    nRegions = nrow(userRegions.df)
    
    if (verbose) message("  Finished parsing. Number of entries processed: ", nRegions)
    
    
    
    if (par.l$linesToParse != -1) {
        warning("The parameter \"linesToParse\" has been set to a non-default value. Only a subset of all lines have been parsed\n")
    }
    
    # Normalize everything to fully closed one-based coordinates. GRange coordinates are fully closed, and 1 based! See http://master.bioconductor.org/help/course-materials/2014/SeattleFeb2014/Ranges.pdf.
    
    # Adjust positions by one if zero-based coordinates
    if (par.l$zeroBasedCoordinates) {
        if (verbose) message("  Increase start and end coordinates by one because they are zero-based")
        userRegions.df$start = userRegions.df$start + 1
        userRegions.df$end   = userRegions.df$end  + 1
    }
    
    # Similar to the BED standard with half open intervals, we also need those
    if (par.l$endOpen & nColOrig != 2) { # only do it when the user did specify end coordinates already
        if (verbose) message("  Decrease end coordinates by one because coordinates have to be closed at the end")  
        userRegions.df$end   = userRegions.df$end  - 1
    }
    
    if (par.l$startOpen) {
        if (verbose) message("  Decrease start coordinates by one because coordinates have to be closed at both start and end according to the BAM standard)")
        userRegions.df$start   = userRegions.df$start  + 1
    }
    
    
    # Check if some of the positions are larger than the chromosome sizes
    index.vec = which(userRegions.df$start > chrSizes.df[userRegions.df$chr,]$size)
    if (length(index.vec) > 0) {    
        stop(length(index.vec), " regions (starting with annotations ", paste(userRegions.df$annotation[head(index.vec)], collapse = ","),") have invalid start positions that exceed the lengths of the corresponding chromosomes. Has the genome assembly version be specified correctly?")
    }
    index.vec = which(userRegions.df$end > chrSizes.df[userRegions.df$chr,]$size)
    if (length(index.vec) > 0) {    
        stop(length(index.vec), " regions (starting with annotations ", paste(userRegions.df$annotation[head(index.vec)], collapse = ","),") have invalid end positions that exceed the lengths of the corresponding chromosomes. Has the genome assembly version be specified correctly?")
    }
    
    
    if (!all((userRegions.df$end - userRegions.df$start + 1) == 1)) {
        stop("Only SNPs are supported, not regions of length > 1. At least one region has a length >1.")
    }
    
    
    userRegions.df$SNPPos  = userRegions.df$start
    
    # Add the regionSize Parameters to the boundaries of the regions
    userRegions.df$start  = userRegions.df$start - par.l$regionSize
    userRegions.df$end    = userRegions.df$end + par.l$regionSize
    
    index1 = which(userRegions.df$start < 0)
    if (length(index1) > 0) {
        warning("For ", length(index1), " user regions, the chromosome starts have been reached after incorporating the regionSize parameter. The start positions have therefore been automatically set to 0.")  
        userRegions.df$start[which(userRegions.df$start < 0)] = 0
    }
    
    index2 = which((userRegions.df$end + par.l$regionSize) > chrSizes.df[userRegions.df$chr,]$size)
    if (length(index2) > 0) {
        warning("For ", length(index2), " user regions, the chromosome ends have been reached after incorporating the regionSize parameter. The end positions have therefore been automatically set to the size of the chromosome.")  
        userRegions.df$end[which(userRegions.df$end > chrSizes.df[userRegions.df$chr,]$size)] = chrSizes.df[userRegions.df$chr[which(userRegions.df$end > chrSizes.df[userRegions.df$chr,]$size)],]$size
        
    }
    
    userRegions.df$length = userRegions.df$end - userRegions.df$start + 1 # Because intervals are fully-closed  
    
    # Check for duplicate rows and filter them
    nDuplicateRows = length(which(duplicated(userRegions.df)))
    if (nDuplicateRows > 0) {
        userRegions.df = userRegions.df[!duplicated(userRegions.df),]
        warning(nDuplicateRows, " duplicate region removed in file ", par.l$path_userRegions," out of ", nRegions," regions. New number of regions: ", nrow(userRegions.df))
        
    }
    
    # Check if the IDs are unique. If not, throw an error
    if (length(unique(userRegions.df$annotation)) != length(userRegions.df$annotation)) {
        warning("All IDs from the user regions file must be unique. However, some of the IDs were identical and have been changed by adding \"_1\", \"_2\", etc to make them unique")
        duplicatedIDs = duplicated(userRegions.df$annotation)
        
        # Change the IDs so that they are all unique
        currentDuplicationIndex = 1
        for (i in 1:nrow(userRegions.df)) {

            if (duplicatedIDs[i]) {
                
                currentDuplicationIndex = currentDuplicationIndex + 1
                userRegions.df$annotation[i] = paste0(userRegions.df$annotation[i],"_",currentDuplicationIndex)
                
            } else {
                currentDuplicationIndex = 1
            }
         
        }
        
    }
    
    assertCharacter(userRegions.df$annotation, unique = TRUE, any.missing = FALSE)
    
    # Finally, order the user regions by chromosome and strand. This is important for the scanBam function so that it does not change the order of the regions
    #userRegions.df = userRegions.df[order(userRegions.df$chr, userRegions.df$start),] TODO: test if working and compare results
    
    gr = GRanges(seqnames = userRegions.df$chr, 
                 ranges = IRanges(start = userRegions.df$start, end = userRegions.df$end), 
                 strand = userRegions.df$strand, 
                 annotation = userRegions.df$annotation, 
                 SNPPos = userRegions.df$SNPPos
                 )
    
    gr = sortSeqlevels(gr)
    gr = sort(gr)
    
    gr
    
    
}

#' @import checkmate SummarizedExperiment
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts sizeFactors
.scaleLibraries <- function(rawCounts.m, verbose = TRUE) {
    
    # Check types and validity of arguments     
    assertMatrix(rawCounts.m, min.rows = 1, min.cols = 1)
    assertFlag(verbose)
    
    # Use DeSeq to compute the size factors and scale the libraries
    # Make sure we have column names
    colnames(rawCounts.m) = paste0("sample",1:ncol(rawCounts.m))

    colData = data.frame(condition = colnames(rawCounts.m)) # column names don't matter because of the specific design matrix (only intercept fitting)
    dds = DESeqDataSetFromMatrix(countData = rawCounts.m, colData = colData, design = ~1)
    
    # Test how many genes will be skipped for determinging the size factors
    log.vec = apply(rawCounts.m,1, function(x) {any(x == 0)})
    nGenesNorm = length(which(!log.vec))
    if (verbose) message("Number of genes on which library normalization using DeSeq2 is based on out of ", nrow(rawCounts.m),":", nGenesNorm)
    
    if (nGenesNorm == 0) {
        stop("All genes have at least one zero count in one of the samples. Size factor normalization cannot be performed.")
    }
    
    # Print a warning if the number of rows with at least one zero is too high, as the size factors are only based on rows with no zeroes altogether
    percThreshold = 0.5
    if (nGenesNorm < (percThreshold*nrow(rawCounts.m)) | nGenesNorm < 50 ) {
        warning("The proportion of genes that can be used for normalization among all samples is lower than ", percThreshold," or 50 altogether. Size factors will still be determined but may be inaccurate.")
    }
    
    
    # Regular DESeq2 workflow not needed here
    dds  =  estimateSizeFactors(dds)
    
    res.l = list()
    res.l$adjCounts = counts(dds, normalized = TRUE)
    res.l$sizeFactors =  sizeFactors(dds)
    
    #TODO: Check here if the names are correct - seems to be the wrong name. 
    names(res.l$sizeFactors) = colnames(rawCounts.m)
    
    res.l

}


# 
# checkAndInstallMissingPackages <- function(packages.vec, verbose = verbose) {
#     
#     if (verbose) cat ("Trying to automatically install missing packages. If this fails, install them manually\n")
#     
#     packagesToInstall = setdiff(packages.vec, rownames(installed.packages()))
#     if (verbose) cat ("Could not find the following packages: ", paste(packagesToInstall , collapse=","),"\n")
#     
#     if (length(packagesToInstall) > 0) {
#         
#         # Try the normal installation routine
#         install.packages(packagesToInstall)  
#         
#         # Try Bioconductor
#         source("http://bioconductor.org/biocLite.R")
#         for (packageCur in packagesToInstall) {
#             biocLite(packageCur)
#         }
#     }
#     
# }



# coverageSummary <- function(bamFile, regions, par.l, ...) {
#     
#     param  =  ScanBamParam(what=c('pos', 'qwidth'),
#                           which = regions,
#                           flag=constructScanBamFlags(par.l))
#     x  =  scanBam(bamFile, ..., param=param)[[1]]
#     
#     cvg = coverage(IRanges(x[["pos"]], width=x[["qwidth"]]))
#     
#     res = list()
#     res$mean   = mean(cvg)
#     res$median = median(cvg)
#     res$max    = min(cvg)
#     res$min    = max(cvg)
#     
#     res
# }
# 
# 
# createDefaultIDCurTime <- function() {
#     format(Sys.time(),"%Y-%m-%d-%H:%M:%S")
# }

# optimizeMemUsage <- function(filtered.vec) {
#     names(filtered.vec) = NULL
#     if (length(unique(filtered.vec)) == 1) {
#         filtered.vec = unique(filtered.vec)
#     }
#     
#     filtered.vec
# }


# .overrideParameters <- function(par.l, overridePar) {
#     
#     # Check types and validity of arguments     
#     assertList(par.l, min.len = 1, all.missing = FALSE)
#     assert(checkNull(overridePar), checkCharacter(overridePar, any.missing = FALSE, min.chars = 1, min.len = 1))
#     
#     # Implement a flexible way to allow overriding of the parameters for the configuration file. 
#     if (!is.null(overridePar)) {
#         args.vec = unlist(strsplit(overridePar, ",", fixed=T))
#         
#         args.l = strsplit(args.vec, "=", fixed=T)
#         
#         for (i in seq_len(length(args.l))) {
#             
#             # Remove white spaces
#             args.trailed.vec = sub("^\\s+","", args.l[[i]])
#             args.trailed.vec = sub("\\s+$","", args.trailed.vec)
#             
#             parName  =  args.trailed.vec[1]
#             parValue =  .changeVariableType(args.trailed.vec[2],.getParameters(.getListOfSupportedParameters(), parName))
#             
#             if (!exists(parName, where=par.l)) {
#                 
#                 if (verbose) message("  The parameter ", parName," has been added. New value: ", parValue)
#                 par.l[parName] = parValue
#                 
#             } else {
#                 
#                 if (par.l[parName] != parValue) {
#                     par.l[parName] = parValue
#                     if (verbose) message("  The parameter value for the parameter ", parName," from the configuration file has been changed. New value: ", parValue)
#                     
#                 }
#             }
#             
#         }
#         
#     }
#     
#     par.l
# }

# parseConfigFile <- function(file_configFile, verbose= FALSE, overridePar = NULL) {
#     
#     # Check types and validity of arguments         
#     assertFile(file_configFile, access = "r")
#     assertFlag(verbose, len = 1, any.missing = FALSE) 
#     assert(checkNull(overridePar), checkCharacter(overridePar, any.missing = FALSE, min.chars = 1, len = 1)) 
#     
#     config.vec = readLines(file_configFile)
#     
#     # Delete lines in the config line that are not of relevance
#     config.vec = config.vec[-which(config.vec == "")]
#     config.vec = config.vec[-which(grepl("^\\s*#", config.vec, perl = TRUE))]
#     
#     par.tmp.l = strsplit(config.vec, "=", fixed=T)
#     
#     # Test a few things to make sure everything was parsed correctly and that everything was specified as supposed to
#     parName.vec = sapply(seq_len(length(par.tmp.l)), FUN = function(x) {gsub(" ","", par.tmp.l[[x]][1])})
#     
#     pars.df = .getParameters(.getListOfSupportedParameters())
#     par.l = list()
#     for (i in 1:length(par.tmp.l)) {
#         
#         if (length(par.tmp.l[[i]]) != 2) {
#             stop("Line ", config.vec[i]," has an invalid format.\n", sep = "")
#         }
#         
#         # For the value, we only want to eliminate white spaces if it is not within a path
#         parValue = sub("^\\s+","", par.tmp.l[[i]][2])
#         parValue = sub("\\s+$","", parValue)
#         parValue = gsub("\\\"","", parValue)
#         
#         # If the parameter name already exists, overwrite
#         if (exists(parName.vec[i], where=par.l)) {
#             warning("Parameter \"", parName.vec[i] , "\" defined multiple times. Old value: ", par.l[[parName.vec[i]]],", new value: ", parValue,". Overwriting old value.\n")
#             par.l[[which(names(par.l) == parName.vec[i])]] = parValue
#             next
#         }
#         
#         if (parName.vec[i]  %in% pars.df$parName) {
#             
#             par.l[[parName.vec[i]]] = .changeVariableType(parValue, pars.df[parName.vec[i], ]$parType)
#             
#         } 
#         
#     }
#     
#     par.l = .overrideParameters  (par.l, overridePar)
#     
#     
#     # Check the validity of the number of cores parameters and initialize
#     multicore.param = .initParallelisation(par.l$nCores)
#     
#     par.l
#     
# }


# .changeVariableType <- function(var, type) {
#     
#     # Check types and validity of arguments     
#     assertChoice(type, c("character", "logical", "numerical", "integer"))
#     
#     if (type == "character") {
#         
#         return(var)
#         
#     } else if (type == "logical") {
#         
#         return(as.logical(var))
#         
#     } else if (type == "numerical") {
#         
#         return(as.numeric(var))
#         
#     } else if (type == "integer") {
#         
#         return(as.integer(var))
#         
#     } 
#     
# }



## Open a connection separately to hold the cursor position
# file.in <- file(file_vcf, 'rt')
## Use a loop to read in the rest of the lines
# line <- readLines(file.in, n=bufSize)
# 
# while (length(line)) {
#     
#     # Skip comment lines
#     line = line[-grepl("^#", line)]
#     
#     line.split <- strsplit(line, '\t')
#     
#     lengthLine = sapply(line.split,length)
#     assertIntegerish(lengthLine, lower = nColumnsMin, any.missing = F)
#     assertNumber(unique(lengthLine))
#     
#     #todo: parse all columns from the 10th on, use individual ID as elem for data frame, also parse REF and ALT allele
#     
#     ids = sapply(line.split,"[[",3)
#     rowsToKeep = which(ids %in% IDsToKeep)
#     
#     if (length(rowsToKeep) > 0) {           
#         mapping = which(IDsToKeep %in% ids[rowsToKeep])
#         genotype.vec[mapping] = sapply(line.split[rowsToKeep],"[[",10)
#     }
#     
#     line<- readLines(file.in, n=bufSize)
# }
# close(file.in)
# 
# return(genotype.vec)
# 
# }

detachAllPackages <- function() {
    
    basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    
    package.list <- search()[ifelse(unlist(gregexpr("package:",search())) == 1,TRUE,FALSE)]
    
    package.list <- setdiff(package.list,basic.packages)
    
    if (length(package.list) > 0)  for (package in package.list) detach(package, character.only = TRUE)
    
}


#' @import checkmate
#' @importFrom stats sd
.normalizeMatrixForClustering <- function(target.m, verbose = TRUE) {
    
    assertMatrix(target.m, any.missing = FALSE, min.rows = 1, min.cols = 1)
    assertFlag(verbose)
    
    valueToAdd = 0.01
    
    SdS = apply(target.m, 1, sd)
    zeroSds = length(which(SdS == 0))
    
    if (zeroSds > 0) {
        if (verbose) message(zeroSds," regions had a standard deviation of 0 across all bins. For normalization and plotting purposes, a small value of ", 
                             valueToAdd," was added")
        
        SdS = replace(SdS, SdS == 0,  valueToAdd)
    }
    
    target.m = (target.m - rowMeans(target.m)) / SdS
    
    target.m
}


#' @import checkmate 
#' @importFrom cluster clara silhouette sortSilhouette
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
.pamClustering  = function(SNPhood.o, target.m, nClustersVec,...){
    
    .checkObjectValidity(SNPhood.o)
    assertMatrix(target.m, any.missing = FALSE, min.rows = 1, min.cols = 1)
    assertIntegerish(nClustersVec, lower = 2, any.missing = FALSE, min.len = 1)   
    
    
    target.m_new = target.m[complete.cases(target.m),]
    if (nrow(target.m) > nrow(target.m_new)) {
        warning((nrow(target.m) - nrow(target.m_new)), " rows in the matrix contained NAs and they have been removed.")
    }
    target.m = target.m_new
    
    # Init the objects that are returned
    clusterplot.list = vector("list", length(nClustersVec))
    avgSilInfo = c()
    q = vector("list", length(nClustersVec))
    
    for (i in 1:length(nClustersVec)) {
        
        .pamClustering = clara(target.m, nClustersVec[i])
        
        pam.clusters.m = as.matrix(.pamClustering$clustering)
        
        pam.clusters.m = as.data.frame(pam.clusters.m, stringsAsFactors = FALSE)
        
        clusterplot.df = merge(x = pam.clusters.m, y = target.m, by = "row.names")
        clusterplot.df = clusterplot.df[, -1]
        
        colnames(clusterplot.df) = c("cluster", annotationBins(SNPhood.o))
        clusterplot.df = clusterplot.df[order(clusterplot.df$cluster, decreasing = TRUE),]
        
        # Collect results   
        avgSilInfo = c(avgSilInfo, .pamClustering$silinfo$avg.width)
        clusterplot.list[[i]] =  clusterplot.df
        q[[i]] = .generateClusterPlot(clusterplot.df, SNPhood.o, clustersToPlot = NULL)
        
        
    }
    
    res.l = list(clusterplot.list, avgSilInfo, q)
    names(res.l) = c("clusteringMatrix","averageSilhouette","plots")
    
    return(res.l)
    
}


#' @importFrom lattice levelplot panel.levelplot panel.abline
#' @importFrom grDevices gray
.generateClusterPlot <- function(clusterplot.df, SNPhood.o, clustersToPlot = NULL) {
    
    .checkObjectValidity(SNPhood.o)
    assertDataFrame(clusterplot.df, min.rows = 1, ncols = nBins(SNPhood.o) + 1)
    assert(checkNull(clustersToPlot), checkSubset(clustersToPlot, unique(clusterplot.df$cluster)))
    
    if (testNull(clustersToPlot)) {
        clustersToPlot = unique(clusterplot.df$cluster)
    }
    
    # Select only the clusters the user wants
    
    clusterplot.df = clusterplot.df[which(clusterplot.df$cluster %in% clustersToPlot),]
    clusterplot.df = clusterplot.df[order(clusterplot.df$cluster, decreasing = TRUE),]
    clusterplot.df = t(clusterplot.df[, -1])
    
    rownames(clusterplot.df) = SNPhood.o@internal$plot_labelBins
    colnames(clusterplot.df) = 1:ncol(clusterplot.df)
    
    nNumbersYAxis = 5
    yPos = seq(1, ncol(clusterplot.df), max(floor(ncol(clusterplot.df) / nNumbersYAxis), 1))
    
    colors = gray(100:0/100) # white for low read counts
    #colors = gray(0:100/100) # black for low read counts
    
    p = levelplot(as.matrix(clusterplot.df), xlab = .getBinLabelXAxis(SNPhood.o) , ylab = "User regions",
                  scales = list(y = list(at = yPos)),
                  col.regions = colors, aspect = "fill", 
                  panel = function(...) {
                      panel.levelplot(...)
                      panel.abline(v = SNPhood.o@internal$plot_origBinSNPPosition , type = "o", pch = 22, lty = 2, col = "red")
                  })
    
    p
}


#' @importFrom BiocParallel multicoreWorkers MulticoreParam
initBiocParallel <- function(nWorkers) {
    
    if (nWorkers > multicoreWorkers()) {
        warning("Requested ", nWorkers, " CPUs, but only ", multicoreWorkers(), "are available and can be used.")
        nWorkers = multicoreWorkers()
    }
    
    MulticoreParam(workers = nWorkers, progressBar = TRUE, stop.on.error = TRUE)
}

execInParallelGen <- function(nCores, iteration, functionName, ...) {
    
    res.l = list()
    # Call the function in parallel

    if (nCores > 1) {
        
        res.l = tryCatch( {
            bplapply(iteration, functionName, ..., BPPARAM = initBiocParallel(nCores))

        }, error = function(e) {
            warning("An error occured while executing the function with multiple CPUs. Trying again using only only one CPU...")
            lapply(iteration, functionName, ...)
        }
        )
        
    } else {
        res.l = lapply(iteration, functionName, ...)
    }
    
    res.l
    
}
