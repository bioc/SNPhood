## ---- results='asis',echo=FALSE------------------------------------------

files.df   = data.frame(signal = c("file1.bam", "file2.bam","file3.bam"),
                        input = NA,
                        stringsAsFactors = FALSE)
knitr::kable(files.df)

cat("\n\n")


files.df   = data.frame(signal = c("file1.bam", "file2.bam","file3.bam"),
                        input = c("input1.bam","input1.bam","input2.bam, input3.bam"),
                        individual = c("S1","S1","S2"),
                        genotype = c("file1.vcf:colName1","file1.vcf:colName1","file1.vcf:colName2"),
                        stringsAsFactors = FALSE)
knitr::kable(files.df)

