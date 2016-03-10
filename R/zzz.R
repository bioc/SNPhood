
.onAttach <- function(libname, pkgname) {
    
    lengthDashes = 114
    message = paste0("\n",
                     paste0(rep("-",lengthDashes),collapse = ""),
                     "\n",
                     "|       Welcome to the SNPhood package and thank you for using our software. This is SNPhood version ", packageVersion("SNPhood"),".      |\n",
                     "| See the vignettes (type browseVignettes(\"SNPhood\") or the help pages for how to use SNPhood for your analyses. |\n",
                     "|       Thank you for using our software. Please do not hesitate to contact us if there are any questions.       |\n",
                     paste0(rep("-",lengthDashes),collapse = ""),
                     "\n"
                     )
    packageStartupMessage(message)
    
    # Turn off scientific notation
    options(scipen = 999) 
    
    # Print warnings as they occur. Don't use but instead print all unique warnings at the end of the run
    #options(warn=1)
}
