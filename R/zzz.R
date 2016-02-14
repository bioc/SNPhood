
.onAttach <- function(libname, pkgname) {
    
    lengthDashes = 146
    message = paste0("\n",
                     paste0(rep("-",lengthDashes),collapse = ""),
                     "\n",
                     "| Welcome to the SNPhood package. See the Vignette (type browseVignettes(\"SNPhood\") and the help pages for how to use SNPhood for your analyses. |\n",
                     paste0(rep("-",lengthDashes),collapse = ""),
                     "\n"
                     )
    packageStartupMessage(message)
    
    # Turn off scientific notation
    options(scipen = 999) 
    
    # Print warnings as they occur. Don't use but instead print all unique warnings at the end of the run
    #options(warn=1)
}
