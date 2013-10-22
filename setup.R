#===============================================================================
#
#   This script sets up the environment required for the analysis code.
#   After this script has finised you can run `analyze_tune.R` and
#   `analyze_final.R`.
#
#-------------------------------------------------------------------------------


#-------------------------------o
#   Make directory structure

dir.create("data", showWarnings=FALSE)
dir.create("tuning", showWarnings=FALSE)
dir.create("results", showWarnings=FALSE)


#-------------------------------o
#   Install required packages

required.pkg <- c("pamr", "predict", "roxygen2", "doSNOW")
required.bioc.pkg <- c("GEOquery")
installed.pkg <- rownames(installed.packages())

required.pkg <- required.pkg[!required.pkg %in% installed.pkg]
required.bioc.pkg <- required.bioc.pkg[!required.bioc.pkg %in% installed.pkg]

if(length(required.pkg) > 0) install.packages(required.pkg)
if(length(required.bioc.pkg) > 0){
    source("http://bioconductor.org/biocLite.R")
    biocLite(required.bioc.pkg)
}

# Compile C functions
setwd("package")
source("pack.R")
setwd("..")


#-------------------------------o
#   Download and prepare data

if(!file.exists("data/phenotypes.Rdata")){
    source("process_phenotypes.R")
}

if(!file.exists("data/annotations.Rdata")){
    source("process_annotations.R")
}

if(!file.exists("data/methylation.Rdata")){
    source("process_methylation.R")
}

