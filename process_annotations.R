#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#   Download and process infinium methylation array probe annotations
#-------------------------------------------------------------------------------

in.file <- "data/SupplementalTable1.txt.gz"
if(!file.exists(in.file))
    download.file("http://genomebiology.com/imedia/6808530410895336/supp1.gz", "data/SupplementalTable1.txt.gz")

met.annot <- read.csv(in.file,
    header=TRUE, sep="\t", nrows=485577,
    colClasses=c("character", "factor", "integer", "character",
                 "integer", "character", "character", "factor",
                 "character", "character", "character", "character",
                 "factor", "character", "character",
                 # probe filtering
                 "integer", "integer", "integer",
                 # Histone marks
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer",
                 # DMCS
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer", 
                 "integer", "integer", "integer", "integer"))
# Make the chromosome name a factor manually to get the order right
met.annot$CHR <- factor(met.annot$CHR, levels=c(1:22, c("X", "Y")))
met.annot$CHROMOSOME_36 <- factor(met.annot$CHROMOSOME_36,
    levels=c(levels(met.annot$CHR), "MULTI"))
met.annot[16:38] <- lapply(met.annot[16:38], as.logical)

# The existing dataset on GEO does not include BWA result for the X chromosome
# Update it from the csv-file in the git repo
bwa.X <- read.csv("data/bwa_hits_x.csv", colClasses=c("character", "integer"))
met.annot$bwa.multi.hit[match(bwa.X$TargetID, met.annot$TargetID)] <-
    as.logical(bwa.X$bwa.multi.hit)
met.annot$analysed <-
    with(met.annot, !snp.hit & !bwa.multi.hit & CHR %in% c(1:22, "X"))

save(met.annot, file="data/annotations.Rdata")

