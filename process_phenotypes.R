#!/opt/apps/R/3.0.1/bin/Rscript

#===============================================================================
#   Download and process phenotypes
#
#   This dataset was originally uploaded for the Nordlund & Bäcklin et al.
#   paper. Since then the subtypes of the samples have been revised, and the
#   dataset is updated accordingly in the final section of this file.
#
#   References:
#    - Genome-wide signatures of differential DNA methylation in pediatric acute
#      lymphoblastic leukemia, Nordlund & Bäcklin et al., Genome Biology 2013,
#      14:r105.
#-------------------------------------------------------------------------------

require(GEOquery)

in.file <- "data/GSE49031_series_matrix.txt.gz"
if(!file.exists(in.file))
    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49031/matrix/GSE49031_series_matrix.txt.gz", in.file)
gse <- getGEO(filename=in.file, destdir="data")


sample.type.patterns <- c(
    control   = "^MethDNA(Pos|Neg)",
    reference = "^(\\d{3}_[BT]|cd34_e1)",
    diagnosis = "^ALL_\\d+$",
    replicate = "^ALL_\\d+rep\\d",
    remission = "^Constitutional_\\d+",
    `relapse 1` = "^ALL_\\d+r1", 
    `relapse 2` = "^ALL_\\d+r2")

char.fields <- c(disease.state = "characteristics_ch1",
    immunophenotype = "characteristics_ch1.1",
    subtype = "characteristics_ch1.2",
    cell.type = "characteristics_ch1.3")

met.pheno <- with(phenoData(gse)@data, data.frame(
    id = sub(" .*$", "", description),
    geo.accession = as.character(geo_accession),
    title = as.character(title),
    sample.type = factor(NA, levels=names(sample.type.patterns)),
    lapply(char.fields, function(x) factor(sub("^.*: ", "", get(x)))),
    tissue = factor(source_name_ch1),
    stringsAsFactors=FALSE))

for(i in seq(sample.type.patterns))
    met.pheno$sample.type[grep(sample.type.patterns[i], met.pheno$id)] <-
        names(sample.type.patterns)[i]

idx <- sapply(met.pheno, function(x) "NA" %in% levels(x))
met.pheno[idx] <- lapply(met.pheno[idx], function(x){
    factor(as.character(x), levels=setdiff(levels(x), "NA"))
})

# Update subtypes
new.sub <- read.table("data/subtypes.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
new.sub$subtype <- factor(new.sub$subtype, levels=c("reference", "T-ALL", "HeH", "t(12;21)",
        "11q23/MLL", "t(1;19)", "dic(9;20)", "t(9;22)", "iAMP21", "<45chr",
        ">67chr", "biclone", "other", "normal", "no result"))
met.pheno$subtype <- new.sub$subtype[match(met.pheno$id, new.sub$id)]
met.pheno$sex <- factor(new.sub$sex)[match(met.pheno$id, new.sub$id)]

save(met.pheno, file="data/phenotypes.Rdata")

