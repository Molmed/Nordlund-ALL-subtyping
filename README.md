ALL subtyping based on 450k methylation assay
======================

Original pubilcation: [Nordlund](http://scholar.google.se/citations?user=ZztFeTEAAAAJ) et al. 2014, 
[DNA methylation-based subtype prediction for pediatric acute lymphoblastic leukemia](http://www.clinicalepigeneticsjournal.com/content/7/1/11/abstract)

This code was produced by the groups of [Cancer Pharmacology and Computational Medicine](http://www.medsci.uu.se/research/Cancer/Cancer+Pharmacology+and+Computational+Medicine/) and [Molecular Medicine](http://www.molmed.medsci.uu.se/) at the [Department of Medical Sciences](http://www.medsci.uu.se) at [Uppsala University](http://www.uu.se).

System requirements
-------------------
### Software
The code is written for a Unix or Linux operating system with [R](http://r-project.org) version 3.0.1 or later, but can be modified to run under Windows fairly easily. Packages 
`doMC` (1.3-0),
[`GEOquery`](http://www.bioconductor.org/packages/2.12/bioc/html/GEOquery.html) (2.26.2),
[`pamr`](http://www-stat.stanford.edu/~tibs/PAM/Rdist/doc/readme.html),
`predict` (2.1-8) and 
[`roxygen2`](http://roxygen.org/) (2.2.2), 
are required, but are installed automatically. The `predict` package was developed inhouse and will soon be released on CRAN (manuscript in preparation).

The `analyse450k` package appearing in the end of `analyze_final.R` is an inhouse package for data management, that will not be distributed. Instead, when the manuscript is accepted for publication and the validation dataset is made publicly available on [GEO](http://www.ncbi.nlm.nih.gov/geo/), it will be incorporated into `setup.R`.

### Hardware
At peak memory, 20 GB of RAM is required (in `process_methylation.R`). If you wish to run the analysis in multicore mode 13 GB/core is required.

Instructions
------------
The most convenient way to run the analysis is to clone the repo to your computer and run the files as shown below (commands for unix/linux).

    git clone git@github.com:Molmed/Nordlund-ALL-subtyping.git

### Replication of the study
The replicate the entire training procedure to create the classifier used for
final prediction run the following commands.

    cd Nordlund-ALL-subtyping
    R -f setup.R
    R -f analyze_tune.R
    R -f analyze_final.R

Notice that `analyze_tune.R` and `analyze_final.R` are designed to run on multiple CPU cores, and that you manually need to specify how many to use by editing the files, setting the variable `number.of.cores`. `analyze_tune.R` can also be run on multiple machines to reduce computation time further, see the comments in beginnig of the file for instructions.

- `setup.R` will download all data from [GEO](http://www.ncbi.nlm.nih.gov/geo/) and prepare it for use in R.
- `analyze_tune.R` will perform the doubly cross validated feature selection routine necessary for model parameter tuning and performance estimation.
- `analyze_final.R` will perform the model tuning, estimate performance and build the the final model that was presented in the paper. It will also produce some tables of the results.

### Use the final classifier directly
To use the trained classifier on your own data follow the instructions in the
script [classifier/classifier.R](classifier/classifier.R).

