Manuscript: ALL subtyping based on 450k methylation assay
======================

[Nordlund J](http://scholar.google.se/citations?user=ZztFeTEAAAAJ), [Bäcklin C](http://scholar.google.se/citations?user=ZMtuZXsAAAAJ), ...

Original publication: [to be added](#)

This code was produced by the groups of [Cancer Pharmacology and Computational Medicine](http://www.medsci.uu.se/research/Cancer/Cancer+Pharmacology+and+Computational+Medicine/) and [Molecular Medicine](http://www.molmed.medsci.uu.se/) at the [Department of Medical Sciences](http://www.medsci.uu.se) at [Uppsala University](http://www.uu.se).

System requirements
-------------------
R version 3.0.1 or later. Packages 
`doSNOW`,
[`GEOquery`](http://www.bioconductor.org/packages/2.12/bioc/html/GEOquery.html)
[`pamr`](http://www-stat.stanford.edu/~tibs/PAM/Rdist/doc/readme.html),
`predict` and 
[`roxygen2`](http://roxygen.org/), 
are required, but are installed automatically. The `predict` package was developed inhouse and will soon be released on CRAN (manuscript in preparation).

Figures were annotated using the [`biomaRt`](http://www.bioconductor.org/packages/2.12/bioc/html/biomaRt.html) and [`GenomeGraphs`](http://www.bioconductor.org/packages/2.12/bioc/html/GenomeGraphs.html) packages. However, code for producing the figures is not included.

Instructions
------------
The most convenient way to run the analysis is to clone the repo to your computer and run the files in the following order (commands for unix/linux):

    git clone git@github.com:Molmed/Nordlund-2013.git
    cd Nordlund-2013
    R -f setup.R
    R -f analyze_tune.R
    R -f analyze_final.R

- `setup.R` will download all data from [GEO](http://www.ncbi.nlm.nih.gov/geo/) and prepare it for use in R.
- `analyze_tune.R` will tune model parameters and estimate performance.
- `analyze_final.R` will build the the final model that was presented in the paper, and produce some tables of the results.

Notice that `analyze_tune.R` is designed for parallelization on a computer cluster or powerful multicore machine, whereas `analyze_final.R` is designed to be run on a single core on a computer with at least 15 GB RAM. Please review the parallelization settings in `analyze_tune.R` for best performance (it is single core by default).