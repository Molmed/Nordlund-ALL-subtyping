Manuscript: ALL subtyping based on 450k methylation assay
======================

[Nordlund J](http://scholar.google.se/citations?user=ZztFeTEAAAAJ), [Bäcklin C](http://scholar.google.se/citations?user=ZMtuZXsAAAAJ), ...

Original publication: [to be added](#)

This code was produced by the groups of [Cancer Pharmacology and Computational Medicine](http://www.medsci.uu.se/research/Cancer/Cancer+Pharmacology+and+Computational+Medicine/) and [Molecular Medicine](http://www.molmed.medsci.uu.se/) at the [Department of Medical Sciences](http://www.medsci.uu.se) at [Uppsala University](http://www.uu.se).

System requirements
-------------------
R version 3.0.1 or later. Packages `pamr`, `predict`, `roxygen2` and [`GEOquery`](http://www.bioconductor.org/packages/2.12/bioc/html/GEOquery.html) are required, but are installed automatically. 

Figures were annotated using the [`biomaRt`](http://www.bioconductor.org/packages/2.12/bioc/html/biomaRt.html) and [`GenomeGraphs`](http://www.bioconductor.org/packages/2.12/bioc/html/GenomeGraphs.html) packages. However, code for producing the figures is not included.

Instructions
------------
Download all scripts in this repo to a new directory. The most convenient way to do this on a linux/unix system is to clone the whole repo. Then run the files `setup.R` and `analyse.R`.

    git clone git@github.com:Molmed/Nordlund-2013.git
    cd Nordlund-Backlin-2013
    R -f setup.R
    R -f analyze.R

`setup.R` will download all data from [GEO](http://www.ncbi.nlm.nih.gov/geo/), prepare it for use in R and store it in a new subfolder called `data`. `analyse.R` will run the analyses, produce the results and save them in a new subfolder called `results`. Plots are not produced.
