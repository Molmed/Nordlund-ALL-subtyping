library("roxygen2")

roxygenize("nordlund2013")
system("R CMD check nordlund2013")
system("R CMD INSTALL nordlund2013")

