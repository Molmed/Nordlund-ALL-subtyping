library("roxygen2")

roxygenize("nordlund2013", "nordlund2013.roxygen", unlink.target = TRUE)
system("rm -rf nordlund2013.roxygen/inst")
system("R CMD check nordlund2013.roxygen")
system("R CMD INSTALL nordlund2013.roxygen")

