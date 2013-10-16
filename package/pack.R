library("roxygen2")

roxygenize("nordlund2014", "nordlund2014.roxygen", unlink.target = TRUE)
system("rm -rf nordlund2014.roxygen/inst")
system("R CMD check nordlund2014.roxygen")
system("R CMD INSTALL nordlund2014.roxygen")

