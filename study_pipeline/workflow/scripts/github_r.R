# Title     : r_github
# Objective : Install r packages that are not on conda
# Created by: Evangelos A. Dimopoulos"
# Created on: 08/04/2022


# allow for command line args
args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]

#load libs
library(devtools)

#define function
github_libs <- function(outfile) {

  # check if treeWAS has already been installed
  if ("treeWAS" %in% rownames(installed.packages())) {
    print("Package treeWAS has already been installed")
   } else {
    install_github("caitiecollins/treeWAS", build_vignettes = TRUE)
   }

   # check if hierfstat has already been installed
  if ("hierfstat" %in% rownames(installed.packages())) {
    print("Package hierfstat has already been installed")
   } else {
    install_github("jgx65/hierfstat")
   }

  if ("poppr" %in% rownames(installed.packages())) {
    print("Package poppr has already been installed")
   } else {
    install.packages('poppr', repos = "http://cran.us.r-project.org")
   }

  # create empty file when done
  file.create(outfile)
}

# run the install function
github_libs(out_file)