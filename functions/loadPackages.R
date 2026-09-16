# This script load all the packages needed to fit the models

# Get a list of installed packages
installedPackages <- installed.packages()


# Now we bring int he rest of the necessary packages
necessaryPackages <- c("dismo",
                       "terra",
                       "effects", 
                       "MuMIn", 
                       "TR8", 
                       "sf", 
                       "dplyr", 
                       "inlabru",
                       "PointedSDMs",
                       "giscoR",
                       "bipartite",
                       "ggpubr",
                       "tidyterra",
                       "viridis",
                       "rgbif",
                       "purrr",
                       "tidyverse",
                       "tidyjson",
                       "xml2",
                       "rnaturalearth",
                       "vegan",
                       "dplyr"#,
                       #"LaplaceDemon"
                       )
uninstalledPackages <- necessaryPackages[!(necessaryPackages %in% row.names(installedPackages))]
if("giscoR" %in% uninstalledPackages){
  install.packages("giscoR",
                   repos = c("https://ropengov.r-universe.dev", "https://cloud.r-project.org")
  )
}

install.packages(uninstalledPackages)


library("dismo")
library("terra")
library(effects)
library(MuMIn)
library(TR8)
library(sf)
#library(intSDM)
 library(dplyr)
library(inlabru)
library(PointedSDMs)
library(bipartite)
library(ggpubr)
library(purrr)
 library(dplyr)
library(rgbif)
library(stringi)
library(stringr)
library(vegan)
library(viridis)
library(tidyterra)
require(tidyverse)
require(tidyjson)
require(xml2)
library(doParallel)
library(foreach)
library(sf)
library(rnaturalearth)
#library(LaplacesDemon)
