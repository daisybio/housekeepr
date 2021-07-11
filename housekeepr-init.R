## environmental variables
# environment in which to run in.
# - possible values: "production", "development"
houseekpr_env <- Sys.getenv("HOUSEKEEPR_ENV")
if (houseekpr_env == "") houseekpr_env <- "development"
# external url for absolute links
housekeepr_external_url <- Sys.getenv("HOUSEKEEPR_EXTERNAL_URL")
# if (housekeepr_external_url == "") housekeepr_external_url <- "/"
# data directory for downloaded files
hkgDataDir <- Sys.getenv("HOUSEKEEPR_DATADIR")
if (hkgDataDir == "") hkgDataDir <- "hkg-data"
# HISTORIC
# mode: should this app instance be used for computation or only presentation? 
# - possible values: "computation", "presentation"
housekeepr_mode <- Sys.getenv("HOUSEKEEPR_MODE")
if (housekeepr_mode == "") housekeepr_mode <- "computation"

# moved to server.R
# if(housekeepr_mode == "computation") library(GEOquery)
library(pool)
library(plyr)
library(dplyr)
library(rentrez)
library(hash)
#library(shiny)
#library(shinyjs)
library(V8)
#library(shinycssloaders)
library(bsplus)
library(devtools)
library(ggthemes)
#library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
# added due to error because require(ensembldb)
library(ensembldb)
library(limma)
library(data.table)
library(Biobase)
library(ggplot2)
library(calibrate)
library(ggrepel)
library(matrixStats)
library(rbokeh)
library(biomaRt)
library(splitstackshape)
library(stringr)
library(annotate)
library(sva)

# resizable plots
library(shinyjqui)
library(stringi)

# execute withTimeout()
library(R.utils)
# fix R.utils maskings
setProgress <- shiny::setProgress
validate <- shiny::validate

# workaround for:
#Warning: Error in curl::curl_fetch_memory: Error in the HTTP2 framing layer
#52: curl::curl_fetch_memory
#51: request_fetch.write_memory
#49: request_perform
#48: httr::GET
#47: <Anonymous>
#45: entrez_summary
#44: <observer> [/media/chris/zfs/owncloud/Work/projects/house-keeping-genes/app3.R#857]
# 1: runApp
# httr::set_config(httr::config(http_version = 0L))
library(curl)
# chandle <- new_handle()
# handle_setopt(chandle, http_version = 0L)
# options(download.file.method="curl")

library(AnnotationHub)
#GLOBAL VARIABLES
setAnnotationHubOption('CACHE', 'myHub')
ah <- AnnotationHub(ask = FALSE)

