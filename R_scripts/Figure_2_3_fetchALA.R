# Taylor Diagrams for ENM: Fetch Atlas of Living Australia data for a target species.
#
# This script requires the following packages to tbe installed:
#    galah:
#       Command to install: install.packages("galah")
#
#    processALA:
#       Pre-requisite: Installation of hte paclkage "remotes" from CRAN (i.e.
#                      install.packages("remotes") if it is not already installed
#       Command to install: remotes::install_github("peterbat1/processALA")
#
# NOTE: If you have just installed package galah you should run the following
# command in the R console before using it:
#
#      galah::galah_config(email = "xxx")
#
# where "xxx" is the email you wish to register with ALA to receive information
# about download issues.
#
# The package processALA is designed to make resolving taxonomic names,
# downloading occurrence data and implementing basic data cleaning a simple
# process. It relies heavily (but not exclusively) on the package galah to
# manage calls on ALA API functions.
#
# The parameter 'baseOutputFolder' in the processALA function fetchALAdata
# should be a folder within which the function will create a sub-folder using
# the name of each taxon in the paramater 'taxonList'. This species sub-folder
# will contain the four files resulting from the fetch and preliminary cleaning
# of occurrence data for each taxon.
#
# Peter D. Wilson
# Adjunct Fellow
# School of Natural Sciences
# Faculty of Science and Engineering
# Macquarie University, Sydney, Australia
#
# 2022-03-24 

library(processALA)

baseFolder <- "/a/path/to/a/folder"

targetTaxon <- "Argyrodendron trifoliolatum"
target_Taxon <- gsub(" ", "_", targetTaxon)

fetchALAdata(taxonList = targetTaxon,
             baseOutputPath = baseFolder)

# The data used for figures 2 & 3 is created by combining the specimen
# occurrence data with the human (incidental) observation data:
specData <- read.csv(paste0(baseFolder, "/", targetTaxon, "/", target_Taxon, "_herbariumRecords.csv"),
                     stringsAsFactors = FALSE)
humanData <- read.csv(paste0(baseFolder, "/", targetTaxon, "/", target_Taxon, "_humanObservations.csv"),
                      stringsAsFactors = FALSE)
comboData <- rbind(specData, humanData)

# Modify column names for latitude and longitude to make it neat and tidy - no
# other reason!
colnames(comboData) <- gsub("decimalL", "l", colnames(comboData))

# Remove records (rows) with missing coordinates: We only need to check one of
# latitude or longitude as if one is missing the other is too
badRows <- which(is.na(comboData$decimalLatitude))
if (length(badRows) > 0) comboData <- comboData[-badRows, ]

write.csv(comboData, 
          paste0(baseFolder, "/", targetTaxon, "/", target_Taxon, "_combined_data.csv"), 
          row.names = FALSE)

