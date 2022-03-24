# Taylor Diagram study: Fit an ensemble of ENMs
#
#
# Peter D. Wilson
# Adjunct Fellow
# School of Natural Sciences
# Faculty of Science and Engineering
# Macquarie University, North Ryde, NSW, Australia 2109
#
# 2022-03-19 

library(raster)
library(hmeasure)
library(maxnet)

source("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/R-scripts/taylor_diagram_ggplot2_ver2.R")

baseModelFolder <- "/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/"

#################################################
# Load occurrence data for model fitting
#################################################
occ_data <- read.csv("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Test data/Occurrence_data/Argyrodendron_trifoliolatum.csv", stringsAsFactors = FALSE)

#################################################
# Prepare environmental data for model fitting
# and projection
#################################################
env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/Climate/CHELSA/Current/Taylor/",
                           "*.tif",
                           full.names = TRUE))

template_ras <- env_stack[[1]]

# Determine grid cell indices which are not NA; this is very handy for fast
# insertion of values into rasters
goodCells <- which(!is.na(template_ras[]))

# Find cell indices of occurrence records (and remove duplicate records in
# cells, and records falling outside the raster extent)
occCellInd <- na.omit(unique(raster::cellFromXY(template_ras, occ_data[, c("longitude", "latitude")])))
cat("Number of unique occurrence records:", length(occCellInd), "\n")

# Environment matrix for occurrence cells
occ_env <- raster::extract(env_stack, occCellInd)

# Set occurrence cells to NA so they will not be selected as background cells
template_ras[occCellInd] <- NA

# Determine the pool of cell indices available for background point selection
availCells <- which(!is.na(template_ras[]))

# Randomly select 10,000 cells for the background sample
set.seed(1953)
bkgCells <- sample(availCells, 10000)

# Generate background environment matrix
bkg_env <- raster::extract(env_stack, bkgCells)

# Make the training environment table
train_env <- data.frame(rbind(occ_env, bkg_env))

#################################################
# MaxEnt using maxnet
#################################################
maxnet_occ <- c(rep(1, nrow(occ_env)), rep(0, nrow(bkg_env)))

cat("Fitting maxnet model\n")
maxnet_model <- maxnet::maxnet(p = maxnet_occ,
                               data = train_env,
                               regmult = 1,
                               maxnet.formula(maxnet_occ, train_env, classes = "lpq"))

cat("Saving maxnet model object\n")
saveRDS(maxnet_model, file = paste0(baseModelFolder, "Argyrodendron_model/maxnet_Argyrodendron_trifoliolatum_", as.character(Sys.Date()), ".rds"))

# Current projection
cat("Projecting onto Current climate\n")
env_table <- raster::extract(env_stack, goodCells)
current_proj <- predict(maxnet_model, env_table, type = "cloglog")[, 1]

# RCP4.5 2050 projection
cat("Projecting onto RCP4.5 2050 climate\n")
env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/Climate/CHELSA/AR5/bioclim_2041-2060/Taylor/rcp45_2050_mean/",
                                      "*.tif",
                                      full.names = TRUE))
env_table <- raster::extract(env_stack, goodCells)
rcp45_2050_proj <- predict(maxnet_model, env_table, type = "cloglog")[, 1]

# RCP8.5 2050 projection
cat("Projecting onto RCP8.5 2050 climate\n")
env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/Climate/CHELSA/AR5/bioclim_2041-2060/Taylor/rcp85_2050_mean/",
                                      "*.tif",
                                      full.names = TRUE))
env_table <- raster::extract(env_stack, goodCells)
rcp85_2050_proj <- predict(maxnet_model, env_table, type = "cloglog")[, 1]

# RCP4.5 2070 projection
cat("Projecting onto RCP4.5 2070 climate\n")
env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/Climate/CHELSA/AR5/bioclim_2061-2080/Taylor/rcp45_2070_mean/",
                                      "*.tif",
                                      full.names = TRUE))
env_table <- raster::extract(env_stack, goodCells)
rcp45_2070_proj <- predict(maxnet_model, env_table, type = "cloglog")[, 1]

# RCP8.5 2070 projection
cat("Projecting onto RCP8.5 2070 climate\n")
env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/Climate/CHELSA/AR5/bioclim_2061-2080/Taylor/rcp85_2070_mean/",
                                      "*.tif",
                                      full.names = TRUE))
env_table <- raster::extract(env_stack, goodCells)
rcp85_2070_proj <- predict(maxnet_model, env_table, type = "cloglog")[, 1]

# Prepare for Taylor Diagram generation
cat("Produce and save Taylor Diagram\n")
occInd <- match(occCellInd, goodCells)
taylor_data <- data.frame(current_proj,
                         #lgm_proj,
                         rcp45_2050_proj,
                         rcp45_2070_proj,
                         rcp85_2050_proj,
                         rcp85_2070_proj)


# Taylor Diagram at training occurrence locations
# Figure 2 of the manuscript
ans <- taylor_diagram(taylor_data[occInd, ],
                      plot_type = "half",
                      show_labels = TRUE,
                      model_labels = c("Current", "RCP4.5 2050",
                                       "RCP4.5 2070", "RCP8.5 2050",
                                       "RCP8.5 2070"))

# Make a preview plot...if you dare ;)
# plot(ans)

ggsave(paste0("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/climate_change_Diagram_occ_locations.png"),
       ans,
       width = 12,
       height = 12,
       units = "cm",
       dpi = 150,
       bg = "white")

# Taylor Diagram over full extent
# Figure 3 of the manuscript
ans <- taylor_diagram(taylor_data,
                      plot_type = "half",
                      show_labels = TRUE,
                      model_labels = c("Current", "RCP4.5 2050",
                                       "RCP4.5 2070", "RCP8.5 2050",
                                       "RCP8.5 2070"))

# plot(ans)

ggsave(paste0("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/climate_change_Diagram_full_extent.png"),
       ans,
       width = 12,
       height = 12,
       units = "cm",
       dpi = 150,
       bg = "white")
