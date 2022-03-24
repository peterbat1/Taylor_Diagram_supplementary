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

library(disdat)
library(dismo)
library(raster)
library(hmeasure)
library(mgcv)
library(randomForest)
library(maxnet)

baseModelFolder <- "/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/"

target_group <- "ba"

target_species <- "nsw04"

#################################################
# Prepare full extent environmental data for
# model projection
#################################################

env_stack <- raster::stack(list.files("/home/peterw/Data_and_Projects/ENM_env_data/Standard_test_data/data/Environment/NSW/",
                                      "*.tif",
                                      full.names = TRUE))

template_ras <- env_stack[[1]]
template_ras[] <- NA

ne_nsw <- as.data.frame(env_stack)
goodCells <- which(!is.na(rowSums(ne_nsw)))

# show unique values in each variable
ans <- apply(ne_nsw, 2, function(x) {length(unique(x))})


#################################################
# Prepare occurrence data for model fitting
#################################################

nsw_data <- disdat::disPo("NSW")

# Preliminary review of the data suggests that bat species "nsw04" or "nsw06" may be good test subjects
# Try "nsw04"

spp_stuff <- which(nsw_data$group == target_group)

nsw_spp_xtab <- as.data.frame.matrix(table(nsw_data$spid[spp_stuff], nsw_data$siteid[spp_stuff]))
spp_prevalence <- rowSums(nsw_spp_xtab)
occ_data <- nsw_data[nsw_data$spid == target_species, ]


nsw_bg <- disdat::disBg("NSW")

test_env <- disEnv("NSW", group = target_group)
test_occ_all <- disPa("NSW", group = target_group)
test_occ <- test_occ_all[, target_species]

train_data <- rbind(occ_data, nsw_bg)


env_vars <- c("cti", "mi", "rainann", "raindq", "rugged", "soildepth", "soilfert", "solrad",  "tempann", "topo", "vegsys")   # , "disturb""soilfert", 


#################################################
# GAM
#################################################

n_pres <- nrow(occ_data)
n_bkg <- nrow(nsw_bg)

case_wt <- ifelse(train_data$occ == 1, 1, n_pres/n_bkg)

#model_formula <- occ ~ s(cti) + s(rainann) + s(rugged) + s(soildepth) + s(solrad) + s(tempann) + s(topo) # + s(disturb) + s(soilfert)
model_formula <- occ ~ s(cti) + s(rainann) + s(raindq) + s(rugged) + s(soildepth) + s(solrad) + s(tempann) + s(topo) # + s(soilfert) + s(vegsys)

gam_model <- mgcv::gam(formula = model_formula,
                       data = train_data,
                       family = binomial, #(link = "logit"),
                       weights = case_wt,
                       method = "REML")

saveRDS(gam_model, file = paste0(baseModelFolder, "ensemble_models/GAM_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".rds"))

gam_perf_train <- plogis(predict(gam_model, train_data[, env_vars], link = "response"))
gam_perf_test <- plogis(predict(gam_model, test_env[, env_vars], link = "response"))

# Model performance
gam_perf_train <- hmeasure::HMeasure(train_data$occ, as.vector(gam_perf_train))
gam_AUC_train <- gam_perf_train$metrics[, "AUC"]

gam_perf_test <- hmeasure::HMeasure(test_occ, as.vector(gam_perf_test))
gam_AUC_test <- gam_perf_test$metrics[, "AUC"]

# Project model onto full extent
pred_vals <- mgcv::predict.gam(gam_model, ne_nsw, link = "logit")
output_ras <- template_ras
output_ras[goodCells] <- as.vector(pred_vals)[goodCells] #mgcv::predict.gam(gam_model, ne_nsw)

raster::writeRaster(output_ras,
                    filename = paste0(baseModelFolder, "model_projections/GAM_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".tif"),
                    options = "COMPRESS=DEFLATE",
                    overwrite = TRUE)


#################################################
# Random Forest
#################################################

# Transform occ to factor to use "classification" mode RF fitting
rf_data <- train_data
rf_data$occ <- as.factor(rf_data$occ)

rf_data <- rf_data[, c("occ", env_vars)]

rf_model <- randomForest::randomForest(formula = occ ~ .,
                                       data = rf_data,
                                       ntree = 500)

saveRDS(rf_model, file = paste0(baseModelFolder, "ensemble_models/RF_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".rds"))

# Prediction at test locations
rf_pred <- predict(rf_model, test_env[, env_vars], type = "prob")

# Model performance
rf_perf_train <- hmeasure::HMeasure(as.numeric(as.character(rf_data$occ)), as.numeric(rf_model$predicted))
rf_AUC_train <- rf_perf_train$metrics[, "AUC"]

rf_perf_test <- hmeasure::HMeasure(test_occ, rf_pred[, "1"])
rf_AUC_test <- rf_perf_test$metrics[, "AUC"]

# Project model onto full extent
pred_vals <- predict(rf_model, ne_nsw, type = "prob")

output_ras <- template_ras
output_ras[goodCells] <- pred_vals[goodCells, "1"] #mgcv::predict.gam(gam_model, ne_nsw)

raster::writeRaster(output_ras,
                    filename = paste0(baseModelFolder, "model_projections/RF_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".tif"),
                    options = "COMPRESS=DEFLATE",
                    overwrite = TRUE)


#################################################
# Random Forest with down-sampling to address
# class imbalance
#################################################

num_pres_train <- sum(rf_data$occ == 1)
#num_bkg_train <- sum(rf_data$occ == 0)
sample_sizes <- c("0" = num_pres_train, "1" = num_pres_train)

rf_model_downsampled <- randomForest::randomForest(formula = occ ~ .,
                                                   data = rf_data,
                                                   ntree = 1000,
                                                   sampsize = sample_sizes,
                                                   replace = TRUE)

saveRDS(rf_model_downsampled,
        file = paste0(baseModelFolder, "ensemble_models/RF_downsampled_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".rds"))


rf_ds_pred_test <- predict(rf_model_downsampled, test_env[, env_vars], type = "prob")

# Model performance
rf_ds_perf_train <- hmeasure::HMeasure(as.numeric(as.character(rf_data$occ)), as.numeric(rf_model_downsampled$predicted))
rf_ds_AUC_train <- rf_ds_perf_train$metrics[, "AUC"]

rf_ds_perf_test <- hmeasure::HMeasure(test_occ, rf_ds_pred_test[, "1"])
rf_ds_AUC_test <- rf_ds_perf_test$metrics[, "AUC"]

# Project model onto full extent
pred_vals <- predict(rf_model_downsampled, ne_nsw, type = "prob")

output_ras <- template_ras
output_ras[goodCells] <- pred_vals[goodCells, "1"] #mgcv::predict.gam(gam_model, ne_nsw)

raster::writeRaster(output_ras,
                    filename = paste0(baseModelFolder, "model_projections/RF_downsampled_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".tif"),
                    options = "COMPRESS=DEFLATE",
                    overwrite = TRUE)


#################################################
# Boosted Regression Tree
#################################################

brt_model <- gbm.step(data = train_data,
                      gbm.x = 7:ncol(train_data),
                      gbm.y = 5,
                      family = "bernoulli",
                      tree.complexity = 5,
                      learning.rate = 0.001,
                      bag.fraction = 0.75,
                      max.trees = 10000,
                      n.trees = 50,
                      n.folds = 5,
                      site.weights = case_wt,
                      silent = TRUE)

saveRDS(brt_model, file = paste0(baseModelFolder, "ensemble_models/BRT_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".rds"))


brt_pred_test <- predict(brt_model, test_env[, 5:ncol(test_env)], n.trees = brt_model$gbm.call$best.tress, type = "response")

# Model performance
brt_perf_train <- hmeasure::HMeasure(train_data$occ, brt_model$fitted)
brt_AUC_train <- brt_perf_train$metrics[, "AUC"]

brt_perf_test <- hmeasure::HMeasure(test_occ, brt_pred_test)
brt_AUC_test <- brt_perf_test$metrics[, "AUC"]

# Project model onto full extent
pred_vals <- predict(brt_model, ne_nsw, n.trees = brt_model$gbm.call$best.tress, type = "response")

# Refresh output_ras
output_ras <- template_ras
output_ras[goodCells] <- as.vector(pred_vals)[goodCells] #mgcv::predict.gam(gam_model, ne_nsw)

raster::writeRaster(output_ras,
                    filename = paste0(baseModelFolder, "model_projections/BRT_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".tif"),
                    options = "COMPRESS=DEFLATE",
                    overwrite = TRUE)


#################################################
# MaxEnt using maxnet
#################################################

maxnet_pres <- train_data$occ
maxnet_env <- train_data[, env_vars]

maxnet_model <- maxnet::maxnet(p = maxnet_pres,
                               data = maxnet_env,
                               regmult = 1,
                               maxnet.formula(maxnet_pres, maxnet_env, classes = "lpq"))


saveRDS(maxnet_model, file = paste0(baseModelFolder, "ensemble_models/maxnet_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".rds"))

maxnet_perf_train <- predict(maxnet_model, train_data[, env_vars], type = "cloglog")
maxnet_perf_test <- predict(maxnet_model, test_env[, env_vars], type = "cloglog")

# Model performance
maxnet_perf_train <- hmeasure::HMeasure(train_data$occ, maxnet_perf_train)
maxnet_AUC_train <- maxnet_perf_train$metrics[, "AUC"]

maxnet_perf_test <- hmeasure::HMeasure(test_occ, maxnet_perf_test)
maxnet_AUC_test <- maxnet_perf_test$metrics[, "AUC"]


# Project model onto full extent
pred_vals <- predict(maxnet_model, ne_nsw[, env_vars], type = "cloglog")

# Refresh otuput_ras
output_ras <- template_ras

output_ras[goodCells] <- pred_vals[, 1]

raster::writeRaster(output_ras,
                    filename = paste0(baseModelFolder, "model_projections/maxnet_NSW_", target_group,"_", target_species, "_", as.character(Sys.Date()), ".tif"),
                    options = "COMPRESS=DEFLATE",
                    overwrite = TRUE)


#################################################
# Make a summary performance table
#################################################

perf_table <- data.frame(model = c("GAM", "RF", "RF_downsampled", "BRT", "maxnet"),
                         AUC_train = c(gam_AUC_train, rf_AUC_train, rf_ds_AUC_train, brt_AUC_train, maxnet_AUC_train),
                         AUC_test = c(gam_AUC_test, rf_AUC_test, rf_ds_AUC_test, brt_AUC_test, maxnet_AUC_test))

write.csv(perf_table, paste0("/home/peterw/Data_and_Projects/Personal Projects/Taylor Diagrams and ENMs/Results/model_perf_table_", target_group, "_", target_species,".csv"), row.names = FALSE)



