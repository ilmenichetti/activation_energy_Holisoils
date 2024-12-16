
names(soil_data)

units = read.x

predictor_list = c("No_trees_ha", "Mean_DBH", "Basal_area", "treatment" ,
                   "pH","bulk_den", "phosphorous","Ntot","C_stocks","fungi_biomass",
                   "bacteria_biomass","actinobacteria_biomass","gram_pos_biomass","gram_neg_biomass","total_biomass","fungi_bact_rate","CO2_flux_norm",
                   "T1_soil", "Soil_moist" ,"month" )

# Create stratification group
soil_strat_group <- interaction(soil_data$treatment, soil_data$site)

# Create stratified partition
train_indices <- createDataPartition(soil_strat_group, p = 0.8, list = FALSE)

# Split the data into training and validation sets
soil_data_train <- soil_data[train_indices, ]
soil_data_valid <- soil_data[-train_indices, ]

rf_soil_data_train <- na.omit(soil_data_train[,predictor_list])
rf_soil_data_valid <- na.omit(soil_data_valid[,predictor_list])

rf_soil_data_train$treatment <- as.factor(rf_soil_data_train$treatment)
rf_soil_data_valid$treatment <- as.factor(rf_soil_data_valid$treatment)


rf_soil_data_train_T1 <- na.omit(soil_data_train[soil_data_train$time == "T1",predictor_list])
rf_soil_data_valid_T1 <- na.omit(soil_data_valid[soil_data_valid$time == "T1",predictor_list])



dim(rf_soil_data_train)
dim(rf_soil_data_valid)



# Define the grid of hyperparameters to search
tune_grid <- expand.grid(
  mtry = c(3, 4, 5, 6) # Adjust based on the number of predictors
)

# Define cross-validation method
control <- trainControl(
  method = "cv",
  number = 5,
  search = "grid"
)

# Perform the tuning
set.seed(123) # For reproducibility
tuned_rf <- train(
  CO2_flux_norm ~ .,
  data = rf_soil_data_train,
  method = "rf",
  trControl = control,
  tuneGrid = tune_grid,
  ntree = 500 # or another fixed value
)


rf_model_CO2_flux <- randomForest(CO2_flux_norm ~ . , data = rf_soil_data_train,
                                  mtry = tuned_rf$bestTune$mtry)


rf_model_CO2_flux_T1 <- randomForest(CO2_flux_norm ~ . , data = rf_soil_data_train_T1,
                                  mtry = tuned_rf$bestTune$mtry)
varImpPlot(rf_model_CO2_flux_T1)

plot(predict(rf_model_CO2_flux_T1, newdata = rf_soil_data_valid_T1), rf_soil_data_valid_T1$CO2_flux_norm)
plot(predict(rf_model_CO2_flux_T1), rf_soil_data_train_T1$CO2_flux_norm)
summary(lm(predict(rf_model_CO2_flux_T1) ~ rf_soil_data_train_T1$CO2_flux_norm))



range_plot = range(c(predict(rf_model_CO2_flux, newdata = rf_soil_data_valid),
                     rf_soil_data_valid$CO2_flux_norm))

png("./Figures/variance_decomposition_edaphics.png", height = 3000, width = 1500, res = 300)
par(mfrow=c(2,1))
# the object range_plot comes from running the ML_model_benchmark.R file first
plot(predict(rf_model_CO2_flux, newdata = rf_soil_data_valid), rf_soil_data_valid$CO2_flux_norm, ylab="predicted", xlab="observed", main="RF model",
     col=palette_treat[as.numeric(rf_soil_data_valid$treatment)], pch=as.numeric(rf_soil_data_valid$treatment), xlim=range_plot, ylim=range_plot)
lm<-lm(rf_soil_data_valid$CO2_flux_norm ~ predict(rf_model_CO2_flux, newdata = rf_soil_data_valid))
summary(lm)
abline(a=0, b=1, col="black", lty=2)
legend("topright", levels(rf_soil_data_valid$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)

varimp_model_CO2_flux <- varImp(rf_model_CO2_flux)
varimp_model_CO2_flux[order(varimp_model_CO2_flux$Overall),]
par(mar=c(3,10,2,2))
barplot(varimp_model_CO2_flux[order(varimp_model_CO2_flux$Overall),], las=2,
        names.arg = rownames(varimp_model_CO2_flux)[order(varimp_model_CO2_flux$Overall)], horiz = T)
dev.off()



