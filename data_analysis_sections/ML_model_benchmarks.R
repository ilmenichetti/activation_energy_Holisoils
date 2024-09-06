

## stratified data split for validation
set.seed(123) # Set seed for reproducibility

# Create stratification group
strat_group <- interaction(processed_data_filtered_preprocess$gen_treatment, processed_data_filtered_preprocess$country)

# Create stratified partition
train_indices <- createDataPartition(strat_group, p = 0.8, list = FALSE)

# Split the data into training and validation sets
processed_data_filtered <- processed_data_filtered_preprocess[train_indices, ]
validation_data <- processed_data_filtered_preprocess[-train_indices, ]

dim(validation_data)
dim(processed_data_filtered)


#Regression Models
linear_model = lm(CO2_flux_norm ~ T1_soil * Soil_moist * treatment, processed_data_filtered)
summary(linear_model)$r.squared
plot(predict(linear_model, newdata = validation_data), validation_data$CO2_flux_norm)

rf_model = randomForest(CO2_flux_norm ~ T1_soil * Soil_moist * treatment, processed_data_filtered)
summary(rf_model)
rf_fit <- lm(predict(rf_model) ~ processed_data_filtered$CO2_flux_norm)
plot(predict(rf_model), processed_data_filtered$CO2_flux_hour)
summary(rf_fit)$r.squared





png("./Figures/fit_ML_validation.png", height = 2000, width = 4000, res=300)

range_plot = range(c(predict(rf_model, newdata = validation_data),
                     validation_data$CO2_flux_norm,
                     predict(linear_model, newdata = validation_data)))

par(mfrow=c(1,2))
plot(predict(rf_model, newdata = validation_data), validation_data$CO2_flux_norm, ylab="predicted", xlab="observed", main="RF model",
     col=palette_treat[as.numeric(validation_data$treatment)], pch=as.numeric(validation_data$treatment), xlim=range_plot, ylim=range_plot)
lm<-lm(validation_data$CO2_flux_norm ~ predict(rf_model, newdata = validation_data))
summary(lm)
abline(a=0, b=1, col="black", lty=2)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
legend("bottomleft", "a)", bty="n", pch=NA, col = NA)

plot(predict(linear_model, newdata = validation_data), validation_data$CO2_flux_norm, ylab="predicted", xlab="observed", main="Multiple linear model",
     col=palette_treat[as.numeric(validation_data$treatment)], pch=as.numeric(validation_data$treatment), xlim=range_plot, ylim=range_plot)
lm<-lm(validation_data$CO2_flux_norm ~ predict(linear_model, newdata = validation_data))
summary(lm)
abline(a=0, b=1, col="black", lty=2)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
legend("bottomleft", "b)", bty="n", pch=NA, col = NA)

dev.off()


res_linear <- validation_data$CO2_flux_norm - predict(linear_model, newdata = validation_data)

png("./Figures/Checks/flinear_model_residuals.png", height = 2000, width = 4000, res=300)
boxplot(res_linear ~ validation_data$treatment, ylab = "residuals linear model")
dev.off()


varImp(rf_model)

