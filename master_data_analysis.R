

library(rstan)
library(paletteer)
library(prettyGraphs)
library(randomForest)
library(dplyr)
library(lubridate)
library(readxl)
library(caret)
library(MASS)
library(parallel)

set.seed(123)

#TODO Use basal area instead of DBH in the RF model

#TODO NORMALIZE THE TEMPERATURE SCENARIOS!!

#TODO: E_o and E_a are mixed up in the manuscript

#TODO plot SOC stocks by site

#TODO: correlation. between E_a and enzyme activities

#NOTE: CUE and aciivation energy: in reality they relate, modelwise no, it's a shortcoming of having a model not distinguishing for C quality

#TODO normalization of respiration by SOC when data are there

#TODO test variable E_a over time, remove seasonality scaling and use E_a by time

#TODO Enzyme relationships

#TODO convert volumetric water content to theta_s? When data come

#DONE interpolated line of temperature by treatment (average)
#DONE add the same graph for respiration
#DONE try to make them side by side by site (3 plots)

#DONE plot the reduction functions in appendix

#DONE reorder the treatments with control first
#DONE in the boxplots treatments control thinning clearcut with different shades, and slah with 45 degrees lines





# edaphics <- read.csv("edaphics.csv")
#data <- read.csv("Database_co2_28_05_2024.csv")
data <- read_xlsx("Database_co2_03_09_2024_format_Lorenzo.xlsx", sheet=1)
str(data)
soil_data <- read_xlsx("Soil_data_Lorenzo_03_09_2024_enzymes.xlsx", sheet=1)


data$treatment <- as.factor(data$treatment)
data$site <- as.factor(data$site)
data$Soil_moist <- as.numeric(data$Soil_moist)/100
data$T1_soil <- as.numeric(data$T1_soil)
data$CO2_flux_hour <- as.numeric(data$CO2_flux_hour)
data$seq <- 1:dim(data)[1]

soil_data

data$Date <- as.Date(data$Date, format = "%d.%m.%Y")
soil_data$date  <- as.Date(soil_data$date, format = "%d.%m.%Y")

#rearrange the levels putting control first
data$treatment <- factor(data$treatment, levels = c("control", "clear_cut_no_slash", "clear_cut_slash", "thinning_no_slash", "thinning_slash"))



#####################

palette_treat <-  c(paletteer::paletteer_c("ggthemes::Orange-Gold", n=5),
                    paletteer::paletteer_c("ggthemes::Green", n=5),
                    paletteer::paletteer_c("ggthemes::Blue-Teal", n=5))
palette_treat_simplified <- c(palette_treat[1], palette_treat[3], palette_treat[3], palette_treat[5], palette_treat[5],
                              palette_treat[6], palette_treat[8], palette_treat[8], palette_treat[10], palette_treat[10],
                              palette_treat[11], palette_treat[13], palette_treat[13], palette_treat[15], palette_treat[15])


palette_plot <-  c(paletteer::paletteer_c("ggthemes::Orange-Gold", n=length(unique(data[data$site=="france",]$plot_ID))),
                   paletteer::paletteer_c("ggthemes::Green", n=length(unique(data[data$site=="romania",]$plot_ID))),
                   paletteer::paletteer_c("ggthemes::Blue-Teal", n=length(unique(data[data$site=="spain",]$plot_ID))))





#loop to parse the SOC data and average respiration by SOC stocks

# Initialize the CO2_flux_norm column
data$CO2_flux_norm <- NA
max_Ctot = max(soil_data$Ctot)

# Loop through each row in the data frame
for (i in 1:nrow(data)) {
  local_ID <- data[i,]$unique_ID
  local_SOC <- mean(soil_data[soil_data$unique_ID == local_ID,]$Ctot)
  local_normalized_flux <- data[i,]$CO2_flux_hour * (local_SOC/max_Ctot)
  data[i,]$CO2_flux_norm <- local_normalized_flux
}



########## ML decomposition for enzyme activities

# Building the table by adding the closest measurtements from the flux sampling campaign

# Add four new empty columns to store the values
soil_data$CO2_flux_hour <- NA
soil_data$CO2_flux_norm <- NA
soil_data$T1_soil <- NA
soil_data$Soil_moist <- NA
soil_data$month <- NA

#looping throught the values
for (i in 1:nrow(soil_data)) {

local_ID <- soil_data[i,]$unique_ID
local_date <- soil_data[i,]$date

soil_data[i,]$month <- month(local_date, label = TRUE)

local_subset <- data[data$unique_ID == local_ID,]

#Finding the closest date in the local subset
dates <- local_subset$Date
dates <- dates[!is.na(dates)]
differences <- abs(dates - local_date)
closest_index <- which.min(differences)

soil_data[i,]$CO2_flux_hour = local_subset[closest_index,]$CO2_flux_hour
soil_data[i,]$CO2_flux_norm = local_subset[closest_index,]$CO2_flux_norm
soil_data[i,]$T1_soil = local_subset[closest_index,]$T1_soil
soil_data[i,]$Soil_moist = local_subset[closest_index,]$Soil_moist

}


names(soil_data)


# ************************************************************
# ************** ML Enzyme variance decomposition ******************
# ************************************************************
source("./data_analysis_sections/ML_variance_decomposition_enzymes.R")






processed_data <- data.frame(ID = data$seq,
                             ID_old = data$ID,
                             CO2_flux_hour = data$CO2_flux_hour,
                             CO2_flux_norm = data$CO2_flux_norm,
                             T1_soil = data$T1_soil,
                             Soil_moist = data$Soil_moist,
                             treatment = interaction(data$treatment,as.factor(data$site)),
                             site = data$site,
                             plot_id = droplevels(interaction(as.factor(data$plot_ID), as.factor(data$site))), #each plot ID is repeated for the sites
                             date = data$Date,
                             gen_treatment = data$treatment,
                             gen_dist = data$disturbance,
                             plot_id_bycountry = as.factor(data$plot_ID),
                             country =  as.factor(data$site))



# Filter out rows with NAs in the first 8 columns, ignore NAs in the last 4 columns
processed_data_filtered_preprocess <- processed_data %>%
  filter(rowSums(across(3:8, ~ is.na(.))) == 0)

dropped_levels<-levels(processed_data_filtered_preprocess$plot_id)[!levels(processed_data_filtered_preprocess$plot_id) %in% levels(droplevels(processed_data_filtered_preprocess$plot_id))]
processed_data_filtered_preprocess$plot_id <- droplevels(processed_data_filtered_preprocess$plot_id)


# ************************************************************
# ************** Site micrometeorology ******************
# ************************************************************
source("./data_analysis_sections/site_micrometeorology.R")




# ************************************************************
# ************** Model benchmarks ******************
# ************************************************************
source("./data_analysis_sections/ML_model_benchmarks.R")


# ************************************************************
# ************** Model benchmarks ******************
# ************************************************************
source("./data_analysis_sections/parametric_model_fitting.R")



## Working on the residuals

residuals_indA_moy <- predicted_means_indA_Moy - processed_data_filtered$CO2_flux_hour

resdata=data.frame(residuals_indA_moy, processed_data_filtered)
mean_res_bytreat <- resdata %>%
  group_by(plot_id) %>%
  summarise(
    mean = mean(residuals_indA_moy, na.rm = T),  # Calculate mean of Value1
    treat = unique(treatment),     # Calculate sum of Value2
    gen_treat = unique(gen_treatment),     # Calculate sum of Value2
    gen_dist = unique(gen_dist)     # Calculate sum of Value2
  )

meanA_bytreat = colMeans(post_bytreat_indA_Moy$A)



png("./Checks/mean_residuals_bytreat.png", height=1800, width = 2000, res = 300)
par(mar=c(12,4,2,2))
numeric_treat <- as.integer(mean_res_bytreat$treat)
boxplot(mean_res_bytreat$mean ~ mean_res_bytreat$treat, xlab="", ylab = "mean residuals",
     col=palette_treat[unique(numeric_treat)], las=2)
abline(h=0)
dev.off()




residuals_quantiles <- quantile(residuals_indA_moy, c(0.05, 0.95))
negative_biased <- (residuals_indA_moy>residuals_quantiles[2])
positive_biased <- (residuals_indA_moy<residuals_quantiles[1])


write.csv(data.frame(processed_data_filtered,
                     predictions = predicted_means_indA_Moy,
                     underpredicted = negative_biased,
                     overpredicted = positive_biased), "./Checks/predictions.csv")

positive_bias_data<-processed_data_filtered[positive_biased,]
negative_bias_data<-processed_data_filtered[negative_biased,]

png("./Checks/study_residuals.png", height=3000, width=4000, res=300)
par(mfrow=c(2,2))
plot(density(processed_data_filtered[!positive_biased,]$T1_soil), main="Temperature", ylim=c(0,0.15))
polygon(density(processed_data_filtered[!positive_biased,]$T1_soil), col = add.alpha("darkorange", 0.5))
polygon(density(positive_bias_data$T1_soil), col = add.alpha("firebrick1", 0.5))
legend("topright", c("Other points", "Overpredicted points"), bty="n", pch=16, col=c("darkorange", "firebrick1"))

plot(density(processed_data_filtered[!positive_biased,]$Soil_moist), main="Moisture")
polygon(density(processed_data_filtered[!positive_biased,]$Soil_moist), col = add.alpha("cadetblue2", 0.5))
polygon(density(positive_bias_data$Soil_moist), col = add.alpha("dodgerblue3", 0.5))
legend("topright", c("Other points", "Overpredicted points"), bty="n", pch=16, col=c("cadetblue2", "dodgerblue3"))

par(mar=c(13,4,1,1))
boxplot(residuals_indA_moy[positive_biased] ~ positive_bias_data$treatment, las=2, ylab="Residuals", xlab="", main="Residuals by treatment")

barplot(table(positive_bias_data$treatment), las=2, main="Frequency among overpredicted")
box()
dev.off()




### Studying the residuals


png("./Figures/residuals_vs_time.png", height = 1500, width = 2800, res = 300)
plot(as.Date(processed_data_filtered$date), residuals_indA_moy, main = "Lloyd-Taylor + Moyano",ylab="Residuals",
     col=palette_treat[as.numeric(processed_data_filtered$treatment)], pch=as.numeric(processed_data$treatment), xaxt = "n", xlab="")
abline(h=0)
axis.Date(1, at = seq(from = min(as.Date(processed_data_filtered$date)), to = max(as.Date(processed_data_filtered$date)), by = "month"), format = "%b %Y", las=2)
legend("bottomright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
dev.off()



png("./Figures/residuals_vs_climate.png", height=1500, width=3000, res=300)
par(mfrow=c(1,2))
# Set parameters
plot(processed_data_filtered$T1_soil, residuals_indA_moy, xlab="T1 soil",
     ylab="residuals", xaxt="n",  col=palette_treat[as.numeric(processed_data_filtered$treatment)][resample_vector],
     pch=as.numeric(processed_data$treatment))
abline(h=0)
axis.Date(1, at = seq(from = min(as.Date(processed_data_filtered$date)), to = max(as.Date(processed_data_filtered$date)), by = "month"), format = "%b %Y", las=2)
# legend("bottomright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)

plot(processed_data_filtered$Soil_moist, residuals_indA_moy,  xlab="Soil moisture",
     ylab="residuals", xaxt="n", col=palette_treat[as.numeric(processed_data_filtered$treatment)][resample_vector],
     pch=as.numeric(processed_data$treatment))
abline(h=0)
axis.Date(1, at = seq(from = min(as.Date(processed_data_filtered$date)), to = max(as.Date(processed_data_filtered$date)), by = "month"), format = "%b %Y", las=2)

dev.off()




### Testing residuals against soil and tree data, only Spain

soil_data_Spain <- read_xlsx("Soil_data_Spain.xlsx")
names(soil_data_Spain)[2] = "plot_id_bycountry"

resdata[resdata$country =="spain",]$plot_id_bycountry

merged_resdata_soil = merge(soil_data_Spain, resdata[resdata$country =="spain",], by ="plot_id_bycountry" )

plot(merged_resdata_soil$residuals_indA_moy, merged_resdata_soil$`No. trees/ha`)
plot(merged_resdata_soil$residuals_indA_moy, merged_resdata_soil$fungi_bact_rate)
plot(merged_resdata_soil$residuals_indA_moy, merged_resdata_soil$fungi_biomass)

predictor_list<-c("disturbance",
                  "No. trees/ha", "No. Quercus/ha","Mean_DBH","porosity_Robin","bulk_den","pH","phosphorous","Ntot",
                  "Ctot","fungi_biomass","bacteria_biomass","actinobacteria_biomass","gram_pos_biomass","gram_neg_biomass",
                  "total_biomass","fungi_bact_rate",  "residuals_indA_moy")
rf_residuals_dataset = as.data.frame(merged_resdata_soil[,predictor_list])
names(rf_residuals_dataset)[2:3] = c("tree_density", "Quercus_density")
rf_residuals_dataset = na.omit(rf_residuals_dataset)

# Train the random forest model on the full dataset
rf_model_res <- randomForest(residuals_indA_moy ~ ., data = rf_residuals_dataset)


png("./Checks/residuals_prediction_spain.png", height=2500, width = 1500, res= 300)
par(mfrow=c(2,1))
# Predicted residuals from the random forest model
predicted_residuals <- predict(rf_model_res)

# Plot predicted residuals against actual residuals
plot(predicted_residuals, rf_residuals_dataset$residuals_indA_moy,
     xlab = "Predicted Residuals", ylab = "Actual Residuals")
abline(0, 1, col = "red")  # Adding a 45-degree line for reference

# Optional: Calculate R-squared for the fit of predicted vs actual residuals
rf_fit_res <- lm(predicted_residuals ~ rf_residuals_dataset$residuals_indA_moy)
summary(rf_fit_res)$r.squared
legend("topleft", paste("R^2 =", round(summary(rf_fit_res)$r.squared,3)), bty="n", pch=NA, col = NA)


varImpPlot(rf_model_res)

dev.off()


processed_data_Spain_calib = merge(soil_data_Spain, processed_data_filtered[processed_data_filtered$country =="spain",], by ="plot_id_bycountry" )
processed_data_Spain_valid = merge(soil_data_Spain, validation_data[validation_data$country =="spain",], by ="plot_id_bycountry" )
normalized_resp_calib <- processed_data_Spain_calib$CO2_flux_hour/processed_data_Spain_calib$Ctot
normalized_resp_valid <- processed_data_Spain_valid$CO2_flux_hour/processed_data_Spain_valid$Ctot

processed_data_Spain_calib <- cbind(processed_data_Spain_calib, normalized_resp = normalized_resp_calib)
processed_data_Spain_valid <- cbind(processed_data_Spain_valid, normalized_resp = normalized_resp_valid)

names(processed_data_Spain_calib)[13] = "No.trees.ha"
names(processed_data_Spain_valid)[13] = "No.trees.ha"

processed_data_Spain_calib$No.trees.ha <- as.numeric(processed_data_Spain_calib$No.trees.ha)
processed_data_Spain_calib$Mean_DBH <- as.numeric(processed_data_Spain_calib$Mean_DBH)
processed_data_Spain_calib$fungi_bact_rate <- as.numeric(processed_data_Spain_calib$fungi_bact_rate)

processed_data_Spain_valid$No.trees.ha <- as.numeric(processed_data_Spain_valid$No.trees.ha)
processed_data_Spain_valid$Mean_DBH <- as.numeric(processed_data_Spain_valid$Mean_DBH)
processed_data_Spain_valid$fungi_bact_rate <- as.numeric(processed_data_Spain_valid$fungi_bact_rate)

processed_data_Spain_calib[processed_data_Spain_calib$disturbance == "clear_cut",]$Mean_DBH = 0
processed_data_Spain_calib[processed_data_Spain_calib$disturbance == "clear_cut",]$No.trees.ha = 0
processed_data_Spain_valid[processed_data_Spain_valid$disturbance == "clear_cut",]$Mean_DBH = 0
processed_data_Spain_valid[processed_data_Spain_valid$disturbance == "clear_cut",]$No.trees.ha = 0



predictor_list = c("No.trees.ha","Mean_DBH","bulk_den_Robin",
               "porosity_Robin","pH","phosphorous","Ntot",
               "Ctot","fungi_biomass","bacteria_biomass","actinobacteria_biomass", "gram_pos_biomass","gram_neg_biomass",
               "total_biomass","CO2_flux_hour","T1_soil",
               "Soil_moist","gen_treatment")



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
  CO2_flux_hour ~ .,
  data = na.omit(processed_data_Spain_calib[,predictor_list]),
  method = "rf",
  trControl = control,
  tuneGrid = tune_grid,
  ntree = 500 # or another fixed value
)

rf_model_spain <- randomForest(CO2_flux_hour ~ . , data = na.omit(processed_data_Spain_calib[,predictor_list]),
                               mtry = tuned_rf$bestTune$mtry)


png("./Figures/variance_decomposition.png", height = 3000, width = 1500, res = 300)
par(mfrow=c(2,1))
plot(predict(rf_model_spain, newdata = processed_data_Spain_valid), processed_data_Spain_valid$CO2_flux_hour ,ylab="predicted", xlab="observed", main="",
     col=palette_treat[as.numeric(processed_data_Spain_valid$treatment.y)], pch=as.numeric(processed_data_Spain_valid$treatment.y), xlim=c(0,3), ylim=c(0,3))
rf_model_spain_fit <- summary(lm(processed_data_Spain_valid$CO2_flux_hour ~ predict(rf_model_spain, newdata = processed_data_Spain_valid[,predictor_list])))$r.squared
legend("topleft", paste("R^2 =", round(rf_model_spain_fit,3)), bty="n", pch=NA, col = NA)
abline(a=0, b=1, lty=2)
varimp_spain <- varImp(rf_model_spain)
varimp_spain[order(varimp_spain$Overall),]
par(mar=c(1,10,2,2))
barplot(varimp_spain[order(varimp_spain$Overall),], las=2, names.arg = rownames(varimp_spain)[order(varimp_spain$Overall)], horiz = T)
dev.off()


png("./Checks/residuals_vs_pH.png", height = 1500, width = 2800, res = 300)
plot(rf_residuals_dataset$pH, rf_residuals_dataset$residuals_indA_moy, main = "Lloyd-Taylor + Moyano",ylab="Residuals", xlab="pH",
     col=palette_treat[as.numeric(processed_data$treatment)], pch=as.numeric(processed_data$treatment))
abline(h=0)
legend("bottomright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
dev.off()








###########################################################################
#               extrapolation
###########################################################################


modeled_extrapolations <- temp_moist_season(temp = stan_data$temp+5,
                                        M = stan_data$M+0.05,
                                        day_year = stan_data$day_year,
                                        a = colMeans(post_bytreat_indA_Moy$a),
                                        b = colMeans(post_bytreat_indA_Moy$b),
                                        Ea = colMeans(post_bytreat_indA_Moy$Ea),
                                        A = colMeans(post_bytreat_indA_Moy$A),
                                        peak_day = colMeans(post_bytreat_indA_Moy$peak_day),
                                        amplitude = colMeans(post_bytreat_indA_Moy$amplitude),
                                        treatment = stan_data$treatment,
                                        plot_id = stan_data$plot_id)
modeled_zero <- temp_moist_season(temp = stan_data$temp,
                                            M = stan_data$M,
                                            day_year = stan_data$day_year,
                                            a = colMeans(post_bytreat_indA_Moy$a),
                                            b = colMeans(post_bytreat_indA_Moy$b),
                                            Ea = colMeans(post_bytreat_indA_Moy$Ea),
                                            A = colMeans(post_bytreat_indA_Moy$A),
                                            peak_day = colMeans(post_bytreat_indA_Moy$peak_day),
                                            amplitude = colMeans(post_bytreat_indA_Moy$amplitude),
                                            treatment = stan_data$treatment,
                                            plot_id = stan_data$plot_id)
plot(modeled_extrapolations,
     modeled_zero)


temp_effect = modeled_extrapolations -  modeled_zero
png("./Figures/scenario_extrapolation.png", height=1800, width = 2000, res = 300)
par(mar=c(12,4,2,2))
bp <- boxplot(temp_effect ~ stan_data$treatment, names = names_treats, las=2, col = palette_treat_simplified, main = "Climate change vulnerability", ylab= "increase in mineralization", xlab="")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
dev.off()
