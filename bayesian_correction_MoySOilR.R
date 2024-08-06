

library(rstan)
library(paletteer)
library(prettyGraphs)
library(randomForest)
library(dplyr)
library(lubridate)
library(readxl)
library(caret)
library(MASS)

#TODO DISCUSSION is there an interaction between texture and the TMS calibration?
#TODO Q10?

#TODO normalization of respiration by SOC when data are there

#DONE interpolated line of temperature by treatment (average)
#DONE add the same graph for respiration
#DONE try to make them side by side by site (3 plots)

#TODO plot the reduction functions in appendix

#TODO reorder the treatments with control first
#TODO treatments control thinning clearcut with different shades, and slah with 45 degrees lines


palette_treat <-  c(paletteer::paletteer_c("ggthemes::Orange-Gold", n=5),
                    paletteer::paletteer_c("ggthemes::Green", n=5),
                    paletteer::paletteer_c("ggthemes::Blue-Teal", n=5))

palette_plot <-  c(paletteer::paletteer_c("ggthemes::Orange-Gold", n=length(unique(data[data$site=="france",]$plot_ID))),
                   paletteer::paletteer_c("ggthemes::Green", n=length(unique(data[data$site=="romania",]$plot_ID))),
                   paletteer::paletteer_c("ggthemes::Blue-Teal", n=length(unique(data[data$site=="spain",]$plot_ID))))



edaphics <- read.csv("edaphics.csv")
#data <- read.csv("Database_co2_28_05_2024.csv")
data <- read_xlsx("Database_co2_28_05_2024.xlsx", sheet=1)
data$treatment <- as.factor(data$treatment)
data$site <- as.factor(data$site)
data$Soil_moist <- as.numeric(data$Soil_moist)/100
data$T1_soil <- as.numeric(data$T1_soil)
data$CO2_flux_hour <- as.numeric(data$CO2_flux_hour)
data$seq <- 1:dim(data)[1]


data$Date <- as.Date(data$Date, format = "%d.%m.%Y")


processed_data <- data.frame(ID = data$seq,
                             ID_old = data$ID,
                             CO2_flux_hour = data$CO2_flux_hour, T1_soil = data$T1_soil,
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
  filter(rowSums(is.na(select(., 3:8))) == 0)

dim(processed_data_filtered_preprocess)
dim(processed_data)

dropped_levels<-levels(processed_data_filtered_preprocess$plot_id)[!levels(processed_data_filtered_preprocess$plot_id) %in% levels(droplevels(processed_data_filtered_preprocess$plot_id))]
processed_data_filtered_preprocess$plot_id <- droplevels(processed_data_filtered_preprocess$plot_id)




########### site micrometeorology

png("./Figures/climatedata.png", height=2500, width = 2000, res=300)
par(mfrow=c(3,2), mar=c(2,5,1,1))

data = data.frame(Date = processed_data_filtered_preprocess$date,
                  Site= processed_data_filtered_preprocess$site,
                  Temperature = processed_data_filtered_preprocess$T1_soil,
                  Moisture = processed_data_filtered_preprocess$Soil_moist * 100)

sites <- levels(data$Site)
capitalized_sites <- paste0(toupper(substr(sites, 1, 1)), substr(sites, 2, nchar(sites)))
data$DateNumeric <- as.numeric(data$Date)
extended_date_range <- seq(min(data$DateNumeric) - 30, max(data$DateNumeric) + 30, length.out = 300)

for(i in 1:length(sites)){
subset<-data[data$Site== sites[i],]
# Remove NAs
data_clean <- subset %>% filter(!is.na(Temperature))
data_clean$DateNumeric <- as.numeric(data_clean$Date)

# Calculate the 2D kernel density estimate
dens_temp <- kde2d(data_clean$DateNumeric, data_clean$Temperature, n = 300, lims = c(range(extended_date_range), c(0,30)))
dens_moist <- kde2d(data_clean$DateNumeric, data_clean$Moisture, n = 300, lims = c(range(extended_date_range), c(0,60)))

custom_palette_temp <- colorRampPalette(c("white", "goldenrod1", "orangered", "red4"))(100)
custom_palette_moist <- colorRampPalette(c("white", "cadetblue2", "dodgerblue", "darkblue"))(100)

average_temperatures <- tapply(data_clean$Temperature, data_clean$DateNumeric, mean, na.rm = TRUE)
average_temperatures_df <- data.frame(DateNumeric = as.numeric(names(average_temperatures)), AverageTemperature = as.numeric(average_temperatures))
loess_model_temp <- loess(AverageTemperature ~ DateNumeric, data = average_temperatures_df, span = 0.25)
date_seq_temp <- seq(min(data_clean$DateNumeric), max(data_clean$DateNumeric), by = 1)
predicted_temperatures <- predict(loess_model_temp, newdata = data.frame(DateNumeric = date_seq_temp))
interpolated_df_temp <- data.frame(DateNumeric = date_seq, InterpolatedTemperature = predicted_temperatures)

average_moisture <- tapply(data_clean$Moisture, data_clean$DateNumeric, mean, na.rm = TRUE)
average_moisture_df <- data.frame(DateNumeric = as.numeric(names(average_moisture)), AverageMoisture = as.numeric(average_moisture))
loess_model_moist <- loess(AverageMoisture ~ DateNumeric, data = average_moisture_df, span = 0.25)
date_seq_moist <- seq(min(data_clean$DateNumeric), max(data_clean$DateNumeric), by = 1)
predicted_moisture <- predict(loess_model_moist, newdata = data.frame(DateNumeric = date_seq_moist))
interpolated_df_moist <- data.frame(DateNumeric = date_seq, InterpolatedMoisture = predicted_moisture)


# Plot temperature density over time with extended date range
image(dens_temp, xlab = "Date", ylab = "Temperature (°C)", main = "", axes = FALSE,
      xlim = range(extended_date_range), col = custom_palette_temp, ylim=c(0,30))
axis(1, at = seq(min(extended_date_range), max(extended_date_range), length.out = 10),
     labels = as.Date(seq(min(extended_date_range), max(extended_date_range), length.out = 10), origin = "1970-01-01"))
axis(2)
lines(interpolated_df_temp$DateNumeric, interpolated_df_temp$InterpolatedTemperature, col="firebrick")
legend("topright", paste("Average temperature", capitalized_sites[i]), pch=NA, bty="n")
if(i ==1){legend("topleft", "a)", pch=NA, bty="n")
} else if (i == 2){legend("topleft", "c)", pch=NA, bty="n")
    } else if (i == 3) {legend("topleft", "e)", pch=NA, bty="n")}
box()

# Plot moisture density over time with extended date range
image(dens_moist, xlab = "Date", ylab = "Moisture (% vol)", main = "", axes = FALSE,
      xlim = range(extended_date_range), col = custom_palette_moist, ylim=c(0,60))
axis(1, at = seq(min(extended_date_range), max(extended_date_range), length.out = 10),
     labels = as.Date(seq(min(extended_date_range), max(extended_date_range), length.out = 10), origin = "1970-01-01"))
axis(2)
lines(interpolated_df_moist$DateNumeric, interpolated_df_moist$InterpolatedMoisture, col="darkblue")
legend("topright", paste("Average moisture", capitalized_sites[i]), pch=NA, bty="n")
if(i ==1){legend("topleft", "b)", pch=NA, bty="n")
} else if (i == 2){legend("topleft", "d)", pch=NA, bty="n")
} else if (i == 3) {legend("topleft", "f)", pch=NA, bty="n")}

box()
}

dev.off()




png("./Figures/climatedata_means.png", height=2000, width = 2000, res=300)

layout( matrix(c(1,1,1,1,1,1,
                 2,2,2,2,2,2,
                 3,3,3,3,3,3,
                 4), ncol=1) )
par(mar=c(0,5,0,1))

data = data.frame(Date = processed_data_filtered_preprocess$date,
                  Site= processed_data_filtered_preprocess$site,
                  Temperature = processed_data_filtered_preprocess$T1_soil,
                  Moisture = processed_data_filtered_preprocess$Soil_moist * 100,
                  Treat = processed_data_filtered_preprocess$treatment,
                  Gen_treat = processed_data_filtered_preprocess$gen_treatment,
                  Resp = processed_data_filtered_preprocess$CO2_flux_hour)

sites <- levels(data$Site)
capitalized_sites <- paste0(toupper(substr(sites, 1, 1)), substr(sites, 2, nchar(sites)))
data$DateNumeric <- as.numeric(data$Date)
extended_date_range <- seq(min(data$DateNumeric) - 30, max(data$DateNumeric) + 30, length.out = 300)


# Load dplyr for data manipulation
library(dplyr)

  data_clean <- data %>% filter(!is.na(Temperature))
  data_clean$DateNumeric <- as.numeric(data_clean$Date)

  # Step 1: Compute average temperatures by DateNumeric and Treat
  average_temperatures <- data_clean %>%
    group_by(DateNumeric, Treat) %>%
    summarise(AverageTemperature = mean(Temperature, na.rm = TRUE))

  average_moistures <- data_clean %>%
    group_by(DateNumeric, Treat) %>%
    summarise(AverageMoisture = mean(Moisture, na.rm = TRUE))

  average_resp <- data_clean %>%
    group_by(DateNumeric, Treat) %>%
    summarise(AverageResp = mean(Resp, na.rm = TRUE))

treats = levels(average_temperatures$Treat)
treats_gen=c(1,2,3,4,5,
        1,2,3,4,5,
        1,2,3,4,5)

plot(average_temperatures[average_temperatures$Treat == treats[1],]$DateNumeric, pch=NA, xlab = "Date", ylab = "Temperature (°C)",
     average_temperatures[average_temperatures$Treat == treats[1],]$AverageTemperature, main = "", axes = FALSE,
     xlim = range(extended_date_range), col = custom_palette_temp, ylim=c(0,23))
axis(2)
box()
for(i in 1:length(treats)){
loess_model_temp <- loess(average_temperatures[average_temperatures$Treat == treats[i],]$AverageTemperature ~
                            average_temperatures[average_temperatures$Treat == treats[i],]$DateNumeric, span = 0.5)
date_seq_temp <- seq(min(average_temperatures[average_temperatures$Treat == treats[i],]$DateNumeric),
                     max(average_temperatures[average_temperatures$Treat == treats[i],]$DateNumeric), by = 1)
predicted_temperatures <- predict(loess_model_temp, newdata =  date_seq_temp)
lines(date_seq_temp, predicted_temperatures, col = palette_treat[i], lty=treats_gen[i], lwd=2)
}
legend("topright", paste("Average temperature"), pch=NA, bty="n")
legend("bottomleft", levels(data$Gen_treat), lty=1:5, bty="n")

plot(average_moistures[average_moistures$Treat == treats[1],]$DateNumeric, pch=NA, xlab = "Date",  ylab = "Moisture (% vol)",
     average_moistures[average_moistures$Treat == treats[1],]$AverageMoisture, main = "", axes = FALSE,
     xlim = range(extended_date_range), col = custom_palette_temp, ylim=c(0,50))
axis(2)
box()
for(i in 1:length(treats)){
  loess_model_temp <- loess(average_moistures[average_moistures$Treat == treats[i],]$AverageMoisture ~
                              average_moistures[average_moistures$Treat == treats[i],]$DateNumeric, span = 0.5)
  date_seq_temp <- seq(min(average_moistures[average_moistures$Treat == treats[i],]$DateNumeric),
                       max(average_moistures[average_moistures$Treat == treats[i],]$DateNumeric), by = 1)
  predicted_temperatures <- predict(loess_model_temp, newdata =  date_seq_temp)
  lines(date_seq_temp, predicted_temperatures, col = palette_treat[i], lty=treats_gen[i], lwd=2)
}
legend("topright", paste("Average moisture"), pch=NA, bty="n")

plot(average_resp[average_resp$Treat == treats[1],]$DateNumeric, pch=NA, xlab = "Date", ylab = expression(Respiration ~ (g ~ C ~ kg^{-1})),
     average_resp[average_resp$Treat == treats[1],]$AverageResp, main = "", axes = FALSE,
     xlim = range(extended_date_range), col = custom_palette_temp, ylim=c(0,0.85))
axis(1, at = seq(min(extended_date_range), max(extended_date_range), length.out = 10),
     labels = as.Date(seq(min(extended_date_range), max(extended_date_range), length.out = 10), origin = "1970-01-01"))
axis(2)
box()
for(i in 1:length(treats)){
  loess_model_temp <- loess(average_resp[average_resp$Treat == treats[i],]$AverageResp ~
                              average_resp[average_resp$Treat == treats[i],]$DateNumeric, span = 0.5)
  date_seq_temp <- seq(min(average_resp[average_resp$Treat == treats[i],]$DateNumeric),
                       max(average_resp[average_resp$Treat == treats[i],]$DateNumeric), by = 1)
  predicted_temperatures <- predict(loess_model_temp, newdata =  date_seq_temp)
  lines(date_seq_temp, predicted_temperatures, col = palette_treat[i], lty=treats_gen[i], lwd=2)
}
legend("topright", paste("Average respiration"), pch=NA, bty="n")


dev.off()



## stratified data split for validation
set.seed(123) # Set seed for reproducibility

# Create stratification group
strat_group <- interaction(processed_data_filtered$gen_treatment, processed_data_filtered$country)

# Create stratified partition
train_indices <- createDataPartition(strat_group, p = 0.8, list = FALSE)

# Split the data into training and validation sets
processed_data_filtered <- processed_data_filtered[train_indices, ]
validation_data <- processed_data_filtered[-train_indices, ]

dim(validation_data)
dim(processed_data_filtered)


#Regression Models
linear_model = lm(CO2_flux_hour ~ T1_soil * Soil_moist * treatment, processed_data_filtered)
summary(linear_model)$r.squared
plot(predict(linear_model, newdata = validation_data), validation_data$CO2_flux_hour)

rf_model = randomForest(CO2_flux_hour ~ T1_soil * Soil_moist * treatment, processed_data_filtered)
summary(rf_model)
rf_fit <- lm(predict(rf_model) ~ processed_data_filtered$CO2_flux_hour)
plot(predict(rf_model), processed_data_filtered$CO2_flux_hour)
summary(rf_fit)$r.squared


png("./Figures/fit_ML_validation.png", height = 2000, width = 4000, res=300)
par(mfrow=c(1,2))
plot(predict(rf_model, newdata = validation_data), validation_data$CO2_flux_hour, ylab="predicted", xlab="observed", main="RF model",
     col=palette_treat[as.numeric(validation_data$treatment)], pch=as.numeric(validation_data$treatment), xlim=c(0,3), ylim=c(0,3))
lm<-lm(validation_data$CO2_flux_hour ~ predict(rf_model, newdata = validation_data))
summary(lm)
abline(a=0, b=1, col="black", lty=2)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
legend("bottomleft", "a)", bty="n", pch=NA, col = NA)

plot(predict(linear_model, newdata = validation_data), validation_data$CO2_flux_hour, ylab="predicted", xlab="observed", main="Multiple linear model",
    col=palette_treat[as.numeric(validation_data$treatment)], pch=as.numeric(validation_data$treatment), xlim=c(0,3), ylim=c(0,3))
lm<-lm(validation_data$CO2_flux_hour ~ predict(linear_model, newdata = validation_data))
summary(lm)
abline(a=0, b=1, col="black", lty=2)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
legend("bottomleft", "b)", bty="n", pch=NA, col = NA)

dev.off()


varImp(rf_model)


library(parallel)

# Set the number of cores to use
options(mc.cores = parallel::detectCores())


# Data list for Stan
stan_data <- list(N = length(processed_data_filtered$CO2_flux_hour),
                  Tr = length(levels(processed_data_filtered$treatment)),
                  treatment = as.numeric(processed_data_filtered$treatment),
                  resp = processed_data_filtered$CO2_flux_hour,
                  temp = processed_data_filtered$T1_soil,
                  M = processed_data_filtered$Soil_moist,
                  Pl = length(levels(processed_data_filtered$plot_id)),
                  plot_id = as.numeric(processed_data_filtered$plot_id),
                  day_year = yday(as.Date(processed_data_filtered$date))
)

stan_data_valid <- list(N = length(validation_data$CO2_flux_hour),
                  Tr = length(levels(validation_data$treatment)),
                  treatment = as.numeric(validation_data$treatment),
                  resp = validation_data$CO2_flux_hour,
                  temp = validation_data$T1_soil,
                  M = validation_data$Soil_moist,
                  Pl = length(levels(validation_data$plot_id)),
                  plot_id = as.numeric(validation_data$plot_id),
                  day_year = yday(as.Date(validation_data$date))
)


# define the number of runs of the MCMC
N_RUNS = 5000

str(stan_data)

# fitting the model
start <- Sys.time()
fit_bytreat_indA_Moy_sin <- stan(file = 'temp_moist_season_model.stan',
                                 data = stan_data,
                                 iter = N_RUNS,
                                 chains = 10,
                                 control = list(adapt_delta = 0.99, max_treedepth = 15))
end <- Sys.time()
tot_time <- end - start
tot_time

#post_bytreat_indA_Moy <- extract(fit_bytreat_indA_Moy)
post_bytreat_indA_Moy <- extract(fit_bytreat_indA_Moy_sin)

dim(post_bytreat_indA_Moy$amplitude)
dim(post_bytreat_indA_Moy$peak_day)
dim(post_bytreat_indA_Moy$Ea)
dim(post_bytreat_indA_Moy$A)
dim(post_bytreat_indA_Moy$a)
dim(post_bytreat_indA_Moy$b)

mean(post_bytreat_indA_Moy$a)
mean(post_bytreat_indA_Moy$b)
mean(post_bytreat_indA_Moy$Ea)
mean(post_bytreat_indA_Moy$A)
mean(post_bytreat_indA_Moy$xi_temp)
mean(post_bytreat_indA_Moy$xi_moist)


posteriors_summary <- data.frame(
  mode_Ea = numeric(15),
  min_Ea = numeric(15),
  max_Ea = numeric(15),
  mode_peak = numeric(15),
  min_peak = numeric(15),
  max_peak = numeric(15),
  mode_ampl = numeric(15),
  min_ampl = numeric(15),
  max_ampl = numeric(15),
  mode_a = numeric(15),
  min_a = numeric(15),
  max_a = numeric(15),
  mode_b = numeric(15),
  min_b = numeric(15),
  max_b = numeric(15)
)

rownames(posteriors_summary) = levels(processed_data_filtered$treatment)


for(i in 1:15){
  dens = density(post_bytreat_indA_Moy$Ea[,i])
  posteriors_summary[i,1] = dens$x[which.max(dens$y)]
  posteriors_summary[i,2] = quantile(post_bytreat_indA_Moy$Ea[,i], 0.025)
  posteriors_summary[i,3] = quantile(post_bytreat_indA_Moy$Ea[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy$peak_day[,i])
  posteriors_summary[i,4] = dens$x[which.max(dens$y)]
  posteriors_summary[i,5] = quantile(post_bytreat_indA_Moy$peak_day[,i], 0.025)
  posteriors_summary[i,6] = quantile(post_bytreat_indA_Moy$peak_day[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy$amplitude[,i])
  posteriors_summary[i,7] = dens$x[which.max(dens$y)]
  posteriors_summary[i,8] = quantile(post_bytreat_indA_Moy$amplitude[,i], 0.025)
  posteriors_summary[i,9] = quantile(post_bytreat_indA_Moy$amplitude[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy$a[,i])
  posteriors_summary[i,10] = dens$x[which.max(dens$y)]
  posteriors_summary[i,11] = quantile(post_bytreat_indA_Moy$a[,i], 0.025)
  posteriors_summary[i,12] = quantile(post_bytreat_indA_Moy$a[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy$b[,i])
  posteriors_summary[i,13] = dens$x[which.max(dens$y)]
  posteriors_summary[i,14] = quantile(post_bytreat_indA_Moy$b[,i], 0.025)
  posteriors_summary[i,15] = quantile(post_bytreat_indA_Moy$b[,i], 0.975)
}

write.csv(posteriors_summary, file="./Tables/posteriors_summary.csv")



png("./Figures/posteriors_temp.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
E_0_densities_indA_Moy <- list()
E_0_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$Ea[,1]), xlim=c(range(post_bytreat_indA_Moy$Ea)[1]*0.9, range(post_bytreat_indA_Moy$Ea)[2]*1.1),
      ylim=c(0,0.07), xlab = expression('E'[0]), col=NA, main="Temperature scaling")
for(i in 1:stan_data$Tr){
  E_0_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$Ea[,i])
  E_0_box_indA_Moy[[i]] <- E_0_densities_indA_Moy[[i]]$x
  polygon(E_0_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.4), border = add.alpha(palette_treat[i],0.8))

  which_max_y  <- which.max(E_0_densities_indA_Moy[[i]]$y)
  text(E_0_densities_indA_Moy[[i]]$x[which_max_y], E_0_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}

A_densities_indA_Moy <- list()
A_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$A[,1]), xlim=c(range(post_bytreat_indA_Moy$A)[1]*0.9, range(post_bytreat_indA_Moy$A)[2]*1.1),
      ylim=c(0, 0.03), xlab = expression('A'), col=NA, main="Temperature scaling")
for(i in 1:stan_data$Pl){
  A_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$A[,i])
  A_box_indA_Moy[[i]] <- A_densities_indA_Moy[[i]]$x
  polygon(A_densities_indA_Moy[[i]], col=add.alpha(palette_plot[i],0.25), border = add.alpha(palette_plot[i],0.6))

  which_max_y  <- which.max(A_densities_indA_Moy[[i]]$y)
  text(A_densities_indA_Moy[[i]]$x[which_max_y], A_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$plot_id)[i], col=palette_plot[i], cex=0.8)
}
dev.off()



png("./Figures/posteriors_seasonality.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
peak_densities_indA_Moy <- list()
peak_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$peak_day[,1]), xlim=c(range(post_bytreat_indA_Moy$peak_day)[1]*0.9, range(post_bytreat_indA_Moy$peak_day)[2]*1.1),
      ylim=c(0,0.05), xlab = "peak_day", col=NA, main="Seasonality scaling")
for(i in 1:stan_data$Tr){
  peak_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$peak_day[,i])
  polygon(peak_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.4), border = add.alpha(palette_treat[i],0.8))
  peak_box_indA_Moy[[i]] <- peak_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(peak_densities_indA_Moy[[i]]$y)
  text(peak_densities_indA_Moy[[i]]$x[which_max_y], peak_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}

amplitude_densities_indA_Moy <- list()
amplitude_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$amplitude[,1]), xlim=c(range(post_bytreat_indA_Moy$amplitude)[1]*0.9, range(post_bytreat_indA_Moy$amplitude)[2]*1.1),
      ylim=c(0, 70), xlab = expression('amplitude'), col=NA, main="Seasonality scaling")
for(i in 1:stan_data$Tr){
  amplitude_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$amplitude[,i])
  polygon(amplitude_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.25), border = add.alpha(palette_treat[i],0.6))
  amplitude_box_indA_Moy[[i]] <- amplitude_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(amplitude_densities_indA_Moy[[i]]$y)
  text(amplitude_densities_indA_Moy[[i]]$x[which_max_y], amplitude_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}
dev.off()



png("./Figures/posteriors_moisture.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
a_densities_indA_Moy <- list()
a_box_indA_Moy <- list()
plot(density(post_bytreat_indA_Moy$a[,1]), xlim=c(range(post_bytreat_indA_Moy$a)[1]*0.9, range(post_bytreat_indA_Moy$a)[2]*1.1),
      ylim=c(0,4.95), xlab = "a", col=NA, main="Moisture scaling")
for(i in 1:stan_data$Tr){
  a_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$a[,i])
  polygon(a_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.4), border = add.alpha(palette_treat[i],0.8))
  a_box_indA_Moy[[i]] <- a_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(a_densities_indA_Moy[[i]]$y)
  text(a_densities_indA_Moy[[i]]$x[which_max_y], a_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}

b_densities_indA_Moy <- list()
b_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$b[,1]), xlim=c(range(post_bytreat_indA_Moy$b)[1]*0.9, range(post_bytreat_indA_Moy$b)[2]*1.1),
      ylim=c(0, 4.4), xlab = expression('b'), col=NA, main="Moisture scaling")
for(i in 1:stan_data$Tr){
  b_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$b[,i])
  polygon(b_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.25), border = add.alpha(palette_treat[i],0.6))
  b_box_indA_Moy[[i]] <- b_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(b_densities_indA_Moy[[i]]$y)
  text(b_densities_indA_Moy[[i]]$x[which_max_y], b_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}
dev.off()




names_treats <- c("CC no-slash France", "CC slash France", "control France", " thin no-slash France", " thin slash France",
"CC no-slash Romania", "CC slash Romania", "control Romania", " thin no-slash Romania"," thin slash Romania",
"CC no-slash Spain", "CC slash Spain", "control Spain", " thin no-slash Spain", " thin slash Spain")


png("./Figures/posteriors_boxplots.png", height = 3000, width = 3000, res= 300)
par(mar=c(10,3,1.5,0), mfrow=c(3,2))
# Set up the layout matrix
layout_matrix <- matrix(c(
  1, 1,
  1, 1,
  1, 1,
  2, 3,
  2, 3,
  4, 5,
  4, 5
), nrow = 7, byrow = TRUE)

# Set up the layout
layout(layout_matrix)
boxplot(E_0_box_indA_Moy, names = names_treats, las=2, col = palette_treat, main = expression(E[0]))
legend("topright", "(a)", pch=NA, bty="n")
boxplot(peak_box_indA_Moy, names = names_treats, las=2, col = palette_treat, main = "peak")
legend("topright", "(b)", pch=NA, bty="n")
boxplot(amplitude_box_indA_Moy, names = names_treats, las=2, col = palette_treat, main = "amplitude")
legend("topright", "(c)", pch=NA, bty="n")
boxplot(a_box_indA_Moy, names = names_treats, las=2, col = palette_treat, main = "a")
legend("topright", "(d)", pch=NA, bty="n")
boxplot(b_box_indA_Moy, names = names_treats, las=2, col = palette_treat, main = "b")
legend("topright", "(e)", pch=NA, bty="n")
dev.off()

png("./Figures/posteriors_boxplots_A.png", height = 1500, width = 3500, res= 300)
par(mar=c(6,3,1.5,0))
boxplot(A_box_indA_Moy, names = levels(processed_data_filtered$plot_id), las=2, col = palette_plot, main = "A", cex.axis =0.5)
dev.off()

# Rhat
max(summary(fit_bytreat_indA_Moy_sin)$summary[,"Rhat"], na.rm = T)



# Calculate mean and standard deviation of predictions for each data point
predicted_means_indA_Moy <- apply(post_bytreat_indA_Moy$model_resp, 2, mean)
predicted_sds_indA_Moy <- apply(post_bytreat_indA_Moy$model_resp, 2, sd)

png("./Figures/fit_model.png", height = 2000, width = 2000, res=300)
plot(predicted_means_indA_Moy, processed_data_filtered$CO2_flux_hour, ylab="observed", xlab="predicted",
     col=palette_treat[as.numeric(processed_data_filtered$treatment)],
     pch=as.numeric(processed_data_filtered$treatment), xlim=c(0,3), ylim=c(0,3))
lm<-lm(processed_data_filtered$CO2_flux_hour ~ predicted_means_indA_Moy)
summary(lm)
abline(a=0, b=1, col="black", lty=2)
#legend("topleft", levels(processed_data$site), bty="n", pch=1:3)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
#legend("bottomright", paste("mean sigma =", round(mean(post_bytreat_indA_Moy$sigma),3)), bty="n", pch=NA, col = NA)
dev.off()


# independent validations


dim(post_bytreat_indA_Moy$amplitude)
dim(post_bytreat_indA_Moy$peak_day)
dim(post_bytreat_indA_Moy$Ea)
dim(post_bytreat_indA_Moy$A)
dim(post_bytreat_indA_Moy$a)
dim(post_bytreat_indA_Moy$b)


str(post_bytreat_indA_Moy)

mean(post_bytreat_indA_Moy$amplitude)/
mean(stan_data$resp) * 100

library(dplyr)
averaged_data <- as.data.frame(stan_data) %>%
  group_by(treatment) %>%
  summarise(
    avg_resp = mean(resp, na.rm = TRUE),
    avg_temp = mean(temp, na.rm = TRUE),
    avg_M = mean(M, na.rm = TRUE)
  )

png("./Figures/Amplitude_ratio.png", width = 2000, height=2000, res=300)
par(mar=c(12,5,2,2))
barplot(colMeans(post_bytreat_indA_Moy$amplitude)/averaged_data$avg_resp*100,
        names.arg = names_treats, las=2, col = palette_treat, ylab= "mean amplitude (as % of mean respiration)", ylim=c(0,15))
box()
dev.off()


dim(post_bytreat_indA_Moy$res)
dim(post_bytreat_indA_Moy$Ea)
dim(post_bytreat_indA_Moy$A)

### decomposition of the variance






#load the model for forward runs
source("temp_moist_season_model.R")

modeled_validation <- temp_moist_season(temp = stan_data_valid$temp,
                                          M = stan_data_valid$M,
                                          day_year = stan_data_valid$day_year,
                                          a = colMeans(post_bytreat_indA_Moy$a),
                                          b = colMeans(post_bytreat_indA_Moy$b),
                                          Ea = colMeans(post_bytreat_indA_Moy$Ea),
                                          A = colMeans(post_bytreat_indA_Moy$A),
                                          peak_day = colMeans(post_bytreat_indA_Moy$peak_day),
                                          amplitude = colMeans(post_bytreat_indA_Moy$amplitude),
                                          treatment = stan_data_valid$treatment,
                                          plot_id = stan_data_valid$plot_id)


png("./Figures/fit_model_validation.png", height = 2000, width = 2000, res=300)
plot(modeled_validation, validation_data$CO2_flux_hour, ylab="observed", xlab="predicted",
     col=palette_treat[as.numeric(validation_data$treatment)],
     pch=as.numeric(validation_data$treatment), xlim=c(0,3), ylim=c(0,3))
lm<-lm(validation_data$CO2_flux_hour ~ modeled_validation)
summary(lm)
abline(a=0, b=1, col="black", lty=2)
#legend("topleft", levels(processed_data$site), bty="n", pch=1:3)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
#legend("bottomright", paste("mean sigma =", round(mean(post_bytreat_indA_Moy$sigma),3)), bty="n", pch=NA, col = NA)
dev.off()


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
boxplot(temp_effect ~ stan_data$treatment, names = names_treats, las=2, col = palette_treat, main = "Climate change vulnerability", ylab= "increase in mineralization", xlab="")
dev.off()
