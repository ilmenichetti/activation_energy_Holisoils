
########### site micrometeorology

png("./Figures/Appendix/climatedata.png", height=2500, width = 2000, res=300)
par(mfrow=c(3,2), mar=c(2,5,1,1))

data_meteo = data.frame(Date = processed_data_filtered_preprocess$date,
                  Site= processed_data_filtered_preprocess$site,
                  Temperature = processed_data_filtered_preprocess$T1_soil,
                  Moisture = processed_data_filtered_preprocess$Soil_moist * 100)

sites <- levels(data_meteo$Site)
capitalized_sites <- paste0(toupper(substr(sites, 1, 1)), substr(sites, 2, nchar(sites)))
data_meteo$DateNumeric <- as.numeric(data_meteo$Date)
extended_date_range <- seq(min(data_meteo$DateNumeric) - 30, max(data_meteo$DateNumeric) + 30, length.out = 300)

for(i in 1:length(sites)){
  subset<-data_meteo[data_meteo$Site== sites[i],]
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
  loess_model_temp <- loess(AverageTemperature ~ DateNumeric, data = average_temperatures_df)
  date_seq_temp <- seq(min(data_clean$DateNumeric), max(data_clean$DateNumeric), by = 1)
  predicted_temperatures <- predict(loess_model_temp, newdata = data.frame(DateNumeric = date_seq_temp))
  interpolated_df_temp <- data.frame(DateNumeric = date_seq_temp, InterpolatedTemperature = predicted_temperatures)

  average_moisture <- tapply(data_clean$Moisture, data_clean$DateNumeric, mean, na.rm = TRUE)
  average_moisture_df <- data.frame(DateNumeric = as.numeric(names(average_moisture)), AverageMoisture = as.numeric(average_moisture))
  loess_model_moist <- loess(AverageMoisture ~ DateNumeric, data = average_moisture_df)
  date_seq_moist <- seq(min(data_clean$DateNumeric), max(data_clean$DateNumeric), by = 1)
  predicted_moisture <- predict(loess_model_moist, newdata = data.frame(DateNumeric = date_seq_moist))
  interpolated_df_moist <- data.frame(DateNumeric = date_seq_temp, InterpolatedMoisture = predicted_moisture)


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
     xlim = range(extended_date_range), col = custom_palette_temp, ylim=c(0,1))
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

