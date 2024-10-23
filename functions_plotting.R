

png("./Figures/Appendix/Temp_function.png", width = 2000, height = 2000, res=300)
temp = seq(0,45)
xi_temp = exp(-311.9924/ ((temp+ 273.15) - 227.13));

plot(temp, xi_temp, type = "l", main = "Loyd-Taylor",
     xlab = "Temperature", ylab = "Activity", col = "firebrick")

dev.off()


png("./Figures/Appendix/Moist_function.png", width = 2000, height = 2000, res=300)

### From Sierra SoilR
fW.Moyano <- function(theta, a = 3.11, b = 2.42) {
  a * theta - b * theta^2
}

# Example usage
th <- seq(0, 1, 0.01)
xi <- fW.Moyano(theta = th)

plot(th, xi, type = "l", main = "Moyano 2013",
     xlab =expression("Volumetric soil water content (cm"^3~"cm"^-3~")"),  ylab = "Activity", col = "dodgerblue3")
dev.off()


png("./Figures/Appendix/Season_function.png", width = 2000, height = 2000, res=300)

amplitude=0.1183339
day_year = seq(1:365)
sine_wave =  amplitude * cos((2 * pi / 365) * day_year + (2 * pi/ 365) * (176.3293 - 1) - pi / 2);

plot(day_year, sine_wave, type = "l", main = "Sine seasonality",
     xlab = "Day of the year", ylab = "Additive seasonality pattern", col = "darkorange")

dev.off()

ampl=0
peak=175


ampl * cos((2 * pi / 365) * day_year + (2 * pi / 365) * (peak - 1) - pi / 2)


