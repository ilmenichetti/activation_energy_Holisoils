
gc()
library(parallel)
library(lubridate)

# Set the number of cores to use
options(mc.cores = parallel::detectCores())

interactions_years = droplevels(interaction(processed_data_filtered$treatment, (year(processed_data_filtered$date)-2021)))
interactions_years_val = droplevels(interaction(validation_data$treatment, (year(validation_data$date)-2021)))

# Data# Data# Data list for Stan
stan_data <- list(N = length(processed_data_filtered$CO2_flux_norm),
                  Tr = length(levels(processed_data_filtered$treatment)),
                  treatment = as.numeric(processed_data_filtered$treatment),
                  resp = processed_data_filtered$CO2_flux_norm,
                  temp = processed_data_filtered$T1_soil,
                  M = processed_data_filtered$Soil_moist,
                  Pl = length(levels(processed_data_filtered$plot_id)),
                  plot_id = as.integer(processed_data_filtered$plot_id),
                  day_year = yday(as.Date(processed_data_filtered$date)),
                  elapsed_y = year(processed_data_filtered$date)-2021,
                  treat_year = as.numeric(interactions_years),
                  treat_year_N = length(levels(interactions_years))
)

str(stan_data)

stan_data_valid <- list(N = length(validation_data$CO2_flux_norm),
                        Tr = length(levels(validation_data$treatment)),
                        treatment = as.numeric(validation_data$treatment),
                        resp = validation_data$CO2_flux_norm,
                        temp = validation_data$T1_soil,
                        M = validation_data$Soil_moist,
                        Pl = length(levels(validation_data$plot_id)),
                        plot_id = as.integer(validation_data$plot_id),
                        day_year = yday(as.Date(validation_data$date)),
                        elapsed_y = year(validation_data$date)-2021,
                        treat_year = as.numeric(interactions_years_val),
                        treat_year_N = length(levels(interactions_years_val))
)

#quick check of the data...
str(stan_data)
boxplot(stan_data$resp ~ stan_data$treat_year)


# define the number of runs of the MCMC
N_RUNS = 5000

# Determine the number of cores available
num_cores <- detectCores()

# fitting the model (in multicore mode)
start <- Sys.time()
fit_bytreat_indA_Moy_sin_timevar <- stan(
  file = 'temp_moist_season_model_timevariant.stan',
  data = stan_data,
  iter = N_RUNS,
  chains = 5,
  cores = num_cores-1,  # Use all available cores
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
end <- Sys.time()
tot_time <- end - start
tot_time

post_bytreat_indA_Moy_timevar <- extract(fit_bytreat_indA_Moy_sin_timevar)

dim(post_bytreat_indA_Moy_timevar$amplitude)
dim(post_bytreat_indA_Moy_timevar$peak_day)
dim(post_bytreat_indA_Moy_timevar$Ea)
dim(post_bytreat_indA_Moy_timevar$A)
dim(post_bytreat_indA_Moy_timevar$a)
dim(post_bytreat_indA_Moy_timevar$b)

mean(post_bytreat_indA_Moy_timevar$a)
mean(post_bytreat_indA_Moy_timevar$b)
mean(post_bytreat_indA_Moy_timevar$Ea)
mean(post_bytreat_indA_Moy_timevar$A)
mean(post_bytreat_indA_Moy_timevar$xi_temp)
mean(post_bytreat_indA_Moy_timevar$xi_moist)


posteriors_summary_timevar <- data.frame(
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

rownames(posteriors_summary_timevar) = levels(processed_data_filtered$treatment)


for(i in 1:15){
  dens = density(post_bytreat_indA_Moy_timevar$Ea[,i])
  posteriors_summary_timevar[i,1] = dens$x[which.max(dens$y)]
  posteriors_summary_timevar[i,2] = quantile(post_bytreat_indA_Moy_timevar$Ea[,i], 0.025)
  posteriors_summary_timevar[i,3] = quantile(post_bytreat_indA_Moy_timevar$Ea[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy_timevar$peak_day[,i])
  posteriors_summary_timevar[i,4] = dens$x[which.max(dens$y)]
  posteriors_summary_timevar[i,5] = quantile(post_bytreat_indA_Moy_timevar$peak_day[,i], 0.025)
  posteriors_summary_timevar[i,6] = quantile(post_bytreat_indA_Moy_timevar$peak_day[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy_timevar$amplitude[,i])
  posteriors_summary_timevar[i,7] = dens$x[which.max(dens$y)]
  posteriors_summary_timevar[i,8] = quantile(post_bytreat_indA_Moy_timevar$amplitude[,i], 0.025)
  posteriors_summary_timevar[i,9] = quantile(post_bytreat_indA_Moy_timevar$amplitude[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy_timevar$a[,i])
  posteriors_summary_timevar[i,10] = dens$x[which.max(dens$y)]
  posteriors_summary_timevar[i,11] = quantile(post_bytreat_indA_Moy_timevar$a[,i], 0.025)
  posteriors_summary_timevar[i,12] = quantile(post_bytreat_indA_Moy_timevar$a[,i], 0.975)
}

for(i in 1:15){
  dens = density(post_bytreat_indA_Moy_timevar$b[,i])
  posteriors_summary_timevar[i,13] = dens$x[which.max(dens$y)]
  posteriors_summary_timevar[i,14] = quantile(post_bytreat_indA_Moy_timevar$b[,i], 0.025)
  posteriors_summary_timevar[i,15] = quantile(post_bytreat_indA_Moy_timevar$b[,i], 0.975)
}

write.csv(posteriors_summary_timevar, file="./Tables/posteriors_summary_timevar.csv")


treat_year_factor <- droplevels(interaction(processed_data_filtered$treatment, (year(processed_data_filtered$date)-2021)))

png("./Figures/Appendix/posteriors_temp_timevar.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
E_0_densities_indA_Moy_timevar <- list()
E_0_box_indA_Moy_timevar <- list()
plot( density(post_bytreat_indA_Moy_timevar$Ea[,1]), xlim=c(range(post_bytreat_indA_Moy_timevar$Ea)[1]*0.9, range(post_bytreat_indA_Moy_timevar$Ea)[2]*1.1),
      ylim=c(0,0.07), xlab = expression('E'[a]), col=NA, main="")
for(i in 1:stan_data$treat_year_N){
  E_0_densities_indA_Moy_timevar[[i]] <- density(post_bytreat_indA_Moy_timevar$Ea[,i])
  E_0_box_indA_Moy_timevar[[i]] <- E_0_densities_indA_Moy_timevar[[i]]$x
  polygon(E_0_densities_indA_Moy_timevar[[i]], col=add.alpha(palette_plot[i],0.4), border = add.alpha(palette_treat[i],0.8))

  which_max_y  <- which.max(E_0_densities_indA_Moy_timevar[[i]]$y)
  text(E_0_densities_indA_Moy_timevar[[i]]$x[which_max_y], E_0_densities_indA_Moy_timevar[[i]]$y[which_max_y]+0.002,
       treat_year_factor[i], col=palette_treat[i], cex=0.8)
}
dev.off()


png("./Figures/Appendix/posteriors_boxplots_Ea_timevariant.png", height = 1500, width = 3500, res= 300)
par(mar=c(8,3,1.5,0))
bp <- boxplot(E_0_box_indA_Moy_timevar,
              names = levels(treat_year_factor), las=2, main = "A", cex.axis =0.5)
dev.off()

# Rhat
max(summary(fit_bytreat_indA_Moy_sin_timevar)$summary[,"Rhat"], na.rm = T)





# Splitting the factor levels by the point (.)
split_factors <- strsplit(as.character(levels(treat_year_factor)), "\\.")
split_treat <- as.data.frame(do.call(rbind, split_factors))
colnames(split_treat) = c("treat", "site", "year")





### calculate mean and CI
# Initialize vectors to store results
means <- numeric(length(E_0_box_indA_Moy_timevar))
ci_lower <- numeric(length(E_0_box_indA_Moy_timevar))
ci_upper <- numeric(length(E_0_box_indA_Moy_timevar))

# Calculate means and 95% CI for each series
for (i in seq_along(E_0_box_indA_Moy_timevar)) {
  series <- E_0_box_indA_Moy_timevar[[i]]
  means[i] <- mean(series)
  stderr <- sd(series) / sqrt(length(series))
  ci <- qt(0.975, df = length(series) - 1) * stderr
  ci_lower[i] <- means[i] - ci
  ci_upper[i] <- means[i] + ci
}

# Combine results into a data frame
Ea_timevar <- data.frame(
  Mean = means,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  treat = split_treat$treat,
  site = split_treat$site,
  year = split_treat$year
)



#Plot 3 panels as if it was a time series
png("./Figures/timeseries_Ea_timevariant.png", height = 1300, width = 3500, res= 300)

sites <- c("france", "romania", "spain")
sites_plot_names=c("France", "Romania", "Spain")
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  # Set up 3 panels for 3 sites
color_index=0
# Loop over each site to plot in separate panels
for (j in 1:length(sites)) {

  site = sites[j]
  # Filter data for the current site
  data_site <- Ea_timevar[Ea_timevar$site == site, ]

  # Find the range for the y-axis
  yrange <- range(c(Ea_timevar$CI_Lower, Ea_timevar$CI_Upper))

  # Create an empty plot
  plot(data_site$year, data_site$Mean, type = "n", ylim = yrange, xlim=c(0,3),
       xlab = "Year", ylab = expression('E'[a]), main = paste(sites_plot_names[j]))

  # Loop over each treatment to add points, lines, and CI bands
  treatments <- unique(data_site$treat)
  for (i in seq_along(treatments)) {
    treat_data <- data_site[data_site$treat == treatments[i], ]

    # Add the confidence interval bands (polygon)
    polygon(c(treat_data$year, rev(treat_data$year)),
            c(treat_data$CI_Upper, rev(treat_data$CI_Lower)),
            col = add.alpha(palette_treat[i+color_index], 0.25), border = NA)

    # Add the points and lines
    points(treat_data$year, treat_data$Mean, col = palette_treat[i+color_index], pch = i)
    lines(treat_data$year, treat_data$Mean, col = palette_treat[i+color_index])
  }
  color_index = color_index + 4

  if(j==1){legend("topright", treatments, pch=1:5, bty="n")}
  legend("bottomleft", paste("(", LETTERS[j], ")", sep=""), pch=NA, bty="n")
}
dev.off()


#to leave memory for following tasks
rm(fit_bytreat_indA_Moy_sin_timevar)
rm(post_bytreat_indA_Moy_timevar)
gc()
