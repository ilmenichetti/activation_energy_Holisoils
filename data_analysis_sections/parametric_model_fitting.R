

library(parallel)

# Set the number of cores to use
options(mc.cores = parallel::detectCores())

#Q10 range, the two temperatures where to calculate Q10
#Q10_range = range(processed_data_filtered$T1_soil)
Q10_range = c(15,25)



#reordering the treatments for plotting uniformity
processed_data_filtered$treatment <- factor(processed_data_filtered$treatment, levels = c("control.france","thinning_slash.france","thinning_no_slash.france",
                                                               "clear_cut_slash.france","clear_cut_no_slash.france",
                                                               "control.romania","thinning_slash.romania", "thinning_no_slash.romania",
                                                               "clear_cut_slash.romania","clear_cut_no_slash.romania",
                                                               "control.spain","thinning_slash.spain", "thinning_no_slash.spain",
                                                               "clear_cut_slash.spain","clear_cut_no_slash.spain"))
validation_data$treatment <- factor(validation_data$treatment, levels = c("control.france","thinning_slash.france","thinning_no_slash.france",
                                                                                          "clear_cut_slash.france","clear_cut_no_slash.france",
                                                                                          "control.romania","thinning_slash.romania", "thinning_no_slash.romania",
                                                                                          "clear_cut_slash.romania","clear_cut_no_slash.romania",
                                                                                          "control.spain","thinning_slash.spain", "thinning_no_slash.spain",
                                                                                          "clear_cut_slash.spain","clear_cut_no_slash.spain"))

# Data list for Stan
stan_data <- list(N = length(processed_data_filtered$CO2_flux_norm),
                  Tr = length(levels(processed_data_filtered$treatment)),
                  treatment = as.numeric(processed_data_filtered$treatment),
                  treatment_txt = (processed_data_filtered$treatment),
                  resp = processed_data_filtered$CO2_flux_norm,
                  temp = processed_data_filtered$T1_soil,
                  M = processed_data_filtered$Soil_moist,
                  Pl = length(levels(processed_data_filtered$plot_id)),
                  plot_id = as.integer(processed_data_filtered$plot_id),
                  day_year = yday(as.Date(processed_data_filtered$date)),
                  Q10_range = Q10_range
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
                        Q10_range = Q10_range
)

#quick check of the data...
str(stan_data)
plot(density(processed_data_filtered$CO2_flux_norm))
polygon(density(validation_data$CO2_flux_norm))


# define the number of runs of the MCMC
N_RUNS = 5000

# Determine the number of cores available
num_cores <- detectCores()

# fitting the model (in multicore mode)
start <- Sys.time()
fit_bytreat_indA_Moy_sin <- stan(
  file = 'temp_moist_season_model.stan',
  data = stan_data,
  iter = N_RUNS,
  chains = 5,
  cores = num_cores-1,  # Use all available cores
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
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
mean(post_bytreat_indA_Moy$q10)
mean(post_bytreat_indA_Moy$xi_temp_high)
mean(post_bytreat_indA_Moy$model_resp_low)
mean(post_bytreat_indA_Moy$model_resp_high)

length(post_bytreat_indA_Moy$q10)
dim(post_bytreat_indA_Moy$q10)
which(is.infinite(post_bytreat_indA_Moy$q10))

i=111
post_bytreat_indA_Moy$xi_temp_high[i]
post_bytreat_indA_Moy$xi_temp_low[i]
post_bytreat_indA_Moy$model_resp_low[i]
post_bytreat_indA_Moy$model_resp_high[i]
post_bytreat_indA_Moy$q10[i]
post_bytreat_indA_Moy$sine_wave[i]

as.numeric(post_bytreat_indA_Moy$q10)
which(is.na((post_bytreat_indA_Moy$model_resp_high / post_bytreat_indA_Moy$model_resp_low)^(10 / (Q10_range[2] - Q10_range[1]))))

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



png("./Figures/Appendix/posteriors_temp.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
E_0_densities_indA_Moy <- list()
E_0_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$Ea[,1]), xlim=c(range(post_bytreat_indA_Moy$Ea)[1]*0.9, range(post_bytreat_indA_Moy$Ea)[2]*1.1),
      ylim=c(0,0.07), xlab = expression('E'[a]), col=NA, main="")
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
      ylim=c(0, 0.03), xlab = expression('A'), col=NA, main="")
for(i in 1:stan_data$Pl){
  A_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$A[,i])
  A_box_indA_Moy[[i]] <- A_densities_indA_Moy[[i]]$x
  polygon(A_densities_indA_Moy[[i]], col=add.alpha(palette_plot[i],0.25), border = add.alpha(palette_plot[i],0.6))

  which_max_y  <- which.max(A_densities_indA_Moy[[i]]$y)
  text(A_densities_indA_Moy[[i]]$x[which_max_y], A_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$plot_id)[i], col=palette_plot[i], cex=0.8)
}
dev.off()

#barplot(rep(1, length(palette_treat)), col = palette_treat, border = NA, space = 0)


png("./Figures/Appendix/posteriors_seasonality.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
peak_densities_indA_Moy <- list()
peak_box_indA_Moy <- list()
plot( density(post_bytreat_indA_Moy$peak_day[,1]), xlim=c(range(post_bytreat_indA_Moy$peak_day)[1]*0.9, range(post_bytreat_indA_Moy$peak_day)[2]*1.1),
      ylim=c(0,0.05), xlab = "peak_day", col=NA, main="")
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
      ylim=c(0, 70), xlab = expression('amplitude'), col=NA, main="")
for(i in 1:stan_data$Tr){
  amplitude_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$amplitude[,i])
  polygon(amplitude_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.25), border = add.alpha(palette_treat[i],0.6))
  amplitude_box_indA_Moy[[i]] <- amplitude_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(amplitude_densities_indA_Moy[[i]]$y)
  text(amplitude_densities_indA_Moy[[i]]$x[which_max_y], amplitude_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}
dev.off()



png("./Figures/Appendix/posteriors_moisture.png", height = 2000, width = 4000, res=350)
par(mfrow=c(1,2))
a_densities_indA_Moy <- list()
a_box_indA_Moy <- list()
plot(density(post_bytreat_indA_Moy$a[,1]), xlim=c(range(post_bytreat_indA_Moy$a)[1]*0.9, range(post_bytreat_indA_Moy$a)[2]*1.1),
     ylim=c(0,4.95), xlab = "a", col=NA, main="")
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
      ylim=c(0, 4.4), xlab = expression('b'), col=NA, main="")
for(i in 1:stan_data$Tr){
  b_densities_indA_Moy[[i]] <- density(post_bytreat_indA_Moy$b[,i])
  polygon(b_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.25), border = add.alpha(palette_treat[i],0.6))
  b_box_indA_Moy[[i]] <- b_densities_indA_Moy[[i]]$x

  which_max_y  <- which.max(b_densities_indA_Moy[[i]]$y)
  text(b_densities_indA_Moy[[i]]$x[which_max_y], b_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}
dev.off()


levels(processed_data_filtered$treatment)

names_treats <- c( "control France", "thin slash France"," thin no-slash France",  "CC slash France", "CC no-slash France",
                   "control Romania",  "thin slash Romania","thin no-slash Romania",  "CC slash Romania", "CC no-slash Romania",
                   "control Spain", "thin slash Spain",  " thin no-slash Spain",  "CC slash Spain","CC no-slash Spain" )


angle_palette = rep(c(180, 45, NA, 45, NA), 3)
density_palette = rep(c(20, 20, NA, 20, NA), 3)

# define a custom function to add the shading to the boxplots
add_shading_to_boxplot <- function(bp, density, angle, x_offset = 0.4) {
  # Loop through each box in the boxplot
  for (i in 1:ncol(bp$stats)) {
    # Check if density and angle are not NA
    if (!is.na(density[i]) && !is.na(angle[i])) {
      # Calculate the positions for shading
      xleft <- i - x_offset  # left x position for shading
      xright <- i + x_offset # right x position for shading
      ybottom <- bp$stats[2, i]  # lower hinge (25th percentile)
      ytop <- bp$stats[4, i]     # upper hinge (75th percentile)

      # Add shading using the rect function with angle and density
      rect(xleft, ybottom, xright, ytop,
           density = density[i], angle = angle[i],
           col = "black", border = NA)
    }
  }
}


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

bp <- boxplot(E_0_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = expression(E[a]))
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(a)", pch=NA, bty="n")
bp <- boxplot(peak_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = "peak")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(b)", pch=NA, bty="n")
bp <- boxplot(amplitude_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = "amplitude")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(c)", pch=NA, bty="n")
bp <- boxplot(a_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = "a")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(d)", pch=NA, bty="n")
bp <- boxplot(b_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = "b")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(e)", pch=NA, bty="n")
dev.off()

png("./Figures/Appendix/posteriors_boxplots_A.png", height = 1500, width = 3500, res= 300)
par(mar=c(6,3,1.5,0))
bp <- boxplot(A_box_indA_Moy, names = levels(processed_data_filtered$plot_id), las=2, col = palette_plot, main = "A", cex.axis =0.5)
dev.off()

# Rhat
# max(summary(fit_bytreat_indA_Moy_sin)$summary[,"Rhat"], na.rm = T)

bp <- boxplot(b_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified, main = "b")
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
legend("topright", "(e)", pch=NA, bty="n")

# Calculate mean and standard deviation of predictions for each data point
predicted_means_indA_Moy <- apply(post_bytreat_indA_Moy$model_resp, 2, mean)
predicted_sds_indA_Moy <- apply(post_bytreat_indA_Moy$model_resp, 2, sd)

length(predicted_means_indA_Moy)

png("./Figures/fit_model.png", height = 2000, width = 2000, res=300)
plot(predicted_means_indA_Moy, processed_data_filtered$CO2_flux_norm, ylab="observed", xlab="predicted",
     col=palette_treat[as.numeric(processed_data_filtered$treatment)],
     pch=as.numeric(processed_data_filtered$treatment), xlim=c(0,3), ylim=c(0,3))
lm<-lm(processed_data_filtered$CO2_flux_norm ~ predicted_means_indA_Moy)
summary(lm)
abline(a=0, b=1, col="black", lty=2)
#legend("topleft", levels(processed_data$site), bty="n", pch=1:3)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
#legend("bottomright", paste("mean sigma =", round(mean(post_bytreat_indA_Moy$sigma),3)), bty="n", pch=NA, col = NA)
dev.off()


png("./Figures/Appendix/posteriors_q10.png", height = 1500, width = 2500, res= 300)

str(post_bytreat_indA_Moy$q10)
length(stan_data$treatment)

q10_densities_indA_Moy <- list()
q10_box_indA_Moy <- list()

plot(density(na.omit(post_bytreat_indA_Moy$q10[,1])), xlim=c(0, 5),
      ylim=c(0,2), xlab = expression('q'[10]), col=NA, main="")
for(i in 1:stan_data$Tr){
  q10_densities_indA_Moy[[i]] <- density(rowMeans(post_bytreat_indA_Moy$q10[,stan_data$treatment == i], na.rm=T))
  q10_box_indA_Moy[[i]] <- q10_densities_indA_Moy[[i]]$x
  polygon(q10_densities_indA_Moy[[i]], col=add.alpha(palette_treat[i],0.4), border = add.alpha(palette_treat[i],0.8))

  which_max_y  <- which.max(q10_densities_indA_Moy[[i]]$y)
  text(q10_densities_indA_Moy[[i]]$x[which_max_y], q10_densities_indA_Moy[[i]]$y[which_max_y]+0.002,
       levels(processed_data_filtered$treatment)[i], col=palette_treat[i], cex=0.8)
}
dev.off()

png("./Figures/posteriors_boxplots_q10.png", height = 1500, width = 1500, res= 300)
par(mar=c(8.5,3,1.5,0))
bp <- boxplot(q10_box_indA_Moy, names = names_treats, las=2, col = palette_treat_simplified,
              main = bquote('Q'['10,(T=' * .(round(Q10_range[1])) * ',' * .(round(Q10_range[2])) * ')']),  cex.axis =0.8)
add_shading_to_boxplot(bp, density = density_palette, angle = angle_palette)
dev.off()


library(dplyr)

# Convert stan_data to a data frame, excluding Q10_range
stan_data_df <- as.data.frame(stan_data[-which(names(stan_data) == "Q10_range")])
# Group by treatment and calculate averages
averaged_data <- stan_data_df %>%
  group_by(treatment) %>%
  summarise(
    avg_resp = mean(resp, na.rm = TRUE),
    avg_temp = mean(temp, na.rm = TRUE),
    avg_M = mean(M, na.rm = TRUE)
  )

png("./Figures/Amplitude_ratio.png", width = 2000, height=2000, res=300)
par(mar=c(12,5,2,2))
bar_positions <-barplot(colMeans(post_bytreat_indA_Moy$amplitude)/averaged_data$avg_resp*100,
                        names.arg = names_treats, las=2, col = palette_treat_simplified, ylab= "mean amplitude (as % of mean respiration)", ylim=c(-40,40))
# add_shading_to_barplot(bar_positions = bar_positions, bar_heights = colMeans(post_bytreat_indA_Moy$amplitude)/averaged_data$avg_resp*100,
#                        density = density_palette, angle = angle_palette)
box()
dev.off()




# fully independent validation

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
# the object range_plot comes from running the ML_model_benchmark.R file first
plot(modeled_validation, validation_data$CO2_flux_norm, ylab="observed", xlab="predicted",
     col=palette_treat[as.numeric(validation_data$treatment)],
     pch=as.numeric(validation_data$treatment), xlim=range_plot, ylim=range_plot)
lm<-lm(validation_data$CO2_flux_norm ~ modeled_validation)
summary(lm)
abline(a=0, b=1, col="black", lty=2)
#legend("topleft", levels(processed_data$site), bty="n", pch=1:3)
legend("topright", levels(processed_data$treatment), bty="n", pch=1:15, col = palette_treat)
legend("bottomright", paste("R^2 =", round(summary(lm)$r.squared,3)), bty="n", pch=NA, col = NA)
#legend("bottomright", paste("mean sigma =", round(mean(post_bytreat_indA_Moy$sigma),3)), bty="n", pch=NA, col = NA)
dev.off()

