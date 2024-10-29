

png("./Figures/SOC.png", height=1800, width = 2000, res = 300)
par(mar=c(13,5,1,1))
boxplot(soil_data$C_stocks ~ interaction(soil_data$site, soil_data$treatment),names = levels(interaction(soil_data$site, soil_data$treatment)),
        # col = palette_treat_simplified,
        main = "SOC content",
        ylab = expression(SOC~(Mg~ha^-1)),
        xlab="", las=2)
dev.off()



# # define the number of runs of the MCMC
# N_RUNS = 5000
#
# # Determine the number of cores available
# num_cores <- detectCores()
#
# # fitting the model (in multicore mode)
# start <- Sys.time()
# fit_bytreat_indA_indEa_Moy_sin <- stan(
#   file = 'temp_moist_season_model_ind_Ea.stan',
#   data = stan_data,
#   iter = N_RUNS,
#   chains = 4,
#   cores = num_cores-1,  # Use all available cores
#   control = list(adapt_delta = 0.99, max_treedepth = 15)
# )
# end <- Sys.time()
# tot_time <- end - start
# tot_time
#
# post_bytreat_indA_indEa_Moy <- extract(fit_bytreat_indA_indEa_Moy_sin)
#
# Ea_vector<- colMeans(post_bytreat_indA_indEa_Moy$Ea)
#


str(soil_data)
# Calculate the mean of numeric columns by 'treatment', ignoring NAs
mean_enzymes <- soil_data %>%
  group_by(unique_ID) %>%  # Group by the factor (e.g., treatment)
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))  # Apply mean to numeric columns, ignoring NAs
names(mean_enzymes)

plot(mean_enzymes$CO2_flux_norm, Ea_vector)


corr_vectors<-c()
rsquared_vectors<-c()
for(i in 1:(36-5)){
  corr_vectors[i]<- cor(mean_enzymes[,6:36][,i], Ea_vector)
  rsquared_vectors[i]<-summary(lm(unlist(mean_enzymes[,6:36][,i])~ Ea_vector))$r.squared
}
names(corr_vectors) = names(mean_enzymes[,6:36])
names(rsquared_vectors) = names(mean_enzymes[,6:36])

png("./Figures/Appendix/Ea_corr.png", height=1800, width = 3000, res = 300)
par(mar=c(10,4,2,2))
barplot(rsquared_vectors, las=2, ylab="R^2")
dev.off()
