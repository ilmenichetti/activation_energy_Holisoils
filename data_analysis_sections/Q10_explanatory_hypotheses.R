




#the levels at which we need to collapse the datasets:
soil_data_aggregation_level = droplevels(interaction(soil_data$plot_ID, soil_data$site))
processed_data_filtered$plot_id


dim(post_bytreat_indA_Moy$q10)

# means of Q10
Q10_averages <- tapply(colMeans(post_bytreat_indA_Moy$q10), processed_data_filtered$plot_id, mean, na.rm = TRUE)


# means of all other variables
library(dplyr)

# Combine the soil_data with the aggregation level
soil_data <- soil_data %>%
  mutate(aggregation_level = soil_data_aggregation_level)  # Add aggregation level as a new column

# Aggregate by the aggregation level
aggregated_soil_data <- soil_data %>%
  group_by(aggregation_level) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))  # Calculate mean for numeric columns

aggregated_soil_data = cbind(aggregated_soil_data, Q10_averages)



# correlation matrix
names(aggregated_soil_data)


selected_columns <- c(
  "year", "No_trees", "No_trees_ha", "Basal_area", "Mean_DBH","bulk_den",
  "pH", "phosphorous", "Ntot", "Ctot",
  "C_stocks", "fungi_biomass", "bacteria_biomass", "actinobacteria_biomass",
  "gram_pos_biomass", "gram_neg_biomass", "total_biomass", "fungi_bact_rate",
  "beta_gluco", "acid_phos", "beta_xylo", "chitise",
  "cellobiohydrolase", "beta_galacto", "alpha_gluco", "lipase",
  "CO2_flux_hour", #"CO2_flux_norm",
  "T1_soil", "Soil_moist",
  "Q10_averages")

# Filter out columns that have all NA values or only one unique value
non_zero_sd_columns <- aggregated_soil_data[, selected_columns][, sapply(aggregated_soil_data[, selected_columns], function(x) {
  sd(x, na.rm = TRUE) > 0 && length(unique(na.omit(x))) > 1
})]

# Calculate correlation on columns with non-zero standard deviation
correlation_matrix <- cor(aggregated_soil_data[, selected_columns], use = "pairwise.complete.obs")
write.csv(correlation_matrix, "./Checks/Correlation_mat_Q10.csv")


# Plot the heatmap
png("./Checks/Correlation_mat_Q10.png", heigh=2000, width = 2000, res = 300)
heatmap(correlation_matrix,
        Rowv = NA,
        Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none",
        margins = c(6.5,6.5),
        main = "")
dev.off()
