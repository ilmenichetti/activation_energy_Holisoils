




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

# Step 2: Extract categorical variables as separate vectors
site_vector <- soil_data %>%
  group_by(aggregation_level) %>%
  summarise(site = first(site)) %>%
  pull(site)

treatment_vector <- soil_data %>%
  group_by(aggregation_level) %>%
  summarise(treatment = first(treatment)) %>%
  pull(treatment)


treatment_vector = treatment_vector[ rownames(aggregated_soil_data) != "19.romania"]
site_vector = site_vector[ rownames(aggregated_soil_data) != "19.romania"]

treatment_vector <- as.factor(treatment_vector)

# correlation matrix
names(aggregated_soil_data)


selected_columns <- c(
 "No_trees_ha", "Basal_area", "bulk_den", #"Mean_DBH",  removed becuase many NAs in ROmania "year",
  "pH",  "Ntot", "Ctot", "phosphorous", #removed becuase many NAs in Romania
  "C_stocks", "fungi_biomass", "bacteria_biomass", "actinobacteria_biomass",
  "gram_pos_biomass", "gram_neg_biomass", "total_biomass", "fungi_bact_rate",
  "beta_gluco", "acid_phos", "beta_xylo", "chitise",
  "cellobiohydrolase", "beta_galacto", "alpha_gluco", "lipase",
  "CO2_flux_norm",   #"CO2_flux_hour",
  "T1_soil", "Soil_moist",
  "Q10_averages")

### removing 19.romania because it misses soil temp and moist!
reduced_dataset = aggregated_soil_data[ rownames(aggregated_soil_data) != "19.romania", selected_columns]
reduced_dataset_sites = aggregated_soil_data[rownames(aggregated_soil_data) != "19.romania", "aggregation_level"]
negative_coords <- which( reduced_dataset< 0, arr.ind = TRUE)
reduced_dataset[negative_coords] = 0 # there are some negs in the enzymes
aggregated_soil_data$year

data_scaled <- scale(reduced_dataset)

# Calculate correlation on columns with non-zero standard deviation
correlation_matrix <- cor(data_scaled, use = "pairwise.complete.obs")
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

# ---------------------------------------------
#              Mixed models
# ---------------------------------------------

############### mixed model with C

library(lme4)
reduced_dataset_mixed = cbind(as.data.frame(data_scaled), treatment_vector)
mixed_model <- lmer(Q10_averages ~ No_trees_ha + fungi_bact_rate + C_stocks + T1_soil + Soil_moist + treatment_vector + (1 | site_vector), data = reduced_dataset_mixed)
summary(mixed_model)

mixed_model_lm <- summary(lm(predict(mixed_model) ~ reduced_dataset_mixed$Q10_averages))

# Calculate relative importance for fixed effects
library(MuMIn)

mixed_model_ml <- lmer(Q10_averages ~ No_trees_ha + fungi_bact_rate + C_stocks + T1_soil + Soil_moist + (1 | treatment_vector),
                       data = reduced_dataset_mixed, REML = FALSE, na.action = na.fail)

dredge_results <- dredge(mixed_model_ml)
rel_importance <- model.avg(get.models(dredge_results, subset = TRUE))
str(sw(rel_importance))

png("./Figures/mixed_model_relimp.png", height=1400, width = 1400, res=300)
par(mar = c(9, 4, 2, 1))  # Adjust 'mar' to control the space around each plot
barplot(as.numeric(sw(rel_importance)), names.arg = names(sw(rel_importance)),  col = "steelblue", ylab = "Relative importance", las=2)
legend("topright", legend = list(bquote(R^2 == .(round(mixed_model_lm$r.squared, 2)))), bty="n")
dev.off()



############### mixed model with N

mixed_model_N <- lmer(Q10_averages ~ No_trees_ha + fungi_bact_rate + Ntot + T1_soil + Soil_moist + treatment_vector + (1 | site_vector), data = reduced_dataset_mixed)
summary(mixed_model_N)

mixed_model_N_lm <-  summary(lm(predict(mixed_model_N) ~ reduced_dataset_mixed$Q10_averages))

# Calculate relative importance for fixed effects
mixed_model_ml_N <- lmer(Q10_averages ~ No_trees_ha + fungi_bact_rate + Ntot + T1_soil + Soil_moist + (1 | treatment_vector),
                       data = reduced_dataset_mixed, REML = FALSE, na.action = na.fail)

dredge_results_N <- dredge(mixed_model_ml_N)
rel_importance_N <- model.avg(get.models(dredge_results_N, subset = TRUE))

png("./Figures/mixed_model_relimp_N.png", height=1400, width = 1400, res=300)
par(mar = c(9, 4, 2, 1))  # Adjust 'mar' to control the space around each plot
barplot(as.numeric(sw(rel_importance_N)), names.arg = names(sw(rel_importance_N)),  col = "steelblue", ylab = "Relative importance", las=2)
legend("topright", legend = list(bquote(R^2 == .(round(mixed_model_N_lm$r.squared, 2)))), bty="n")
dev.off()





##########

##### PLS analysis
library(pls)
pls_model <- plsr(Q10_averages ~ ., data = as.data.frame(data_scaled), ncomp = 5)
summary(pls_model)
plot(pls_model, "validation", val.type = "R2")
pls_scores <- scores(pls_model)

explained_variance_response <- R2(pls_model, estimate = "train")

# Extract scores and loadings for the first two components
scores <- scores(pls_model)[, 1:2]
loadings <- loadings(pls_model)[, 1:2]

png("./Figures/Q10_PLS_biplot.png", height=2500, width = 1800, res=300)

layout(matrix(c(1, 1, 1, 2,2), nrow = 5, ncol = 1, byrow = TRUE))
par(mar = c(4, 4, 2, 1))  # Adjust 'mar' to control the space around each plot

# Define the biplot
plot(scores, xlab = "Component 1", ylab = "Component 2", main = "PLS Q10 ~ . ",  col = palette_site[as.factor(site_vector)], pch=as.numeric(as.factor(treatment_vector)))
# Add loadings as arrows
scaling=17
arrows(0, 0, loadings[, 1]*scaling, loadings[, 2]*scaling, col = "darkgrey", length = 0.1)
text(loadings[, 1]*scaling, loadings[, 2]*scaling, labels = rownames(loadings), col = "black", pos = 3, cex=0.8)
legend("topleft", levels(as.factor(treatment_vector)), pch=1:5, bty="n")
legend("topright", levels(as.factor(site_vector)), pch=1, bty="n", col=palette_site)
legend("bottomright", paste("Explained variance (2 comps) = ", round(explained_variance_response$val[2],2)), pch=NA, bty="n")

par(mar = c(12, 4, 2, 1))  # Adjust 'mar' to control the space around each plot
abs_coef <- (coef(pls_model))
barplot(abs_coef[order(abs_coef)], main = "", las = 2, col = "steelblue", ylab = "Regression coefficient", names.arg = rownames(abs_coef)[order(abs_coef)])
abline(h=0)
dev.off()

