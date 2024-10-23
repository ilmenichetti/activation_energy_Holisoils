
library("multcompView")

predictor_list

# Function to calculate confidence intervals using quantiles
calc_ci <- function(data) {
  lower <- quantile(data, probs = 0.025)  # 2.5% quantile for lower bound
  upper <- quantile(data, probs = 0.975)  # 97.5% quantile for upper bound
  return(c(lower, upper))
}




plot_names = colnames(soil_data)[30:44]
# reorder the plots
interaction_levels <- with(soil_data, interaction(treatment, site))
reordered_interaction <- factor(interaction_levels, levels = c("control.france","thinning_slash.france","thinning_no_slash.france",
                                                               "clear_cut_slash.france","clear_cut_no_slash.france",
                                                               "control.romania","thinning_slash.romania", "thinning_no_slash.romania",
                                                               "clear_cut_slash.romania","clear_cut_no_slash.romania",
                                                               "control.spain","thinning_slash.spain", "thinning_no_slash.spain",
                                                               "clear_cut_slash.spain","clear_cut_no_slash.spain"
))

enzyme_means_table = mat.or.vec( length(unique(interaction(soil_data$treatment, soil_data$site))), length(plot_names))
colnames(enzyme_means_table) = plot_names
rownames(enzyme_means_table) = levels(reordered_interaction)

png("./Figures/enzymes_boxplots.png", height=4000, width = 5000, res = 300)
par(mar=c(14,5,1,1), mfrow=c(4,4))
for (i in 1:length(plot_names)) {

  name = plot_names[i]


  boxplot(soil_data[[name]] ~ reordered_interaction, xlab="", ylab="",
          col = palette_treat,
          las = 2,
          main = paste(name))

  # Perform ANOVA
  anova_model <- aov(soil_data[[name]] ~ reordered_interaction)

  # Get the means and confidence intervals for each group
  means_table <- aggregate(soil_data[[name]],
                           by = list(interaction(soil_data$treatment, soil_data$site)),
                           FUN = mean)
  colnames(means_table) <- c("Group", "Mean")

  # Apply quantile-based CI calculation for each group
  ci_table <- aggregate(soil_data[[name]],
                        by = list(interaction(soil_data$treatment, soil_data$site)),
                        FUN = calc_ci)

  # Extract lower and upper CIs using apply
  lower_cis <- apply(ci_table$x, 1, function(row) row[1])  # Extract lower bound (2.5%)
  upper_cis <- apply(ci_table$x, 1, function(row) row[2])  # Extract upper bound (97.5%)
  ci_table <- data.frame(Group = ci_table$Group,
                         Lower = lower_cis,
                         Upper = upper_cis)

  # Perform ANOVA again (to match grouping)
  anova_model <- aov(soil_data[[name]] ~ reordered_interaction)

  # Extract letter groupings using multcompLetters from the ANOVA model
  group_letters <- multcompLetters(TukeyHSD(anova_model)$reordered_interaction[, "p adj"])$Letters

  # Merge means and CI into one table
  combined_table <- merge(means_table, ci_table, by = "Group")

  # Add the groupings to the combined table
  combined_table$Grouping <- group_letters

  # Create a final column combining Mean (Lower, Upper) ^Grouping
  combined_table$Result <- paste0(
    round(combined_table$Mean, 2),
    " (", round(combined_table$Lower, 2), "-", round(combined_table$Upper, 2), ")",
    "<sup>", combined_table$Grouping, "</sup>"
  )

  enzyme_means_table[, i] <- combined_table$Result

}
dev.off()

write.csv(enzyme_means_table, "./Tables/Enzyme_means.csv")



calc_ci <- function(x) {
  if (length(na.omit(x)) < 2) return(c(NA, NA))  # Handle cases where not enough non-NA data points
  mean_x <- mean(x, na.rm = TRUE)
  se_x <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))  # Calculate standard error, removing NAs
  ci <- qt(0.975, df = length(na.omit(x)) - 1) * se_x     # 95% CI
  return(c(mean_x - ci, mean_x + ci))  # Return lower and upper bounds
}


plot_names_edaphic = c( "No_trees_ha","Basal_area", "Mean_DBH", "bulk_den",
                       "pH","phosphorous","Ntot","Ctot","C_stocks")


edaphic_means_table = mat.or.vec( length(unique(interaction(soil_data$treatment, soil_data$site))), length(plot_names_edaphic))
colnames(edaphic_means_table) = plot_names_edaphic
rownames(edaphic_means_table) = levels(reordered_interaction)

png("./Figures/edaphic_boxplots.png", height=3000, width = 4000, res = 300)
par(mar=c(14,5,1,1), mfrow=c(3,3))
for (i in 1:length(plot_names_edaphic)) {

  name = plot_names_edaphic[i]

  boxplot(soil_data[[name]] ~ reordered_interaction, xlab="", ylab="",
          col = palette_treat,
          las = 2,
          main = paste(name))

  # Perform ANOVA
  anova_model <- aov(soil_data[[name]] ~ reordered_interaction)

  # Get the means and confidence intervals for each group
  means_table <- aggregate(soil_data[[name]],
                           by = list(interaction(soil_data$treatment, soil_data$site)),
                           FUN = mean)
  colnames(means_table) <- c("Group", "Mean")

  # Apply quantile-based CI calculation for each group
  ci_table <- aggregate(soil_data[[name]],
                        by = list(interaction(soil_data$treatment, soil_data$site)),
                        FUN = function(x) if (all(is.na(x))) c(NA, NA) else calc_ci(x))

  colnames(ci_table) = c("Group", "CI")

  # Extract lower and upper CIs using apply
  lower_cis <- ci_table$CI[,1]  # Extract lower bound (2.5%)
  upper_cis <- ci_table$CI[,2]  # Extract upper bound (97.5%)
  ci_table <- data.frame(Group = ci_table$Group,
                         Lower = lower_cis,
                         Upper = upper_cis)

  # Perform ANOVA again (to match grouping)
  anova_model <- aov(soil_data[[name]] ~ reordered_interaction)

  # Extract letter groupings using multcompLetters from the ANOVA model
  group_letters_full <- rep(NA, length(levels(reordered_interaction)))
  names(group_letters_full) <- levels(reordered_interaction)
  group_letters_full[names(group_letters)] <- group_letters
  combined_table$Grouping <- group_letters_full[combined_table$Group]

    # Merge means and CI into one table
  combined_table <- merge(means_table, ci_table, by = "Group")


  # Create a final column combining Mean (Lower, Upper) ^Grouping
  combined_table$Result <- paste0(
    round(combined_table$Mean, 2),
    " (", round(combined_table$Lower, 2), "-", round(combined_table$Upper, 2), ")",
    "<sup>", combined_table$Grouping, "</sup>"
  )

  edaphic_means_table[, i] <- combined_table$Result

}
dev.off()



write.csv(edaphic_means_table, "./Tables/Edaphics_means.csv")

