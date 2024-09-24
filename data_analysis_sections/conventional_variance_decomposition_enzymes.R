
library("multcompView")



# Function to calculate confidence intervals using quantiles
calc_ci <- function(data) {
  lower <- quantile(data, probs = 0.025)  # 2.5% quantile for lower bound
  upper <- quantile(data, probs = 0.975)  # 97.5% quantile for upper bound
  return(c(lower, upper))
}




plot_names = colnames(soil_data)[30:44]

enzyme_means_table = mat.or.vec( length(unique(interaction(soil_data$treatment, soil_data$site))), length(plot_names))
colnames(enzyme_means_table) = plot_names
rownames(enzyme_means_table) = levels(interaction(soil_data$treatment, soil_data$site))

png("./Figures/Appendix/enzymes_boxplots.png", height=4000, width = 5000, res = 300)
par(mar=c(14,5,1,1), mfrow=c(4,4))
for (i in 1:length(plot_names)) {

  name = plot_names[i]

  boxplot(soil_data[[name]] ~ interaction(soil_data$treatment, soil_data$site), xlab="", ylab="",
          col = palette_treat,
          las = 2,
          main = paste(name))

  # Perform ANOVA
  anova_model <- aov(soil_data[[name]] ~ interaction(soil_data$treatment, soil_data$site))

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
  anova_model <- aov(soil_data[[name]] ~ interaction(soil_data$treatment, soil_data$site))

  # Extract letter groupings using multcompLetters from the ANOVA model
  group_letters <- multcompLetters(TukeyHSD(anova_model)$`interaction(soil_data$treatment, soil_data$site)`[, "p adj"])$Letters

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

