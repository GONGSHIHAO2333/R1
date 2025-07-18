
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("randomForest")
# install.packages("Hmisc")
# install.packages("showtext")

library(ggplot2)
library(dplyr)
library(randomForest)
library(Hmisc)


alpha_diversity_path <- "alpha_diversity_results.csv"
env_factors_path <- "environment.csv"
output_dir <- "************"


alpha_diversity <- read.csv(alpha_diversity_path, row.names = 1)
env_factors <- read.csv(env_factors_path, row.names = 1)


alpha_diversity <- cbind(env_factors, alpha_diversity)


alpha_diversity <- na.omit(alpha_diversity)

env_factors <- c(
  "pH", "Salinity", "DOC", "DO", "Temperature", "NO3", "PO4", "Total_Antibiotic_Stress"
)
diversity_measures <- unique(c(
  "Anammox", "ANR", "Nitrifictaion", "Denitrification", 
  "DNR", "Nitrogen_Fixation","BAB", "Total_Pathway"
))


print(diversity_measures)


results <- data.frame(
  EnvironmentalFactor = character(),
  DiversityMeasure = character(),
  Importance = numeric(),
  Correlation = numeric(),
  PValue = numeric(),
  stringsAsFactors = FALSE
)


for (factor in env_factors) {
  factor_dir <- file.path(output_dir, factor)
  dir.create(factor_dir, showWarnings = FALSE)
  
  for (diversity in diversity_measures) {
    if (diversity %in% colnames(alpha_diversity)) {
      
      if (is.factor(alpha_diversity[[diversity]]) || length(unique(alpha_diversity[[diversity]])) <= 5) {
        rf_model <- randomForest(
          x = alpha_diversity[, factor, drop = FALSE],
          y = alpha_diversity[[diversity]],
          importance = TRUE,
          proximity = TRUE,
          ntree = 1000,
          keep.forest = TRUE
        )
      } else {
        rf_model <- randomForest(
          x = alpha_diversity[, factor, drop = FALSE],
          y = alpha_diversity[[diversity]],
          importance = TRUE,
          proximity = TRUE,
          ntree = 1000,
          keep.forest = TRUE
        )
      }
      
      
      importance_vals <- importance(rf_model)
      importance_score <- importance_vals[1, 1]
      importance_score <- ifelse(importance_score < 0, 0, importance_score)
      
      
      correlation_value <- cor(
        alpha_diversity[[factor]], alpha_diversity[[diversity]],
        method = "spearman", use = "complete.obs"
      )
      cor_test <- cor.test(
        alpha_diversity[[factor]], alpha_diversity[[diversity]],
        method = "spearman"
      )
      p_value <- cor_test$p.value
      
      
      results <- rbind(results, data.frame(
        EnvironmentalFactor = factor,
        DiversityMeasure = diversity,
        Importance = importance_score,
        Correlation = correlation_value,
        PValue = p_value
      ))
      
      
      plot <- ggplot(alpha_diversity, aes_string(x = factor, y = diversity)) +
        geom_point() +
        geom_smooth(method = "lm") +
        theme_minimal(base_family = "SimHei") +
        labs(
          title = paste("Random Forest for", diversity, "vs", factor),
          x = factor,
          y = diversity
        )
      ggsave(
        file.path(factor_dir, paste0(diversity, "_vs_", factor, "_random_forest.png")),
        plot = plot, width = 8, height = 6, bg = "white"
      )
    }
  }
}


results$DiversityMeasure <- factor(results$DiversityMeasure, levels = diversity_measures)
results$EnvironmentalFactor <- factor(results$EnvironmentalFactor, levels = env_factors)


ggplot(results, aes(x = DiversityMeasure, y = EnvironmentalFactor)) +
  
  geom_point(
    data = subset(results, PValue > 0.05),
    aes(size = Importance),
    color = "gray",
    alpha = 0.3
  ) +
  
  geom_point(
    data = subset(results, PValue <= 0.05),
    aes(size = Importance, color = Correlation),
    alpha = 1
  ) +
  
  scale_color_gradient2(
    low = "blue", mid = "white", high = "#8B0000",
    midpoint = 0, name = "Correlation (%)"
  ) +
  
  scale_size(range = c(1, 10), name = "Variable Importance") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,         
      hjust = 1,          
      vjust = 1,          
      size = 14           
    ),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Diversity Index",
    y = "Environmental Factors",
    title = "Correlation and Importance of Environmental Factors"
  )

ggsave(
  file.path(output_dir, "fig3a.pdf"),
  width = 12, height = 6, bg = "white"
)
ggsave(
  file.path(output_dir, "fig3a.tiff"),
  width = 12, height = 6, bg = "white", dpi = 300
)

write.csv(results, file.path(output_dir, "fig3a.csv"), row.names = FALSE)