
library(ggplot2)
library(dplyr)
library(gridExtra)

alpha_diversity_path <- "alpha_diversity_results.csv"
env_factors_path <- "environment.csv"
output_dir <- "fig4"

alpha_diversity <- read.csv(alpha_diversity_path, row.names = 1, na.strings = c("", "NA"))
env_factors <- read.csv(env_factors_path, row.names = 1, na.strings = c("", "NA"))

alpha_diversity <- cbind(env_factors, alpha_diversity)

col_type_status <- list()
for (col_name in colnames(alpha_diversity)) {
  tryCatch({
    alpha_diversity[[col_name]] <- as.numeric(alpha_diversity[[col_name]])
    col_type_status[[col_name]] <- "Succed"
  }, error = function(e) {
    col_type_status[[col_name]] <- paste("failed，season：", e$message)
  })
}


print(col_type_status)

env_factors <- c(
  "Total_Antibiotic_Stress"
)

diversity_measures <- unique(c( 
  "Denitrification_Bacteria", "Denitrification_Archaea", "Denitrification_Fungi", "DNR_Bacteria", "DNR_Archaea", "DNR_Fungi", "DNR_Archaea"
))


results <- data.frame(EnvironmentalFactor = character(),
                      DiversityMeasure = character(),
                      Spearman_Rho = numeric(),
                      P_Value = character(),
                      stringsAsFactors = FALSE)

min_p_value <- 0.001


plots_list <- list()


for (factor in env_factors) {
  dir.create(file.path(output_dir, factor), showWarnings = FALSE, recursive = TRUE)
  
  for (diversity in diversity_measures) {
    
    if (!is.numeric(alpha_diversity[[factor]]) || !is.numeric(alpha_diversity[[diversity]])) {
      warning(paste("1", factor, "2", diversity ))
      next
    }
    
    
    cor_test <- cor.test(
      alpha_diversity[[factor]], 
      alpha_diversity[[diversity]], 
      method = "spearman", 
      use = "complete.obs"
    )
    
    spearman_rho <- cor_test$estimate
    p_value <- cor_test$p.value
    
    if (p_value < min_p_value) {
      p_value_text <- "< 0.001"
    } else {
      p_value_text <- format.pval(p_value, digits = 3)
    }
    
    
    correlation_text <- paste0("ρ = ", round(spearman_rho, 3), ", P = ", p_value_text)
    
    
    plot <- ggplot(alpha_diversity, aes(x = .data[[factor]], y = .data[[diversity]])) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      annotate("text", 
               x = mean(range(alpha_diversity[[factor]], na.rm = TRUE)), 
               y = max(alpha_diversity[[diversity]], na.rm = TRUE) + 0.2 * diff(range(alpha_diversity[[diversity]], na.rm = TRUE)), # Y轴顶部上方偏移
               label = correlation_text,
               hjust = 0.5, vjust = 1, size = 4) +
      theme_bw() +
      theme(
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(family = "Times New Roman")  
      ) +
      labs(x = factor, y = diversity)
    
    
    plots_list[[length(plots_list) + 1]] <- plot
    
    ggsave(filename = file.path(output_dir, factor, paste0(diversity, "_vs_", factor, "_Spearman.png")),
           plot = plot, width = 3.2, height = 2.4, units = "in")
    
    ggsave(filename = file.path(output_dir, factor, paste0(diversity, "_vs_", factor, "_Spearman.pdf")),
           plot = plot, width = 3.2, height = 2.4, units = "in")
    
    results <- rbind(results, data.frame(EnvironmentalFactor = factor,
                                         DiversityMeasure = diversity,
                                         Spearman_Rho = spearman_rho,
                                         P_Value = p_value_text))
  }
}

write.csv(results, file.path(output_dir, "Fig4.csv"), row.names = FALSE)

combined_plot <- grid.arrange(grobs = plots_list, nrow = 2, ncol = 3)

ggsave(filename = file.path(output_dir, "Fig4.pdf"), plot = combined_plot, width = 10, height = 5, units = "in")

ggsave(filename = file.path(output_dir, "Fig4.tiff"), plot = combined_plot, width = 10, height = 5, units = "in", dpi = 300)

