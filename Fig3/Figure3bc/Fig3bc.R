library(ggplot2)
library(dplyr)

alpha_diversity_path <- "alpha_diversity_results.csv"
env_factors_path <- "environment.csv"
output_dir <- "..........."

alpha_diversity <- read.csv(alpha_diversity_path, row.names = 1)
env_factors <- read.csv(env_factors_path, row.names = 1)

alpha_diversity <- cbind(env_factors, alpha_diversity)
alpha_diversity <- na.omit(alpha_diversity)

env_factors <- c("Total_Antibiotic_Stress")

diversity_measures <- unique(c("Denitrification", "DNR"))

results <- data.frame(Diversity = character(), Factor = character(), P_value = numeric(), R_squared = numeric(), AIC_Linear = numeric(), AIC_Polynomial = numeric(), Equation = character(), stringsAsFactors = FALSE)

for (diversity in diversity_measures) {
  if (diversity %in% colnames(alpha_diversity) &&
      is.numeric(alpha_diversity[[diversity]])) {
    
    for (factor in env_factors) {
      if (factor %in% colnames(alpha_diversity)) {
        
        lm_model   <- lm(as.formula(paste(diversity, "~", factor)), data = alpha_diversity)
        poly_model <- lm(as.formula(paste(diversity, "~ poly(", factor, ", 2)")), data = alpha_diversity)
        
        use_poly <- AIC(poly_model) < AIC(lm_model)
        best_model <- if (use_poly) poly_model else lm_model
        
        coefs <- coef(best_model)
        p_value   <- summary(best_model)$coefficients[2, 4]
        r_squared <- summary(best_model)$r.squared
        
        eqn <- if (use_poly) {
          sprintf("y = %.3fx² %+.3fx %+.3f", coefs[3], coefs[2], coefs[1])
        } else {
          sprintf("y = %.3fx %+.3f",  coefs[2], coefs[1])
        }
        
        x_pos <- min(alpha_diversity[[factor]]) + 0.1 * diff(range(alpha_diversity[[factor]]))  
        y_pos <- max(alpha_diversity[[diversity]]) - 0.1 * diff(range(alpha_diversity[[diversity]]))  
        
        plot <- ggplot(alpha_diversity, aes_string(x = factor, y = diversity)) +
          geom_point(color = "black", size = 2) +
          geom_smooth(method = "lm",
                      formula = y ~ poly(x, 2), 
                      se = TRUE, color = "blue", linewidth = 1.5) +
          annotate("text",
                   x = x_pos,  
                   y = y_pos,  
                   hjust = 0, vjust = 1, size = 5, 
                   label = paste0(
                     eqn,
                     "\n",
                     ifelse(p_value < 0.001, "P < 0.001",
                            ifelse(p_value < 0.01, "P < 0.01",
                                   ifelse(p_value < 0.05, "P < 0.05",
                                          paste0("P = ", format.pval(p_value, digits = 3))))),
                     "\nR² = ", round(r_squared, 3)
                   )) +
          theme_minimal(base_family = "serif") +
          theme(panel.grid  = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, size = 1),
                axis.text   = element_text(size = 16, color = "black"),
                axis.title  = element_text(size = 20, face = "bold"),
                plot.title  = element_blank()) +  
          labs(x = factor, y = paste(diversity, "Abundance (Scaled TPM)"))
        
        
        factor_dir <- file.path(output_dir, paste0(factor, "_analysis"))
        dir.create(factor_dir, showWarnings = FALSE, recursive = TRUE)
    
        
        ggsave(file.path(factor_dir, paste0(factor, "_vs_", diversity, "_best_model.png")),
               plot, width = 8, height = 6, bg = "white")
        ggsave(file.path(factor_dir, paste0(factor, "_vs_", diversity, "_best_model.pdf")),
               plot, width = 20, height = 15, units = "cm", bg = "white")
        
        
        results <- rbind(results,
                         data.frame(Diversity      = diversity,
                                    Factor         = factor,
                                    P_value        = p_value,
                                    R_squared      = r_squared,
                                    AIC_Linear     = AIC(lm_model),
                                    AIC_Polynomial = AIC(poly_model),
                                    Equation       = eqn,
                                    stringsAsFactors = FALSE))
        
        cat(sprintf("%s vs %s → 1P: %s | P = %.3g | R² = %.3f\n",
                    diversity, factor, eqn, p_value, r_squared))
      }
    }
  }
}


write.csv(results,
          file.path(output_dir, "all_factors_regression_results.csv"),
          row.names = FALSE)

cat("PNG/PDF.\n")
