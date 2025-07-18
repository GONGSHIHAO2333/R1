library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggsignif)
library(patchwork)
library(grid)

input_file <- "environment.csv"
data <- read.csv(input_file, stringsAsFactors = FALSE)

colnames(data) <- c("ANR", "Anammox", "Denitrification", "DNR", 
                    "Nitrification", "Nitrogen_fixation", "BAB", "Total_Pathway", 
                    "Location", "Season")
data <- data %>% mutate(across(c(ANR, Anammox, Denitrification, DNR, 
                                 Nitrification, Nitrogen_fixation, BAB, Total_Pathway), 
                               as.numeric))

variables <- c("ANR", "Anammox", "Denitrification", "DNR", 
               "Nitrification", "Nitrogen_fixation", "BAB", "Total_Pathway")

plots <- list()
test_results <- list()

get_significance <- function(p_value) {
  if (p_value < 0.001) return("***")
  else if (p_value < 0.01) return("**")
  else if (p_value < 0.05) return("*")
  else return("NS.")
}

data$Season <- factor(data$Season, levels = c("HKDry", "HKWet", "QDDry", "QDWet"))

season_to_num <- setNames(1:length(levels(data$Season)), levels(data$Season))

base_theme <- theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    
    axis.title.y = element_blank(), 
    plot.margin = margin(5, 5, 5, 5),  
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold", hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    panel.grid.major = element_line(size = 0.2),  
    panel.grid.minor = element_line(size = 0.1),  
    axis.ticks = element_line(size = 0.2),        
    axis.ticks.length = unit(2, "pt")             
  )


for (var in variables) {
  
  test_HK <- wilcox.test(data[[var]] ~ data$Season, 
                         subset = data$Season %in% c("HKDry", "HKWet"), 
                         exact = FALSE)
  test_QD <- wilcox.test(data[[var]] ~ data$Season, 
                         subset = data$Season %in% c("QDDry", "QDWet"), 
                         exact = FALSE)
  
  sig_HK <- get_significance(test_HK$p.value)
  sig_QD <- get_significance(test_QD$p.value)
  test_results[[var]] <- list("HKDry vs HKWet" = test_HK$p.value, 
                              "QDDry vs QDWet" = test_QD$p.value)
  
  
  min_val <- min(data[[var]], na.rm = TRUE)
  max_val <- max(data[[var]], na.rm = TRUE)
  y_range <- c(min_val - 0.4*(max_val - min_val), 
               max_val + 0.3*(max_val - min_val))
  y_pos <- max_val * 1.1  
  
  
  p <- ggplot(data, aes(x = as.numeric(Season), y = .data[[var]], fill = Season)) +  
    
    geom_boxplot(width = 0.5, size = 0.2, outlier.size = 0.5) + 
    scale_fill_brewer(palette = "RdYlBu") +  
    base_theme +
    labs(fill = var, title = var) +  
    
    scale_x_continuous(
      breaks = 1:length(levels(data$Season)),
      labels = levels(data$Season),
      limits = c(0.5, length(levels(data$Season)) + 0.5),  
      expand = expansion(add = 0)  
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.1)),
      breaks = pretty(data[[var]], n = 5)  
    ) +  
    coord_cartesian(expand = FALSE, ylim = y_range) +  
    geom_signif(
      comparisons = list(c(1, 2), c(3, 4)),
      map_signif_level = FALSE, 
      y_position = y_pos,  
      annotations = c(sig_HK, sig_QD),
      tip_length = 0.01,
      textsize = 2
    )
  
  
  if (var %in% c("ANR", "Nitrification")) {
    p <- p + theme(axis.title.y = element_text(size = 8, face = "bold", 
                                               angle = 90, 
                                               margin = margin(r = 10))) +
      labs(y = "Relative Abundance (Scaled TPM)")
  }
  
  plots[[var]] <- p
}

combined_plot <- wrap_plots(plots, ncol = 4) +
  plot_layout(guides = "collect") +  
  theme(plot.margin = margin(5, 5, 5, 5),  
        axis.ticks.length = unit(2, "pt"),      
        axis.text.y = element_text(margin = margin(r = 2)))  

print(combined_plot)  

output_file_pdf <- "C:/Users/shgong/Desktop/2024宏基因组分析结果/不同通路的总体丰度/输出/combined_boxplot_final.pdf"
ggsave(output_file_pdf, plot = combined_plot, width = 6.4, height = 4.8, units = "in")  

print(test_results)