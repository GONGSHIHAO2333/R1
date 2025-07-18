library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


read_and_preprocess <- function(file_path) {
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data_long <- data %>% 
    pivot_longer(cols = -Functional_Zone, names_to = "Gene", values_to = "Value") %>% 
    filter(!is.na(Value))
  data_summary <- data_long %>% 
    group_by(Functional_Zone, Gene) %>% 
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SE = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  return(data_summary)
}

create_plot <- function(data, plot_title, show_x_axis = FALSE) {
  ggplot(data, aes(x = Functional_Zone, y = Mean, group = Gene, color = Gene)) +
    geom_line(linetype = "dotted", linewidth = 1) + 
    geom_point(size = 3) + 
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, linewidth = 0.8) + 
    theme_minimal(base_size = 12) + 
    theme(
      
      axis.text.x = if (show_x_axis) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black")
      } else {
        element_blank()
      },
      axis.title.x = if (show_x_axis) {
        element_text(size = 12, face = "bold", margin = margin(t = 10), color = "black")
      } else {
        element_blank()
      },
      
      axis.text.y = element_text(size = 10, color = "black"),  
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10), color = "black"),  
      axis.ticks.y = element_line(color = "black", linewidth = 0.6),  
      
      plot.title = element_text(size = 14, face = "bold", hjust = 0, color = "black"), 
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  
      plot.margin = margin(10, 10, 10, 10)  
    ) +
    labs(
      x = if (show_x_axis) "" else NULL,  
      y = "",  
      title = plot_title
    ) +
    facet_wrap(~Gene, nrow = 2, ncol = 4, scales = "free_y")
}


input_files <- list(
  hk_winter = "HKDry.txt",
  hk_summer = "HKWet.txt",
  qd_winter = "QDDry.txt",
  qd_summer = "QDWet.txt"
)

plots <- list()
plots$hk_winter <- create_plot(read_and_preprocess(input_files$hk_winter), "HK Dry", show_x_axis = FALSE)
plots$hk_summer <- create_plot(read_and_preprocess(input_files$hk_summer), "HK Wet", show_x_axis = FALSE)
plots$qd_winter <- create_plot(read_and_preprocess(input_files$qd_winter), "QD Dry", show_x_axis = FALSE)
plots$qd_summer <- create_plot(read_and_preprocess(input_files$qd_summer), "QD Wet", show_x_axis = TRUE)


final_plot <- plots$hk_winter + plots$hk_summer + plots$qd_winter + plots$qd_summer +
  plot_layout(ncol = 1, heights = c(1, 1, 1, 1)) &
  theme(legend.position = "none") &
  plot_annotation(tag_levels = NULL)


ggsave(
  filename = "C:/Users/shgong/Desktop/Figs3.pdf",
  plot = final_plot,
  width = 8, height = 16,
  units = "in",
  dpi = 600
)

