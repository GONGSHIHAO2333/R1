library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(patchwork)

input_file_hk_summer <- "gene_abundance2HKsummer.txt"
data_hk_summer <- read.table(input_file_hk_summer, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data_long_hk_summer <- data_hk_summer %>%
  pivot_longer(cols = -Functional_Zone, names_to = "Gene", values_to = "Value") %>%
  filter(!is.na(Value))

data_summary_hk_summer <- data_long_hk_summer %>%
  group_by(Functional_Zone, Gene) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

p_hk_summer <- ggplot(data_summary_hk_summer, aes(x = Functional_Zone, y = Mean_Value, group = Gene, color = Gene)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean_Value - SE, ymax = Mean_Value + SE), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),  
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
    axis.line = element_line(color = "black", size = 0.1),  
    axis.ticks = element_line(color = "black", size = 0.1),  
    strip.background = element_rect(color = "black", fill = "lightgray", size = 0),  
    strip.text = element_text(size = 8, face = "bold")  
  ) +
  labs(fill = "Relative abundance (Scaled TPM)", title = "HK Wet", y = "Relative abundance (Scaled TPM)") +
  theme(
    axis.title.y = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "italic")
  ) +
  facet_wrap(~Gene, nrow = 2, ncol = 4, scales = "free_y")

input_file_hk_winter <- "gene_abundance2HKwinter.txt"
data_hk_winter <- read.table(input_file_hk_winter, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data_long_hk_winter <- data_hk_winter %>%
  pivot_longer(cols = -Functional_Zone, names_to = "Gene", values_to = "Value") %>%
  filter(!is.na(Value))

data_summary_hk_winter <- data_long_hk_winter %>%
  group_by(Functional_Zone, Gene) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

p_hk_winter <- ggplot(data_summary_hk_winter, aes(x = Functional_Zone, y = Mean_Value, group = Gene, color = Gene)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean_Value - SE, ymax = Mean_Value + SE), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),  
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
    axis.line = element_line(color = "black", size = 0.1),  
    axis.ticks = element_line(color = "black", size = 0.2),  
    strip.background = element_rect(color = "black", fill = "lightgray", size = 0),  
    strip.text = element_text(size = 8, face = "bold")  
  ) +
  labs(fill = "Relative abundance (Scaled TPM)", title = "HK Dry", y = "Relative abundance (Scaled TPM)") +
  theme(
    axis.title.y = element_blank(),  
    legend.title = element_text(size = 8, face = "italic")
  ) +
  facet_wrap(~Gene, nrow = 2, ncol = 4, scales = "free_y")

input_file_qd_summer <- "gene_abundance3QDsummer.txt"
data_qd_summer <- read.table(input_file_qd_summer, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data_long_qd_summer <- data_qd_summer %>%
  pivot_longer(cols = -Functional_Zone, names_to = "Gene", values_to = "Value") %>%
  filter(!is.na(Value))

data_summary_qd_summer <- data_long_qd_summer %>%
  group_by(Functional_Zone, Gene) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

p_qd_summer <- ggplot(data_summary_qd_summer, aes(x = Functional_Zone, y = Mean_Value, group = Gene, color = Gene)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean_Value - SE, ymax = Mean_Value + SE), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),  
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
    axis.line = element_line(color = "black", size = 0.1),  
    axis.ticks = element_line(color = "black", size = 0.2),  
    strip.background = element_rect(color = "black", fill = "lightgray", size = 0),  
    strip.text = element_text(size = 8, face = "bold")  
  ) +
  labs(fill = "Relative abundance (Scaled TPM)", title = "QD Wet", y = "Relative abundance (Scaled TPM)") +
  theme(
    axis.title.y = element_text(size = 8, face = "bold"),  
    legend.title = element_text(size = 8, face = "italic")
  ) +
  facet_wrap(~Gene, nrow = 2, ncol = 4, scales = "free_y")

input_file_qd_winter <- "gene_abundance3QDwinter.txt"
data_qd_winter <- read.table(input_file_qd_winter, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

data_long_qd_winter <- data_qd_winter %>%
  pivot_longer(cols = -Functional_Zone, names_to = "Gene", values_to = "Value") %>%
  filter(!is.na(Value))

data_summary_qd_winter <- data_long_qd_winter %>%
  group_by(Functional_Zone, Gene) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n())
  )

p_qd_winter <- ggplot(data_summary_qd_winter, aes(x = Functional_Zone, y = Mean_Value, group = Gene, color = Gene)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean_Value - SE, ymax = Mean_Value + SE), width = 0.2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold"),  
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
    axis.line = element_line(color = "black", size = 0.1),  
    axis.ticks = element_line(color = "black", size = 0.2),  
    strip.background = element_rect(color = "black", fill = "lightgray", size = 0),  
    strip.text = element_text(size = 8, face = "bold")  
  ) +
  labs(fill = "Relative abundance (Scaled TPM)", title = "QD Dry", y = "Relative abundance (Scaled TPM)") +
  theme(
    axis.title.y = element_blank(),  
    legend.title = element_text(size = 8, face = "italic")
  ) +
  facet_wrap(~Gene, nrow = 2, ncol = 4, scales = "free_y")

combined_plot <- (p_hk_summer / p_qd_summer) | (p_hk_winter / p_qd_winter)

print(combined_plot)


output_file_pdf_combined <- "Figure2b.pdf"
ggsave(output_file_pdf_combined, plot = combined_plot, width = 6.4, height = 4.8, units = "in")
