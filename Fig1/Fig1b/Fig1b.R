setwd("*****************")

library(dplyr)
library(linkET)
library(ggplot2)


speciese <- read.csv("gene.csv", header = TRUE, row.names = 1)
varechem <- read.csv("antibiotic2.csv", header = TRUE, row.names = 1) 

mantel <- mantel_test(speciese, varechem,
                      spec_select = list(Nitrogen_Cycle_Gene = 1:59),
                      mantel_fun = "mantel",
                      env_ctrl = TRUE) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

correlation_plot <- qcorrplot(correlate(varechem, method = "spearman"), type = "upper", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = 0.1) +  
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +  
  labs(title = "Mantel test (Spearman)")

ggsave("Fig1b.pdf", plot = correlation_plot, width = 8, height = 6)

ggsave("Fig1b.tiff", plot = correlation_plot, width = 8, height = 6, dpi = 300)
write.table(mantel, "Fig1b.csv", row.names = FALSE, col.names = TRUE, sep = ",")
