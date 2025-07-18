library(vegan)

denitrification <- read.csv("DNR.csv", row.names = 1)
environment2 <- read.csv("environment.csv", row.names = 1)

denitrification_Hill <- decostand(denitrification, method = "hellinger")

antibiotic_factors <- c("P1", "P2", "P3", "P4", "P5", 
                        "P6", "P7")
environmental_factors <- c("pH", "DO", "Temperature", "DOC", "NO3", "PO4", "Salinity")

antibiotic_data <- environment2[, colnames(environment2) %in% antibiotic_factors, drop = FALSE]
environmental_data <- environment2[, colnames(environment2) %in% environmental_factors, drop = FALSE]

if (ncol(antibiotic_data) == 0) {
  stop("none")
}

if (ncol(environmental_data) == 0) {
  stop("none")
}

var_varpart_combined <- varpart(denitrification_Hill, 
                                antibiotic_data, 
                                environmental_data)

tiff("anrvariance_partitioning_venn223.tiff", 
     width = 2500, height = 1500, res = 600)  # Increased resolution to 600 DPI
plot(var_varpart_combined, 
     bg = c("pink", "lightblue"), 
     main = "Venn Diagram of Variance Partitioning", 
     cex.main = 1.5, cex = 1.2)  # Adjust main title size

# Manually add enlarged axis labels
text(-0.5, 0.2, "Antibiotics", cex = 3, font = 2)  # Adjust position and font size
text(1, 0.2, "Non-Antibiotic Factors", cex = 3, font = 2)  # Adjust position and font size

dev.off()  # Close device

pdf("anrvariance_partitioning_venn223.pdf", 
    width = 10, height = 6)  # Reduce width and height for smaller output
plot(var_varpart_combined, 
     bg = c("pink", "lightblue"), 
     main = "Venn Diagram of Variance Partitioning", 
     cex.main = 2, cex = 1.8)  # Slightly smaller font size to reduce file size

# Manually add enlarged axis labels
text(-0.5, 0.2, "Antibiotics", cex = 1.5, font = 1)  # Adjust position and font size
text(1.5, 0.2, "Non-Antibiotic Factors", cex = 1.5, font = 1)  # Adjust position and font size

dev.off()  # Close device

