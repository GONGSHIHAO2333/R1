
library("FactoMineR")
library("factoextra")
library("corrplot")


mydata <- read.csv("C:/Users/shgong/Desktop/ProClass.csv", header = TRUE, row.names = 1)


selected_columns <- c("Amoxicillin", "Ampicillin", "Azithromycin", "Cefalexin", "Cefotaxime", 
                      "Chlortetracycline", "Ciprofloxacin", "Clarithromycin", "Doxycycline", 
                      "Enrofloxacin", "Erythromycin_H2O", "Lincomycin", "Ofloxacin", "Oxytetracycline", 
                      "Roxythromycin", "Sulfadiazine", "Sulfamethazine", "Sulfamethoxazole", "Trimethoprim")  


data <- mydata[, selected_columns]


data <- scale(data)


res.pca <- PCA(data, graph = FALSE, ncp = 19)  

eig.val <- get_eigenvalue(res.pca)


print(eig.val)

pca_coords <- res.pca$ind$coord


dim(pca_coords) 


write.csv(pca_coords, "C:/Users/shgong/Desktop/PCA19.csv")


fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80))


var <- get_pca_var(res.pca)


corrplot(var$cos2, is.corr = FALSE)


