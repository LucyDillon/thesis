# Read in libraries:
library(dplyr)
# Read in data:
data<- read.csv("Accuracy_precision_recall_vals_r.csv")

# Write function for wilcoxon signed rank test:
pairwise_wilcox_test <- function(data, column, groups){
  
  for (i in 1:length(groups)){
    sub_data <- subset(data, Method %in% groups[[i]])
    result <- wilcox.test(as.formula(paste(column, "~ Method")), data = sub_data,
                          exact = FALSE, paired = TRUE)
    print(paste("Comparison between", groups[[i]][[1]], "and", groups[[i]][[2]], ":", sep = " "))
    print(result)
    print("\n")
  }
}

# Perform stats test on accuracy values:
pairwise_wilcox_test(data, "Accuracy", list(c("OG_RGI", "LR"),
                                            c("OG_RGI", "RGI_specific"),
                                            c("OG_RGI", "RGI_all"),
                                            c("OG_RGI", "Eggnog"),
                                            c("LR", "RGI_specific"),
                                            c("LR", "RGI_all"),
                                            c("LR", "Eggnog"),
                                            c("RGI_specific", "RGI_all"),
                                            c("RGI_specific", "Eggnog"),
                                            c("RGI_all", "Eggnog")))

# Perform stats test on precision values:
pairwise_wilcox_test(data, "average_precision", list(c("OG_RGI", "LR"),
                                            c("OG_RGI", "RGI_specific"),
                                            c("OG_RGI", "RGI_all"),
                                            c("OG_RGI", "Eggnog"),
                                            c("LR", "RGI_specific"),
                                            c("LR", "RGI_all"),
                                            c("LR", "Eggnog"),
                                            c("RGI_specific", "RGI_all"),
                                            c("RGI_specific", "Eggnog"),
                                            c("RGI_all", "Eggnog")))

# Perform stats test on recall values:
pairwise_wilcox_test(data, "average_recall", list(c("OG_RGI", "LR"),
                                            c("OG_RGI", "RGI_specific"),
                                            c("OG_RGI", "RGI_all"),
                                            c("OG_RGI", "Eggnog"),
                                            c("LR", "RGI_specific"),
                                            c("LR", "RGI_all"),
                                            c("LR", "Eggnog"),
                                            c("RGI_specific", "RGI_all"),
                                            c("RGI_specific", "Eggnog"),
                                            c("RGI_all", "Eggnog")))


