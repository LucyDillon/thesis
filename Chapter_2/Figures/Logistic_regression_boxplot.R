#Read in packages:
library(ggplot2)
library(ggsignif)
library(gridExtra)
#Read in data:
data<- read.csv("Logistic_regression_r_input.csv")
#Subset data by method
RGI <-subset(data, Tool == "RGI")
Resfinder<- subset(data, Tool == "ResFinder")
NCBI <- subset(data, Tool == "NCBI")

# Accuracy:

wilcox.test.results.RGRF <- wilcox.test(RGI$Accuracy, Resfinder$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.RGNC <- wilcox.test(RGI$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.RFNC <- wilcox.test(Resfinder$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#000000")

map_signif_level <- function(p.value) {
  if (p.value < 0.0005) return("***")
  if (p.value < 0.005) return("**")
  if (p.value < 0.05) return("*")
  return("ns")
}

p1 <- ggplot(data, aes(x=Tool, y=Accuracy, fill=Tool)) + 
  geom_boxplot(alpha=0.6, outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black") + 
  geom_point(aes(), position = position_jitter(width = 0.05)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("A.") +
  scale_fill_manual(values=cbbPalette, guide="none")


wilcox.test(RGI$Accuracy, Resfinder$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(Resfinder$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

# Precision

p2<- ggplot(data, aes(x=Tool, y=Average_precision, fill=Tool)) + 
  geom_boxplot(alpha=0.6, outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black") + 
  geom_point(aes(), position = position_jitter(width = 0.05))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) + 
  ggtitle("B.") +
  scale_fill_manual(values=cbbPalette, guide="none")


wilcox.test(RGI$Average_precision, Resfinder$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI$Average_precision, NCBI$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(Resfinder$Average_precision, NCBI$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95) 

# Recall:
p3 <- ggplot(data, aes(x=Tool, y=Average_recall, fill=Tool)) + 
  geom_boxplot(alpha=0.6, outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black") + 
  geom_point(aes(), position = position_jitter(width = 0.05))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) + 
  ggtitle("C.") +
  scale_fill_manual(values=cbbPalette, guide="none")


wilcox.test(RGI$Average_recall, Resfinder$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI$Average_recall, NCBI$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(Resfinder$Average_recall, NCBI$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

grid.arrange(p1, p2, p3, nrow=2, ncol=2)
