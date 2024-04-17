#Read in packages:
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(tibble)
#Read in data:
data<- read.csv("OG_analysis.csv")
#Subset data by Tool
ResFinder <-subset(data, Tool == "ResFinder")
RGI <- subset(data, Tool == "RGI")
NCBI <- subset(data, Tool == "NCBI")


# Make sure data is in the correct order
data$Tool <- factor(data$Tool,
                    levels = c('ResFinder', 'RGI','NCBI'),ordered = TRUE)
# Assign which level of significance each star gets 
map_signif_level <- function(p.value) {
  if (p.value < 0.0005) return("***")
  if (p.value < 0.005) return("**")
  if (p.value < 0.05) return("*")
  return("ns")
}
# Make a colour pallete
cbbPalette <- c("#ff006e", "#8338ec", "#3a86ff", "#F0E442", "#CC79A7", "#000000")
# Make a basic box plot with the points on and the outliers and assign to a value
A <- ggplot(data, aes(x=Tool, y=Accuracy, fill=Tool)) + 
  geom_boxplot(alpha=0.6, outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black") + 
  geom_point(aes(), position = position_jitter(width = 0.05))
# using ggplot_build(p) this will give info about the boxplot
gg<-ggplot_build(A)
#
gg$data[[1]]
#
xx<-gg$data[[1]][c("group","outliers")]
# Change the group values to the Tool values
xx$group<-c("ResFinder","RGI", "NCBI")
#
data.new<-merge(data,xx,by.x="Tool",by.y="group")
#
data.new$out<-apply(data.new,1,function(x) x$Accuracy %in% x$outliers)
Plot1 <- ggplot(data.new, aes(factor(Tool), Accuracy, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15)) +
  geom_signif(comparisons = list(c("ResFinder", "RGI")), y_position = 102, tip_length = 0.01, textsize = 2.5, vjust = 0.75, # These following lines add sig. bars and *s 
              test.stat = wilcox.test.results.OLa$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("ResFinder", "NCBI")),y_position = 107, tip_length = 0.01, textsize = 2.5, vjust = 0.75,# vjust changes the gap between text and the bar
              test.stat = wilcox.test.results.OSa$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("RGI", "NCBI")), y_position = 122, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.LSa$p.value, map_signif_level = map_signif_level) +
  scale_shape_manual(values=c(16,17),guide="none") + 
  scale_size_manual(values=c(1.5,2.5),guide="none") + 
  scale_x_discrete(labels=c('1', '2', '3')) +
  labs(y = "Average accuracy (%)", x = "Tool") +
  ylim(0, 137) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) + 
  ggtitle("A.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c("ResFinder","RGI", "NCBIpecific"))

# Average_precision
Average_precision <- subset(data, select = c(Antibiotic, Tool, Average_precision))

ResFinder_P <-subset(Average_precision, Tool == "ResFinder")
RGI_P <- subset(Average_precision, Tool == "RGI" )
NCBI_P <- subset(Average_precision, Tool == "NCBI")

wilcox.test.results.OLp <- wilcox.test(ResFinder_P$Average_precision, RGI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.OSp <- wilcox.test(ResFinder_P$Average_precision, NCBI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.LSp <- wilcox.test(RGI_P$Average_precision, NCBI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)


P <- ggplot(data, aes(x=Tool, y=Average_precision, fill=Tool)) + 
  geom_boxplot(outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black", alpha=0.6) + 
  geom_point(aes(), position = position_jitter(w = 0.1, h = 0)) 

ggp <- ggplot_build(P)
#
ggp$data[[1]]
#
xxp<-ggp$data[[1]][c("group","outliers")]
# Change the group values to the Tool values
xxp$group<-c("ResFinder","RGI", "NCBI")
#
data.newp <-merge(data,xxp,by.x="Tool",by.y="group")
#
data.newp$out<-apply(data.newp,1,function(x) x$Average_precision %in% x$outliers)

Plot2 <- ggplot(data.newp, aes(x=Tool, y=Average_precision, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_signif(comparisons = list(c("ResFinder", "NCBI")),y_position = 102, tip_length = 0.01, textsize = 2.5, vjust = 0.75,# vjust changes the gap between text and the bar
              test.stat = wilcox.test.results.OSp$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("ResFinder", "RGI")), y_position = 107, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.OAp$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("RGI", "NCBI")), y_position = 117, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.LSp$p.value, map_signif_level = map_signif_level) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
  labs(y = "Average Average_precision (%)", x = "Tool") +
  ylim(0, 137) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("B.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c("ResFinder","RGI", "NCBIpecific"))
# Average_recall:

Average_recall <- subset(data, select = c(Antibiotic, Tool, Average_recall))

ResFinder_R <-subset(Average_recall, Tool == "ResFinder")
RGI_R <- subset(Average_recall, Tool == "RGI" )
NCBI_R <- subset(Average_recall, Tool == "NCBI")



wilcox.test.results.ORGI <- wilcox.test(ResFinder_R$Average_recall, RGI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.OSr <- wilcox.test(ResFinder_R$Average_recall, NCBI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test.results.LSr <- wilcox.test(RGI_R$Average_recall, NCBI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

R <- ggplot(data, aes(x=Tool, y=Average_recall, fill=Tool)) + 
  geom_boxplot(outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black", alpha=0.6) + 
  geom_point(aes(), position = position_jitter(w = 0.1, h = 0)) 


ggr <- ggplot_build(R)
#
ggr$data[[1]]
#
xxr<-ggr$data[[1]][c("group","outliers")]
# Change the group values to the Tool values
xxr$group<-c("ResFinder","RGI", "NCBI")
#
data.newr <-merge(data,xxr,by.x="Tool",by.y="group")
#
data.newr$out<-apply(data.newr,1,function(x) x$Average_recall %in% x$outliers)

Plot3 <- ggplot(data.newr, aes(x=Tool, y=Average_recall, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_signif(comparisons = list(c("RGI", "NCBI")), y_position = 102, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.LSr$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("RGI", "RGI_all")), y_position = 107, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.LAr$p.value, map_signif_level = map_signif_level) +
  geom_signif(comparisons = list(c("RGI", "Eggnog")), y_position = 112, tip_length = 0.01, textsize = 2.5, vjust = 0.75,
              test.stat = wilcox.test.results.LEr$p.value, map_signif_level = map_signif_level) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
  labs(y = "Average Average_recall (%)", x = "Tool") +
  ylim(0, 137) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("C.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c('1: Original RGI analysis', 
                                                              '2: Logistic regression of RGI genes', 
                                                              '3: J48 model using RGI specific genes', 
                                                              '4: J48 model using RGI all genes', 
                                                              '5: J48 model using Eggnog gene families'))


grid.arrange(Plot1,Plot2, Plot3, nrow=2, ncol=2)


# Now I will make the plots using a label rather than a significance bar:

Plot4 <- ggplot(data.new, aes(factor(Tool), Accuracy, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15)) +
  geom_text(aes(x = 1, y = 107, label = "AB")) +
  geom_text(aes(x = 2, y = 107, label = "A")) +
  geom_text(aes(x = 3, y = 107, label = "B")) +
  scale_shape_manual(values=c(16,17),guide="none") + 
  scale_size_manual(values=c(1.5,2.5),guide="none") + 
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average accuracy (%)" )+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) + 
  ggtitle("A.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c("ResFinder","RGI", "NCBIpecific"))


Plot5 <- ggplot(data.newp, aes(x=Tool, y=Average_precision, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_text(aes(x = 1, y = 107, label = "A")) +
  geom_text(aes(x = 2, y = 107, label = "A")) +
  geom_text(aes(x = 3, y = 107, label = "A")) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average precision (%)" )+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("B.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c("ResFinder","RGI", "NCBIpecific"))


Plot6 <- ggplot(data.newr, aes(x=Tool, y=Average_recall, fill=Tool)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_text(aes(x = 1, y = 107, label = "AB")) +
  geom_text(aes(x = 2, y = 107, label = "A")) +
  geom_text(aes(x = 3, y = 107, label = "B")) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average recall (%)" )+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("C.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c('1: ResFinder', 
                                                              '2: RGI', 
                                                              '3: NCBI AMRFinder'))


grid.arrange(Plot4,Plot5, Plot6, nrow=2, ncol=2)

# Statistics:
wilcox.test(ResFinder$Accuracy, RGI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(ResFinder$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI$Accuracy, NCBI$Accuracy, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

wilcox.test(ResFinder_P$Average_precision, RGI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(ResFinder_P$Average_precision, NCBI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI_P$Average_precision, NCBI_P$Average_precision, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

wilcox.test(ResFinder_R$Average_recall, RGI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(ResFinder_R$Average_recall, NCBI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)
wilcox.test(RGI_R$Average_recall, NCBI_R$Average_recall, paired = TRUE, alternative = "two.sided", conf.level = 0.95)

library(stats)
Accuracy_p <- c(0.09524, 0.05939, 0.0002439)
Precision_p <- c(0.8438, 0.04637, 0.09375)
Recall_p <- c(0.06053, 0.1305, 0.001829)

p.adjust(Accuracy_p, method = "BH")

p.adjust(Precision_p, method = "BH")

p.adjust(Recall_p, method = "BH")
