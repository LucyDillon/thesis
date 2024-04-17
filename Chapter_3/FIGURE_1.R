#Read in packages:
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(tibble)
#Read in data:
data<- read.csv("Accuracy_precision_recall_vals_r.csv")
#Subset data by method
OG_RGI <-subset(data, Method == "OG_RGI")
LR <- subset(data, Method == "LR")
RGI_S <- subset(data, Method == "RGI_specific")
RGI_A <- subset(data, Method == "RGI_all")
Egg <- subset(data, Method == "Eggnog")
# Make sure data is in the correct order
data$Method <- factor(data$Method,
                      levels = c('OG_RGI', 'LR','RGI_specific', 'RGI_all', 'Eggnog'),ordered = TRUE)
# Assign which level of significance each star gets 
map_signif_level <- function(p.value) {
  if (p.value < 0.0005) return("***")
  if (p.value < 0.005) return("**")
  if (p.value < 0.05) return("*")
  return("ns")
}
# Make a colour pallete
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#000000")
# Make a basic box plot with the points on and the outliers and assign to a value
A <- ggplot(data, aes(x=Method, y=Accuracy, fill=Method)) + 
  geom_boxplot(alpha=0.6, outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black") + 
  geom_point(aes(), position = position_jitter(width = 0.05))
# using ggplot_build(p) this will give info about the boxplot
gg<-ggplot_build(A)
#
gg$data[[1]]
#
xx<-gg$data[[1]][c("group","outliers")]
# Change the group values to the Method values
xx$group<-c("OG_RGI","LR", "RGI_specific", "RGI_all", "Eggnog")
#
data.new<-merge(data,xx,by.x="Method",by.y="group")
#
data.new$out<-apply(data.new,1,function(x) x$Accuracy %in% x$outliers)

# Precision

P <- ggplot(data, aes(x=Method, y=average_precision, fill=Method)) + 
  geom_boxplot(outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black", alpha=0.6) + 
  geom_point(aes(), position = position_jitter(w = 0.1, h = 0)) 

ggp <- ggplot_build(P)
#
ggp$data[[1]]
#
xxp<-ggp$data[[1]][c("group","outliers")]
# Change the group values to the Method values
xxp$group<-c("OG_RGI","LR", "RGI_specific", "RGI_all", "Eggnog")
#
data.newp <-merge(data,xxp,by.x="Method",by.y="group")
#
data.newp$out<-apply(data.newp,1,function(x) x$average_precision %in% x$outliers)

# Recall:

R <- ggplot(data, aes(x=Method, y=average_recall, fill=Method)) + 
  geom_boxplot(outlier.shape = 17, outlier.size = 2.5, outlier.colour= "black", alpha=0.6) + 
  geom_point(aes(), position = position_jitter(w = 0.1, h = 0)) 


ggr <- ggplot_build(R)
#
ggr$data[[1]]
#
xxr<-ggr$data[[1]][c("group","outliers")]
# Change the group values to the Method values
xxr$group<-c("OG_RGI","LR", "RGI_specific", "RGI_all", "Eggnog")
#
data.newr <-merge(data,xxr,by.x="Method",by.y="group")
#
data.newr$out<-apply(data.newr,1,function(x) x$average_recall %in% x$outliers)


# Now I will make the plots using a label rather than a significance bar:

Plot4 <- ggplot(data.new, aes(factor(Method), Accuracy, fill=Method)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15)) +
  geom_text(aes(x = 1, y = 107, label = "A")) +
  geom_text(aes(x = 2, y = 107, label = "B")) +
  geom_text(aes(x = 3, y = 107, label = "C")) +
  geom_text(aes(x = 4, y = 107, label = "D")) +
  geom_text(aes(x = 5, y = 107, label = "CD")) +
  scale_shape_manual(values=c(16,17),guide="none") + 
  scale_size_manual(values=c(1.5,2.5),guide="none") + 
  scale_x_discrete(name = "Method", labels=c('1', '2', '3', '4', '5')) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average accuracy (%)" )+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),) + 
  ggtitle("A.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c('Original RGI analysis', 'Logistic regression of RGI genes', 'J48 model using RGI specific genes', 'J48 model using RGI all genes', 'J48 model using Eggnog gene families'))


Plot5 <- ggplot(data.newp, aes(x=Method, y=average_precision, fill=Method)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_text(aes(x = 1, y = 107, label = "A")) +
  geom_text(aes(x = 2, y = 107, label = "A")) +
  geom_text(aes(x = 3, y = 107, label = "B")) +
  geom_text(aes(x = 4, y = 107, label = "C")) +
  geom_text(aes(x = 5, y = 107, label = "C")) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average precision (%)" )+
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("B.") +
  scale_fill_manual(values=cbbPalette, guide="none", labels=c('Original RGI analysis', 'Logistic regression of RGI genes', 'J48 model using RGI specific genes', 'J48 model using RGI all genes', 'J48 model using Eggnog gene families'))


Plot6 <- ggplot(data.newr, aes(x=Method, y=average_recall, fill=Method)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.6) + 
  geom_point(aes(shape=out,size=out), position = position_jitter(w=0.15))+
  geom_text(aes(x = 1, y = 107, label = "ABC")) +
  geom_text(aes(x = 2, y = 107, label = "A")) +
  geom_text(aes(x = 3, y = 107, label = "B")) +
  geom_text(aes(x = 4, y = 107, label = "C")) +
  geom_text(aes(x = 5, y = 107, label = "BC")) +
  scale_shape_manual(values=c(16,17),guide="none")+ #these numbers relate to the shape numbers (of the points)
  scale_size_manual(values=c(1.5,2.5),guide="none") + #these are the respective sizes
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits=c(0,109), name ="Average recall (%)" )+
  scale_x_discrete(labels=c('1', '2', '3', '4', '5')) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ggtitle("C.") +
  scale_fill_manual(values=cbbPalette, guide=FALSE, labels=c('1: Original RGI analysis', 
                                                              '2: Logistic regression of RGI genes', 
                                                              '3: J48 model using RGI specific genes', 
                                                              '4: J48 model using RGI all genes', 
                                                              '5: J48 model using Eggnog gene families'))


grid.arrange(Plot4,Plot5, Plot6, nrow=2, ncol=2)
                     
