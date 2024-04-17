library(ggplot2)
library(gridExtra)
library(tibble)

overall_cog <- read.csv("Eggnog_overall_cog_cat.csv")
Cellular_processes <- read.csv("Cellular_processes_table.csv")
Metabolism <- read.csv("Metabolism_table.csv")
Information_storage <- read.csv("Information_storage_table.csv")


overall_cog$genome_type <- factor(overall_cog$genome_type,
                     levels = c("core_genome", "softcore","shell", "cloud_genome"),ordered = TRUE)
Cellular_processes$genome_type <- factor(Cellular_processes$genome_type,
                                  levels = c("core_genome", "softcore","shell", "cloud_genome"),ordered = TRUE)
Metabolism$genome_type <- factor(Metabolism$genome_type,
                                  levels = c("core_genome", "softcore","shell", "cloud_genome"),ordered = TRUE)
Information_storage$genome_type <- factor(Information_storage$genome_type,
                                  levels = c("core_genome", "softcore","shell", "cloud_genome"),ordered = TRUE)




Palette1 <- c("#F11A73", "#47D8D4","#8F2CF2")
Palette2 <- c("#393860","#EEA9B0", "#CFF7C4", "#440C57", "#358A31","#41A293" ,"#6A2C30")
Palette3 <- c("#AB599C", "#1E88E5" ,"#DB1943" , "#AADBC2")
Palette4 <- c( "#19162F","#DEA5E7", "#303FD3",  "#DF07BB", "#004D40", "#F2EFC6",  "#ABB6A1", "#8CF379")
# Stacked + percent
P1 <- ggplot(overall_cog, aes(fill=Overall_cog_cat, y=X, x=genome_type)) + 
      geom_bar(position="fill", stat="identity")+
      labs(title = "A.  Overall COG category",
       x = "Pangenome section",
       y = "Percentage (%)") +
      theme_minimal()+
  theme(legend.position = "bottom")+
  scale_x_discrete(name = "Pangenome section", labels = c("core", "softcore", "shell", "cloud")) +
      scale_fill_manual(values=Palette1, guide="none", labels=c("Cellular processes & signaling", "Information storage & processing", "Metabolism"))

P2<- ggplot(Cellular_processes, aes(fill=COG, y=Percentage, x=genome_type)) + 
     geom_bar(position="fill", stat="identity")+
     labs(title = "B.  Cellular processes & signaling",
     x = "Pangenome section",
     y = "Percentage (%)") +
     theme_minimal()+
  theme(legend.position = "bottom")+
  scale_x_discrete(name = "Pangenome section", labels = c("core", "softcore", "shell", "cloud")) + # 
     scale_fill_manual(values=Palette2, guide="none")

P3 <- ggplot(Information_storage, aes(fill = COG, y = Percentage, x = genome_type)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "C. Information storage & processing",
       y = "Percentage(%)") +
  theme_minimal() +
  theme(legend.position = "bottom")+
  scale_x_discrete(name = "Pangenome section", labels = c("core", "softcore", "shell", "cloud")) + # , labels = c("core", "softcore", "shell", "cloud")
  scale_fill_manual(values = Palette3, guide="none")


P4<- ggplot(Metabolism, aes(fill=COG, y=Percentage, x=genome_type)) + 
     geom_bar(position="fill", stat="identity")+
     labs(title = "D.  Metabolism",
       x = "Pangenome section",
       y = "Percentage (%)") +
     theme_minimal()+
  theme(legend.position = "bottom")+
  scale_x_discrete(name = "Pangenome section", labels = c("core", "softcore", "shell", "cloud")) + # , labels = c("core", "softcore", "shell", "cloud")
     scale_fill_manual(values=Palette4, guide="none") #, guide="none"

grid.arrange(P1, P2, P3,P4, nrow=2, ncol=2)






