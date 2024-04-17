# Read in libraries
library(ggplot2)
library(gridExtra)
require(grid) 

# read in RGI genes across all genomes:
RGI_data <- read.csv("RGI_genes_by_drugclass.csv")

# calculate the drug class proportions:
RGI_data$Percentage <- (RGI_data$Genes / sum(RGI_data$Genes)) * 100

ggplot(RGI_data, aes(x = Drug_class, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Percentage of AMR Genes (%)") +
  xlab("Drug Class") +
  ggtitle("Percentage of AMR Genes per Drug Class ") +
  scale_fill_manual(values = c("skyblue", "plum1")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_y_continuous(breaks = c(0, 2, 4, 6 ), limits = c(0, 7))


# Now read in eggNOG data
eggnog_data <- read.csv("COG_cat_by_count_both_species.csv")

# calculate the eggnog proportions across all genomes:
eggnog_data$Percentage <- (eggnog_data$Count / sum(eggnog_data$Count)) * 100


eggnog_data$COG_cat <- factor(eggnog_data$COG_cat, levels = c("A", "B", "J", 
                                                              "K", "L", "D", 
                                                              "Y", "V", "T", "M",
                                                              "N","Z", "W", "U", 
                                                              "O", "C", "G","E", 
                                                              "F", "H","I", "P", 
                                                              "Q","R", "S"))


# Remove NA and poorly characterised from data:
reduced <- subset(eggnog_data, Overall_cog_cat != "NA")
reduced_data <- subset(eggnog_data, Overall_cog_cat != "POORLY CHARACTERIZED")

# Subset data for facet plot:
Cell <- subset(reduced_data, Overall_cog_cat == 'CELLULAR PROCESSES AND SIGNALING')
Metabolism <- subset(reduced_data, Overall_cog_cat == 'METABOLISM' )
information <- subset(reduced_data, Overall_cog_cat == 'INFORMATION STORAGE AND PROCESSING')


P1 <- ggplot(Cell, aes(x = COG_cat, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Percentage of total eggNOG gene families (%)") +
  scale_fill_manual(values = c("skyblue", "plum1"), guide="none") +
  theme(axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5)) +
  facet_wrap(~Overall_cog_cat, scales = "free_y")

P2<- ggplot(Metabolism, aes(x = COG_cat, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("skyblue", "plum1"), guide="none") +
  theme(axis.text.y = element_blank(),     # Remove y-axis labels
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5)) +
  facet_wrap(~Overall_cog_cat, scales = "free_y")

P3<- ggplot(information, aes(x = COG_cat, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("skyblue", "plum1")) +
  theme(
    axis.text.y = element_blank(),     # Remove y-axis labels
    axis.title.y = element_blank(), 
    axis.title.x = element_blank(),# Remove y-axis title
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8)  # Adjust the size as needed
  ) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 5)) +
  facet_wrap(~Overall_cog_cat, scales = "free_y")


grid.arrange(P1, P2, P3, nrow=1, ncol=3, widths = c(9/22, 8/22, 8/22),
             bottom=textGrob("eggNOG COG category", gp=gpar(fontsize=15)))
