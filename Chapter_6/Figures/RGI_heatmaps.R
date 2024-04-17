# Read in libraries 
library(ggplot2)
library(viridis)
library(gridExtra)

# Read in RGI data across the pangenomes:
PA = read.csv("Pseudomonas_RGI_data.csv")
EC =read.csv("E_coli_RGI_data.csv")

P1<- ggplot(PA, aes(x = DrugClass, y = Best_Hit_ARO, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis(breaks = c(0, 15, 30, 45, 60, 75, 90)) +  # You can choose a different color scale if needed
  labs(title = expression(A. ~ italic("P. aeruginosa")),
       y = "Gene",
       x = "Drug Class") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
P2 <-ggplot(EC, aes(x = DrugClass, y = Best_Hit_ARO, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", breaks = c(0, 15, 30, 45, 60, 75, 90))+  # You can choose a different color scale if needed
  labs(title = expression(B. ~ italic("E. coli")),
       y = "Gene",
       x = "Drug Class") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

grid.arrange(P1, P2, nrow=1, ncol=2)



# Now lets read in the core genome data:
PA_core <- read.csv("Pseudomonas_core_rgi_genes.csv")
EC_core <- read.csv("E_coli_core_rgi_genes.csv")

P3 <- ggplot(PA_core, aes(x = Best_Hit_ARO, y = DrugClass, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis(breaks = c(95, 96, 97, 98, 99)) +  
  labs(title = expression(A. ~ italic("P. aeruginosa")),
       y = "Drug Class",
       x = "Gene") +
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = 45,hjust = 1 )) 

P4 <- ggplot(EC_core, aes(x =Best_Hit_ARO , y = DrugClass, fill = Percentage)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", breaks = c(95, 96, 97, 98, 99)) + 
  labs(title = expression(B. ~ italic("E. coli")),
       y = "Drug Class",
       x = "Gene") +
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

grid.arrange(P3, P4, nrow=2, ncol=1)

