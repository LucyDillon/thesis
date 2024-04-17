library(VennDiagram)

grid.newpage() 
draw.triple.venn(area1=1400, area2=1063, area3=570,  
                 n12=793, n23=484, n13=533, n123=475,  
                 category=c("NCBI","RGI(CARD)","ResFinder"), 
                 col="Black",fill=c("#3a86ff","#8338ec","#ff006e"))

