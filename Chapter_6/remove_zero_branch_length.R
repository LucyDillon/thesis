# Load libraries
install.packages("ape")
library(ape)

# Load data
tree<- read.tree("core_gene_alignment.newick")

# Add small value to branch length
tree[["edge.length"]] <- tree[["edge.length"]]+1e-6 # could go even smaller, just demonstrating

write.tree(tree, file='tree_1e-06.nwk')
