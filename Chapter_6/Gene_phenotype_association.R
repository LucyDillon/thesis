library(stats)
library(rstatix)
# Define a function for gene family analysis
perform_gene_family_analysis <- function(data) {
  chisq_result <- chisq.test(data)
  fisher_result <- fisher.test(data)
  pairwise_result <- pairwise_fisher_test(data)
  
  return(list(chisq_result = chisq_result, fisher_result = fisher_result, pairwise_result = pairwise_result))
}

# Read data from file
gene_data <- read.csv("merged_file.csv", header = TRUE)

# Get unique drug pairs
drug_pairs <- unique(gene_data$Drug_pair)

# Loop over drug pairs
results_list <- list()

for (drug_pair in drug_pairs) {
  # Subset data for the current drug pair
  drug_pair_data <- subset(gene_data, Drug_pair == drug_pair)
  
  # Loop over gene names within the drug pair
  gene_family_names <- unique(drug_pair_data$Gene_name)
  drug_pair_results <- list()
  
  for (gene_family_name in gene_family_names) {
    # Add a prefix "COG_" to each gene name
    prefixed_gene_family_name <- paste("COG_", gene_family_name, sep = "")
    
    # Subset data for the current gene family within the drug pair
    gene_family_data <- subset(drug_pair_data, Gene_name == gene_family_name)
    
    # Create a contingency table
    data_matrix <- as.matrix(gene_family_data[, c("SR", "SS", "RR", "RS")])
    rownames(data_matrix) <- gene_family_data$Contribution_type
    
    # Perform analysis
    results <- perform_gene_family_analysis(data_matrix)
    drug_pair_results[[prefixed_gene_family_name]] <- results
  }
  
  results_list[[drug_pair]] <- drug_pair_results
}

# Access results for a specific drug pair and gene family
print(results_list$ciprofloxacin_gentamicin$COG_COG2746)

sink("MDR_DT_STATS.txt")

for (drug_pair in names(results_list)) {
  drug_pair_results <- results_list[[drug_pair]]
  
  for (gene_family_name in names(drug_pair_results)) {
    cat("Results for", drug_pair, "-", gene_family_name, ":\n")
    print(drug_pair_results[[gene_family_name]])
    cat("\n")
  }
}

sink()

