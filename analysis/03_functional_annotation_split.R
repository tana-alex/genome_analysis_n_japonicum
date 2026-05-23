#!/usr/bin/env Rscript
# 03_functional_annotation.R


#==============================================================================
# PACKAGE LOADING
#==============================================================================

library("dplyr")

#==============================================================================
# SETUP
#==============================================================================

in_dir <- "data"
results_dir <- "results"

# Load upregulated genes from previous analysis
load(paste(results_dir, "deseq2_objects.RData", sep="/"))

#==============================================================================
# EXTRACT PROTEIN SEQUENCES
#==============================================================================

# Read FASTA file
fasta_lines <- readLines(paste(in_dir, "braker.aa", sep="/"))

#  upregulated genes list
upregulated_genes <- c("g17.t1", "g21.t1", "g25.t1", "g26.t1", "g30.t1", 
                       "g59.t1", "g64.t1", "g69.t1", "g79.t1", "g97.t1", 
                       "g98.t1", "g102.t1", "g110.t1", "g117.t1", "g123.t1", 
                       "g124.t1", "g128.t1", "g139.t1", "g145.t1", "g183.t1", 
                       "g187.t1", "g189.t1", "g193.t1", "g194.t1", "g207.t1", 
                       "g213.t1", "g217.t1", "g219.t1", "g231.t1", "g232.t1", 
                       "g236.t1", "g254.t1", "g259.t1", "g262.t1", "g268.t1", 
                       "g286.t1")

# Extract sequences
output_lines <- c()
keep_sequence <- FALSE

for(line in fasta_lines) {
  if(grepl("^>", line)) {
    gene_name <- gsub("^>", "", line)
    keep_sequence <- gene_name %in% upregulated_genes
  }
  
  if(keep_sequence) {
    output_lines <- c(output_lines, line)
  }
}

# Write to file
writeLines(output_lines, paste(results_dir, "heat_response_proteins.fa", sep="/"))

# Check results
cat("Number of sequences extracted:", length(grep("^>", output_lines)), "\n")

#==============================================================================
# LOAD ANNOTATION DATA
#==============================================================================

# Read EggNOG annotation file
eggnog_results <- read.delim(paste(in_dir, "annotation.emapper.annotations", sep="/"), 
                             header = FALSE,
                             stringsAsFactors = FALSE,
                             comment.char = "#",
                             sep = "\t")

# Add column names
colnames(eggnog_results) <- c(
  "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
  "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", 
  "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", 
  "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
)

# Read BLAST results
blast_results <- read.csv(paste(in_dir, "blast.tab", sep="/"), 
                          sep = "", header = FALSE)

colnames(blast_results) <- c(
  "query", "subject", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

#==============================================================================
# PROCESS ANNOTATIONS
#==============================================================================

# Get best BLAST hits
blast_best_hits <- blast_results %>%
  group_by(query) %>%
  arrange(evalue) %>%
  slice(1) %>%
  ungroup() %>%
  select(query, subject, pident, evalue, bitscore)

# Add confidence categories
blast_best_hits$confidence <- case_when(
  blast_best_hits$pident >= 80 ~ "High",
  blast_best_hits$pident >= 50 ~ "Good", 
  blast_best_hits$pident >= 25 ~ "Moderate",
  TRUE ~ "Low"
)

# Combine annotations
combined_annotation <- blast_best_hits %>%
  left_join(eggnog_results, by = "query", suffix = c("_blast", "_eggnog")) %>%
  select(query, 
         blast_subject = subject,
         blast_identity = pident,
         blast_evalue = evalue_blast,
         eggnog_name = Preferred_name,
         eggnog_description = Description,
         GO_terms = GOs,
         KEGG_pathways = KEGG_Pathway)

# Preview results
head(combined_annotation)

#==============================================================================
# TOP CANDIDATES ANALYSIS
#==============================================================================

#   top candidates analysis
top_candidates <- combined_annotation %>%
  filter(blast_identity > 70) %>%
  arrange(desc(blast_identity)) %>%
  head(8)

print("TOP CANDIDATES FOR DETAILED ANALYSIS:")
print(top_candidates)

#==============================================================================
# SAVE RESULTS
#==============================================================================

# Save combined results
write.csv(combined_annotation, paste(results_dir, "combined_functional_annotation.csv", sep="/"), 
          row.names = FALSE)

# Save top candidates
write.csv(top_candidates, paste(results_dir, "top_candidates.csv", sep="/"), 
          row.names = FALSE)

cat("Functional annotation complete! Results saved in:", results_dir, "\n")
