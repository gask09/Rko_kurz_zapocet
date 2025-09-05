if (!require("rentrez")) {
  install.packages("rentrez")
}

if (!require("glue")) {
  install.packages("glue")
}

if (!require("ape")) {
  install.packages("ape")
}

if (!require("stringr")) {
  install.packages("rentrez")
}

if (!require("pheatmap")) {
  install.packages("pheatmap")
}

library(rentrez)
library(glue)
library(ape)
library(stringr)
library(pheatmap)

#______________________Variables that can be modified__________________________

# !!! Set this variable according to your preferences !!!
OUTPUT_DIR_FOR_PLASTOMES="~/Desktop/annotated_plastomes"

if (!dir.exists(OUTPUT_DIR_FOR_PLASTOMES)) {
  dir.create(OUTPUT_DIR_FOR_PLASTOMES)
}
# !!! Set this variable according to your preferences. Make sure it is written correctly
# Otherwise the script will fail at some point !!!
TAXON_TO_SEARCH_FOR <- "Pleurothallidinae"

#______________________Function definitions____________________________________

# Queries the database for IDs of complete plastomes
# Override plastome sizes only if you know that they differ 
# from the predefined range i.e. you work with parasitic plants!
query_nucleotide_gb <- function(organism_name, 
                                lower_bound=135000, 
                                higher_bound=250000) {
  parsed_query <- glue("({organism_name}[Organism] OR {organism_name}[All Fields]) AND (ddbj_embl_genbank[filter] AND chloroplast[filter] AND ({lower_bound}[SLEN] : {higher_bound}[SLEN]))")
  queried_entries_count <- entrez_search(db="nucleotide", term=parsed_query, retmax=0)$count
  Sys.sleep(0.1)
  ids <- entrez_search(db="nucleotide", term=parsed_query, retmax=queried_entries_count)$ids
  return(ids)
} 

# Downloads a fasta file of plastome
get_fasta <- function(gb_id) {
  paphiopedilum_plastomes <- entrez_fetch(db="nucleotide", id=gb_id, rettype="fasta")
  filename_with_id <- glue("/tmp/plastome_seq_{gb_id}.fasta") 
  write(paphiopedilum_plastomes, file=filename_with_id)
}

# Before running chloe to annotate plastomes, it needs to be installed!
# See: https://github.com/ian-small/Chloe.jl
annotate_w_chloe <- function(fasta_path, 
                             output_name, 
                             pathway_to_chloe="~/chloe") {
  parsed_command <- glue("julia --project={pathway_to_chloe} {pathway_to_chloe}/chloe.jl annotate {fasta_path} -o {output_name} --sff --use-id")
  system(parsed_command)
}

parse_attribute <- function(gff_attribute) {
  split_attribute <- strsplit(x=gff_attribute, split=";")[[1]][3]
  gene_name <- strsplit(split_attribute, "=")[[1]][2]
  return(gene_name)
}


#________________________________Script itself___________________________________

ids <- query_nucleotide_gb(TAXON_TO_SEARCH_FOR)

for (id in ids) {
  get_fasta(id)
  annotate_w_chloe(fasta_path = glue("/tmp/plastome_seq_{id}.fasta"), output_name = OUTPUT_DIR_FOR_PLASTOMES)
}
# Just an empty list to hold vectors of gene names
list_of_individuals_with_all_genes <- list()

# All the magic happens here
for (filename in dir(OUTPUT_DIR_FOR_PLASTOMES, patter="*.gff")) {
  loaded_gff <- ape::read.gff(glue("{OUTPUT_DIR_FOR_PLASTOMES}/{filename}"))
  gff_genes_only <- loaded_gff[loaded_gff$type=="gene",]
  gff_attributes <- gff_genes_only$attributes
  gene_names <- c()
  for (attribute in gff_attributes) {
    gene_name <- parse_attribute(attribute)
    gene_names <- append(gene_names, gene_name)
  }
  ncbi_entry_id <- stringr::str_remove(filename ,".chloe.gff")
  list_of_individuals_with_all_genes[[ncbi_entry_id]] <- gene_names
}

# Unique set of all genes, definitely not the smartest way to do this, some genes might be missing if they are missing in all plastomes
all_genes <- unique(unlist(list_of_individuals_with_all_genes))

# Maps all present genes for each sample on a reference and recovers number of copies
# rps12 is transsplicing and thus occurs ~4x in genome
genes_presence_table <- sapply(list_of_individuals_with_all_genes, function(genes) {
  table(factor(genes, levels = all_genes))
})

# Appends row names and plots a heatmap
row.names(genes_presence_table) <- all_genes
rownames(genes_presence_table)
genes_presence_table_reordered <- genes_presence_table[order(rownames(genes_presence_table)),]


pheatmap(t(genes_presence_table_reordered), cluster_rows = FALSE, cluster_cols = FALSE)

