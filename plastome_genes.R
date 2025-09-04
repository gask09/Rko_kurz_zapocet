library(rentrez)
library(glue)
library(ape)
library(stringr)
library(ggplot2)
install.packages("pheatmap")
library(pheatmap)
# !!! Set this variable according to your preferences !!!
OUTPUT_DIR_FOR_PLASTOMES="~/Desktop/paph_annotated_plastomes"

if (!dir.exists(OUTPUT_DIR_FOR_PLASTOMES)) {
  dir.create(OUTPUT_DIR_FOR_PLASTOMES)
}

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

get_gene_from_attribute <- function(gff_attribute) {
  split_attribute <- strsplit(x=gff_attribute, split=";")[[1]][3]
  gene_name <- strsplit(split_attribute, "=")[[1]][2]
  return(gene_name)
}


ids <- query_nucleotide_gb("Paphiopedilum")

for (id in ids) {
  get_fasta(id)
  annotate_w_chloe(fasta_path = glue("/tmp/plastome_seq_{id}.fasta"), output_name = OUTPUT_DIR_FOR_PLASTOMES)
}

list_of_individuals_with_all_genes <- list()

for (filename in dir(OUTPUT_DIR_FOR_PLASTOMES, patter="*.gff")) {
  loaded_gff <- ape::read.gff(glue("{OUTPUT_DIR_FOR_PLASTOMES}/{filename}"))
  gff_genes_only <- loaded_gff[loaded_gff$type=="gene",]
  gff_attributes <- gff_genes_only$attributes
  gene_names <- c()
  for (attribute in gff_attributes) {
    gene_name <- get_gene_from_attribute(attribute)
    gene_names <- append(gene_names, gene_name)
  }
  print(gene_names)
  ncbi_entry_id <- stringr::str_remove(filename ,".chloe.gff")
  print(ncbi_entry_id)
  list_of_individuals_with_all_genes[[ncbi_entry_id]] <- gene_names
}

all_genes <- unique(unlist(list_of_individuals_with_all_genes))

test <- sapply(list_of_individuals_with_all_genes, function(genes) {
  table(factor(genes, levels = all_genes))
})

row.names(test) <- all_genes

pheatmap(test,cluster_rows = FALSE, cluster_cols = FALSE, border_color = "black")
