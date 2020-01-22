#
# Use biomaRt to download all the protein coding genes in a several QTL regions
#

# If not installed, install biomart by using:
#   install.packages("BiocManager")
#   BiocManager::install("biomaRt")


# example: genes <- getregion(bio.mart, 9, 107582143, 111667129) # biomart chr start end #äquivalent to the other example
# example: genes <- getregion(bio.mart, regions[1, "chr"], regions[1, "start"], region[1, "end"]) # biomart chr start end
# genes: get 120 genes

library(biomaRt)

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

getregion <- function(bio.mart, chr, startpos, endpos) {
  region <- paste0(chr, ":",startpos, ":", endpos)
  cat("function: ", " has region: ", region, "\n")
  res.biomart <- getBM(attributes = c("ensembl_gene_id",                                                    # Things that we want to get from biomart
                                      "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                       filters = c("chromosomal_region", "biotype"),                                        # Things that we will use to query biomart
                       values = list(region, "protein_coding"),                                             # The thing that we are querying
                       mart = bio.mart)

  cat("function: ", " has ", nrow(res.biomart), "\n")
  return (res.biomart) 
}

# test: getregion(bio.mart, 1, 10000000, 11000000), chromosome 1, start_position and end_position 
