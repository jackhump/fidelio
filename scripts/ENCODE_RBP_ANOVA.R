# redo ENCODE but this time bootstrap genes for confidence intervals
# and perform ANOVA to get P-values for comparing knockdown to controls


library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(fidelio)


#setwd("/Users/Jack/SAN/IoN_RNAseq/Kitty/Encode/RBPs/")

options(echo =TRUE)


hg38_introns <- "/SAN/vyplab/HuRNASeq/GENCODE/gencode.hg38.v29_all_introns.bed.gz"

intron_db <- prepareIntrons(intronList=hg38_introns)


files <- "/SAN/vyplab/HuRNASeq/GTEX/RBPs/all_ENCODE_junctions.txt"
files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

# add experiment column
files <- files %>%
  mutate(
    experiment_code = str_split_fixed(file_id, "_", 4)[,3],
    experiment_target = str_split_fixed(file_id, "_",4)[,1],
    replicate = str_split_fixed(file_id, "_",5)[,5] # watch out for this!
  )

# for each RBP
RBPs <- unique(files$experiment_target)

# for testing

RBPs <- RBPs[2]


fidelioPipeline <- function(files){
  result <-
    files[,1:2] %>%
    purrr::pmap(annotateJunctions, intron_db = intron_db) %>%
    fidelio::downSampleJunctions() %>%
    fidelio::filterJunctions(min_count = 10) %>%
    fidelio::bootstrapProportions(percent_genes = 0.5, times = 100)

  result$ANOVA <- fidelio::junctionANOVA(result)
}


for( rbp in RBPs){

  message(rbp)

  # get metadata for RBP
  experiments <- filter(files, experiment_target == rbp)

  # split by experiment



  # annotate junctions and filter
  annotated <-
    experiments[,1:2] %>%
    purrr::pmap(annotateJunctions, intron_db = intron_db)

  cleaned <-
    annotated %>%
    fidelio::filterJunctions(min_count = 10) %>%
    fidelio::downSampleJunctions() %>%
    fidelio::bootstrapProportions(percent_genes = 0.9, times = 100)

  junctionANOVA(cleaned)

}
#
#   # create proportions
#   annotated <- purrr::map(annotated, createSimpleProportions)
#
#   # look for knockdown of target RBP
#   annotated <- purrr::map(annotated,  knockdownStrength, gene = rbp )
#
#   # write out proportions table
#   proportions <- purrr::map_df( annotated, "proportions")
#   knockdown <- purrr::map_df( annotated, "countData")
#   # create out table
#   out <- proportions %>%
#     left_join( experiments, by = c("ID" = "file_id") ) %>%
#     left_join( knockdown )
#
#   # write out table
#   outFile <- paste0("RBPs/", rbp, "_fidelity.tab")
#   write.table(out, outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#
#   # save annotated as Rdata object
#   outFile <- gsub(".tab$", ".Rdata", outFile)
#   save(annotated, file = outFile)
#}

message( "all RBPs complete!")

