# Builds a list of cryptic exon splice sites
library(GenomicRanges)
library(dplyr)

# Script needs an annotation file, and a list of fidelity junction files
#anno <- read.table("~/Dropbox (Sydney Uni)/References/biomart_annotations_human.tab",
#                   header = TRUE, stringsAsFactors=FALSE)
#refrange <- GRanges(seqnames = anno$chromosome_name, ranges =  IRanges(anno$start_position, anno$end_position),
#                    strand = anno$strand, gene_name = anno$external_gene_name)
#fidelity.filelist <- c()

##############################################################################################
## Three functions:
## findAllJunctions() takes in the list of fidelity files and a reference range and
##                    outputs a data.frame of unique junctions
## findCrypticExons() takes in a list of anchored junctions and the gene names and
##                    output a data.frame of cryptic exons
## countCryptSkip() takes in the list of cryptic exon locations, a list of fidelity files
##                    and outputs a count data.frame of cryptic and skiptic junction counts
###############################################################################################


# Build a list of all junctions from the list of junction files
findAllJunctions <- function(fidelity.files, refrange) {

  for(eachfile in fidelity.files) {
    load(eachfile)

    # older versions have a weird format
    #if( names(annotated)[1] == "file"){
    #  annotated <- annotated$file
    #}

    # assumes that multiple samples are contained in the same rdata file
    # if only 1 sample present?
    #if(all(names(annotated) %in% c("ID", "counts", "all"))){
    #      annotated <- list(file = annotated)
    #}
    n.samples <- length(annotated)

    for(i in 1:n.samples) {

      ID = annotated[[i]]$ID
      message(ID)
      sample.data <- annotated[[i]]$all

      rbp_gr_sample <- GenomicRanges::GRanges(seqnames = sample.data$chr,
                               ranges = IRanges::IRanges(sample.data$start, sample.data$end),
                               strand = sample.data$strand,
                               count = sample.data$count, type = sample.data$type)
      hits <- GenomicRanges::findOverlaps(rbp_gr_sample, refrange, select = "first")
      rbp_gr_sample$gene_name = refrange$gene_name[hits]
      rbp_df_sample <- dplyr::tbl_df(data.frame(rbp_gr_sample) )

      # Remove counts in the intergenic regions
      rbp_df_sample <- rbp_df_sample %>% dplyr::filter(!is.na(gene_name))

      rbp_df_sample <- rbp_df_sample %>% dplyr::select(-count)
      if(!exists("all_junctions")) {
         all_junctions <- rbp_df_sample
      } else {
         all_junctions <- dplyr::bind_rows(all_junctions, rbp_df_sample)
      }

    } # loop over samples in each fidelity file
    all_junctions <- unique(all_junctions)

  } # Now loop over all provided files

  return(all_junctions)
} # end of function


# Find cryptic exons
# Algorithm goes through the 3' and 5' anchored juctions and find pairs where the
# unanchored ends are less than 500 base pairs apart
findCrypticExons <- function(anchored_junctions, genes) {
  for(gene in genes) {
    print(gene)
    a_junctions <- anchored_junctions %>% filter(gene_name == gene)
    anno_junctions <- annotated_junctions %>% filter(gene_name == gene)

    for(i in 1:length(anno_junctions)) {
       anno_start = as.integer(anno_junctions[i,"start"])
       anno_end = as.integer(anno_junctions[i,"end"])

       jmatches_start <- a_junctions %>% filter(start == anno_start & end < anno_end)
       jmatches_end <- a_junctions %>% filter(start > anno_start & end == anno_end)

       if(nrow(jmatches_start) == 0) { next }

       if(nrow(jmatches_end) == 0) { next }

       intron_width <- as.integer(anno_junctions[i, "width"])

       merge.matches <- full_join(jmatches_start, jmatches_end, by = "gene_name")
       merge.matches$exon_length <- intron_width - ( merge.matches$width.x + merge.matches$width.y )
       merge.matches <- merge.matches %>% filter(exon_length > 0 & exon_length < 500)

       if(nrow(merge.matches) == 0) { next }

       if(!exists("cryptic_exons")) {
           cryptic_exons <- merge.matches
       } else {
           cryptic_exons <- rbind(cryptic_exons, merge.matches)
       }
    } # Loop over annotated junctions
  } # Loop over genes

  return(cryptic_exons)
} # end of function

# Count cryptics and skiptics
countCryptSkip <- function(cryptic_exons, fidelity.files, refrange) {
   cryptic_exons_ids <- paste0(cryptic_exons$seqnames.x, ":", cryptic_exons$start.x, "-",
                            cryptic_exons$end.x, ":", cryptic_exons$strand.x)
   cryptic_exons_ids <- c(cryptic_exons_ids, paste0(cryptic_exons$seqnames.y, ":", cryptic_exons$start.y, "-",
                                                 cryptic_exons$end.y, ":", cryptic_exons$strand.y))
   cryptic_exons_ids <- unique(cryptic_exons_ids)


   counts.df <- data.frame(ID=character(),
                        skiptics=integer(),
                        cryptics=integer(),
                        totalc=integer(),
                        stringsAsFactors=FALSE)

   # Using the precomputed listed of cryptic exons, count the cryptics and skiptics
   for(eachfile in fidelity.files) {

      load(eachfile)


     if(all(names(annotated) %in% c("ID", "counts", "all"))){
       annotated <- list(file = annotated)
     }
      # The fidelity file Rdata should have an object called annotated
      n.samples <- length(annotated)
      for(i in 1:n.samples) {

         ID = annotated[[i]]$ID
         message(ID)
         sample.data <- annotated[[i]]$all

         rbp_gr_sample <- GRanges(seqnames = sample.data$chr,
                             ranges = IRanges(sample.data$start, sample.data$end),
                             strand = sample.data$strand,
                             count = sample.data$count, type = sample.data$type)
         hits <- findOverlaps(rbp_gr_sample, refrange, select = "first")
         rbp_gr_sample$gene_name = refrange$gene_name[hits]
         rbp_df_sample <- tbl_df(data.frame(rbp_gr_sample) )

        # Remove counts in the intergenic regions
        rbp_df_sample <- rbp_df_sample %>% dplyr::filter(!is.na(gene_name))
        rbp_df_sample <- rbp_df_sample %>% mutate(seq.id = paste0(seqnames, ":", start, "-", end, ":", strand))

        rbp_with_cryptics <- rbp_df_sample %>% filter(seq.id %in% cryptic_exons_ids)

        cryptic.counts <- sum(rbp_with_cryptics$count)
        skiptic.counts <- sum(rbp_df_sample %>% filter(type == "skiptic") %>% select(count))
        total.counts <- sum(rbp_df_sample %>%  select(count))

        new.row <- data.frame("ID" = ID, "skiptics" = skiptic.counts, "cryptics" = cryptic.counts,
                          "totalc" = total.counts)
        counts.df <- rbind(counts.df, new.row)
      } # Loop over the samples
   } # loop over each fidelity file

   return(counts.df)

} # End function

# Example code of running the two functions to get the cryptic exons
#all_junctions <- findAllJunctions(fidelity.filelist, refrange)
#anchored_junctions <- all_junctions %>% filter(type == "3'-anchored" | type == "5'-anchored")
#annotated_junctions <- all_junctions %>% filter(type == "annotated")
#skiptic_junctions <- all_junctions %>% filter(type == "skiptic")
#genes <- unique(anchored_junctions$gene_name)
#cryptic_exons <- findCrypticExons(anchored_junctions, genes)
#cryptSkipCounts <- countCryptSkip(cryptic_exons, fidelity.filelist, refrange)
