library(dplyr)
library(stringr)
library(data.table)
# for a given junction file:
#setwd("/Users/Jack/SAN")
#setwd("/SAN/vyplab/HuRNASeq/GTEX/")
#file <- "HuRNASeq/GTEX/2016-11-11-fullrun1/SRR1413562.leafcutter_op"
#introns <- "gencode_hg19_all_introns.bed.gz"
#introns <- "/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafviz/annotation_codes/gencode_hg19/gencode_hg19_all_introns.bed.gz"

#

#' Read in intron data and create a tibble
#'
#' @param intronList
#'
#' @return intron_db tibble
#' @export
#'
#' @examples
prepareIntrons <- function(intronList){
  introns <- data.table::fread(paste( "zcat < ", intronList) , data.table=FALSE)
  intron_db <- introns %>%
    dplyr::rename(chr = V1,
         start = V2,
         end = V3,
         gene = V4,
         EnsemblID = V5,
         strand = V6,
         transcript = V7,
         dot = V8,
         transcript_type = V9,
         transcript_annotation = V10) %>%
    dplyr::select( -dot ) %>%
    dplyr::mutate( end = end - 1 ) %>% # due to 0-base/1-base shenanigans
    dplyr::filter( !duplicated( paste(chr, start, end) ) ) # remove duplicate entries
  return(intron_db)
}



#' Annotate a set of junctions
#'
#' @param file
#' @param intron_db
#' @param file_id
#'
#' @return
#' @export
#'
#' @examples
annotateJunctions <- function(file, intron_db, file_id){
  # function to annotate a list of junctions according to a set of introns
  # presumably from GENCODE

  if( !file.exists(file) ){
    message("file doesnt exist")
    return(NULL)
  }

  # check if gzipped
  if( grepl(".gz$", file) ){
    junctions <- data.table::fread(paste("zless", file), data.table=FALSE) # in Rstudio you need, logical01=FALSE)
  }else{
    junctions <- data.table::fread(file, data.table=FALSE) # in Rstudio you need, logical01=FALSE)
  }


  # work out whether a leafcutter junction file or a STAR junction file
  fileType <- "unknown"
  if(
    ncol(junctions) == 9 &
    all( sapply(junctions[,2:ncol(junctions)], is.integer) )
    ){
    fileType <- "STAR"

    sorted <- junctions %>%
      dplyr::rename( chr = V1,
              start = V2,
              end = V3,
              strand = V4,
              intron_motif = V5,
              annotated = V6,
              count = V7,
              multiCount = V8,
              maxOverhang = V9) %>%
      dplyr::arrange(chr, start) %>%
      dplyr::select( chr, start, end, count, strand ) %>%
      dplyr::filter(count > 0) %>% # filter out junctions with only multi-mapping reads
      dplyr::filter( strand != 0) %>%
      dplyr::mutate( start = start - 1,
              strand = ifelse( strand == 1, "+", "-" ))
  }
  if(ncol(junctions) == 6){
  fileType <- "leafcutter"
  sorted <- junctions %>%
    dplyr::arrange(V1, V2) %>%
    dplyr::mutate( V1 = paste0("chr", V1)) %>%
    dplyr::rename( chr = V1, start = V2, end = V3, dot = V4, count = V5, strand = V6) %>%
    dplyr::select( chr, start, end, count, strand )
  }

  if( fileType == "unknown"){
    message("file not recognised")
    return(NULL)
  }

  # find exact matches of both splice sites
  intersect_both <- sorted %>%
    dplyr::left_join(intron_db, by = c("chr","start", "end", "strand")) %>%
    dplyr::select( chr, start,end, strand, count, transcript )

  # get just annotated junctions out
  annotated <- intersect_both %>%
    dplyr::filter( !is.na(transcript) ) %>%
    dplyr::mutate( type = "annotated") %>%
    dplyr::select( chr, start, end, strand, count, type)

  # by definition any junction that can't be found in the GENCODE intron table is novel
  novel_junctions <- filter(intersect_both, is.na(transcript)) %>%
    dplyr::select( chr, start, end, strand, count)


  # split novel junctions into different types
  # skiptics - both ends are annotated but separately
  # anchored cryptics - only one end is annotated
  # cryptic_unanchored cryptics - neither end are annotated

  # semi join only keeps rows in X that match in Y
  # so only keep novel junctions where start and end match separately

  skiptic <- novel_junctions %>%
    dplyr::semi_join( intron_db, by = c("chr", "start", "strand") ) %>%
    dplyr::semi_join( intron_db, by = c( "chr", "end", "strand" ) ) %>%
    dplyr::arrange( chr,start )  %>%
    dplyr::mutate( type = "skiptic")

  # for singly annotated junctions
  # I can find left and right anchored junctions but I need to take strand into account
  # eg a left-anchored + junction is 5'-annotated and a right-anchored + junction is 3'-annotated
  # so I should categorise by this instead

  # left-anchored
  anchored_start <- novel_junctions %>%
    dplyr::semi_join( intron_db, by = c("chr", "start", "strand") ) %>%
    dplyr::mutate( type = ifelse( strand == "+", "5'-anchored", "3'-anchored") )

  # right anchored
  # + are 3'-anchored, - are 5'-anchored
  anchored_end <- novel_junctions %>%
    dplyr::semi_join( intron_db, by = c("chr", "end", "strand") ) %>%
    dplyr::mutate( type = ifelse( strand == "+", "3'-anchored", "5'-anchored") )

  cryptic_anchored <- rbind( anchored_start, anchored_end) %>%
    dplyr::arrange( chr,start ) %>%
    dplyr::anti_join( skiptic, by = c("chr", "start", "end", "count", "strand") ) %>%
    dplyr::filter( !duplicated( paste(chr, start, end) ) )

  cryptic_unanchored <- novel_junctions %>%
    dplyr::anti_join(skiptic, by = c("chr", "start", "end", "count", "strand")) %>%
    dplyr::anti_join(cryptic_anchored, by = c("chr", "start", "end", "count", "strand")) %>%
    dplyr::mutate( type = "cryptic_unanchored")

  # bind all together
  all_junctions <- rbind( annotated, skiptic, cryptic_unanchored, cryptic_anchored ) %>%
    dplyr::arrange(chr,start)

  # message(dim(all_junctions))
  # print(table(all_junctions$type))
  # print(head(intron_db))

  # add back in gene and transcript names for annotated junctions
  all_junctions <- all_junctions %>%
    dplyr::left_join(intron_db, by = c("chr","start", "end", "strand")) %>%
    dplyr::select( chr, start,end, count, strand, type, EnsemblID, gene_name, transcript )

  # create
  summary <- group_by(all_junctions, type) %>%
    dplyr::summarise( n_unique = n(), sum_counts = sum(count)) %>%
    dplyr::mutate( prop_unique = n_unique / sum(n_unique),
            prop_counts = sum_counts / sum(sum_counts) )

  #SRR_code <- gsub(".leafcutter_op", "", basename(file) )

  # test that everything worked
  if(
    ( nrow(cryptic_unanchored) + nrow(cryptic_anchored) + nrow(skiptic) == nrow(novel_junctions) ) &
    nrow(all_junctions) == nrow(intersect_both) ){
    return(
      list(
        ID = file_id,
        counts = summary,
        all = as.tbl(all_junctions) )
      )
  }else{
    message("Error! The sums don't add up")
    return("ERROR")
  }

}



#' Create simple proportion of annotated vs non annotated
#' Adds a new element to the list containing counts and relative proportions of each type of junction
#' @param annotated_df
#'
#' @return annotated_df
#' @export
#'
#' @examples
createSimpleProportions <- function( annotated_df ){
    sample <- annotated_df$counts
    ID <- annotated_df$ID
    # unique junctions
    prop_unique_3anchored <- sample$n_unique[1] / sum(sample$n_unique)
    prop_unique_5anchored <- sample$n_unique[2] / sum(sample$n_unique)
    prop_unique_unanchored <- sample$n_unique[4] / sum(sample$n_unique)
    prop_unique_skiptic <- sample$n_unique[5] / sum(sample$n_unique)
    prop_unique <- 1 - (sample$n_unique[3] / sum(sample$n_unique) )
    # sum total junctions
    prop_sum_3anchored <- sample$sum_counts[1] / sum(sample$sum_counts)
    prop_sum_5anchored <- sample$sum_counts[2] / sum(sample$sum_counts)
    prop_sum_unanchored <- sample$sum_counts[4] / sum(sample$sum_counts)
    prop_sum_skiptic <- sample$sum_counts[5] / sum(sample$sum_counts)
    prop_sum <- 1 - (sample$sum_counts[3] / sum(sample$sum_counts) )
    annotated_df$proportions <- data.frame(ID = ID,
                                           propunique.all = prop_unique,
                                           propsum.all = prop_sum,
                                           propunique.5anchored = prop_unique_5anchored,
                                           propunique.3anchored = prop_unique_3anchored,
                                           propunique.unanchored = prop_unique_unanchored,
                                           propunique.skiptic = prop_unique_skiptic,
                                           propsum.5anchored = prop_sum_5anchored,
                                           propsum.3anchored = prop_sum_3anchored,
                                           propsum.unanchored = prop_sum_unanchored,
                                           propsum.skiptic = prop_sum_skiptic)
  return(annotated_df)
}

#with(all_props, plot( prop_unique, prop_sum) )
#
# chr_sizes <- intron_db %>%
#         group_by(chr) %>%
#         summarise( end = max(end) )
#
#
# # what is the distribution of novel junctions per chromosome?
# chr_sizes <- intron_db %>%
#   group_by(chr) %>%
#   summarise( end = max(end) )
#
# novel_junctions_per_chr <- novel_junctions %>%
#         group_by( chr ) %>%
#         summarise( count = sum(count) ) %>%
#         left_join( chr_sizes, by = "chr" ) %>%
#         filter( !grepl( "GL|MT", chr) ) %>%
#         mutate(prop = (count /  end) * 1E6,
#                chr = factor(chr, levels = paste0("chr", c(1:22,"X","Y") ) ) ) %>%
#         arrange( chr )



