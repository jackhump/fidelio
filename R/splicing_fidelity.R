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
         gene_name = V4,
         EnsemblID = V5,
         strand = V6,
         transcript = V7,
         dot = V8,
         transcript_type = V9,
         transcript_annotation = V10) %>%
    dplyr::select( -dot ) %>%
    dplyr::mutate( end = end - 1 ) %>% # due to 0-base/1-base shenanigans
    dplyr::filter( !duplicated( paste(chr, start, end) ) ) # remove duplicate entries

  #anno <- data.table::fread(annotation, data.table = FALSE)
  #intron_db$gene_name <- anno$external_gene_name[ match( str_split_fixed(intron_db$EnsemblID, "\\.", 2)[,1], anno$EnsemblID)  ]


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
annotateJunctions <- function(file, intron_db, file_id, drop_chrM = FALSE, drop_sex = FALSE){
  # function to annotate a list of junctions according to a set of introns
  # presumably from GENCODE?

  if( !file.exists(file) ){
    message("file doesnt exist")
    return(NULL)
  }

  # check if gzipped
  if( grepl(".gz$", file) ){
    junctions <- data.table::fread(paste("zless", file), data.table=FALSE, sep = "\t") # in Rstudio you need, logical01=FALSE)
  }else{
    junctions <- data.table::fread(file, data.table=FALSE, sep = "\t") # in Rstudio you need, logical01=FALSE)
  }


  # work out whether a leafcutter junction file or a STAR junction file
  fileType <- "unknown"
  if(
    ncol(junctions) == 9 &
    all( sapply(junctions[,2:ncol(junctions)], is.integer) )
    ){
    fileType <- "STAR"
    message("file type is STAR")
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
  # junctions from leafcutter
  if(ncol(junctions) == 6){
  fileType <- "leafcutter"
  message("file type is leafcutter")
  sorted <- junctions %>%
    dplyr::arrange(V1, V2) %>%
    dplyr::mutate( V1 = paste0("chr", V1)) %>%
    dplyr::rename( chr = V1, start = V2, end = V3, dot = V4, count = V5, strand = V6) %>%
    dplyr::select( chr, start, end, count, strand )
  }
  # junctions from RegTools
  if(ncol(junctions) == 12){
    fileType <- "regtools"
    message("file type is regtools")
    sorted <- junctions %>%
      dplyr::arrange(V1,V2) %>%
      dplyr::rename(
        chr = V1,
        start = V2,
        end = V3,
        name = V4,
        count = V5,
        strand = V6,
        thickStart = V7,
        thickEnd = V8,
        itemRgb = V9,
        blockCount = V10,
        blockSizes = V11,
        blockStarts = V12
        ) %>%
      tidyr::separate(col = blockSizes, into = c("anchorStart", "anchorEnd"), sep = ",", remove = FALSE) %>%
      dplyr::mutate( anchorStart = as.numeric(anchorStart), anchorEnd = as.numeric(anchorEnd)) %>%
      dplyr::mutate( start = start + anchorStart, end = end - anchorEnd) %>%
      dplyr::select( chr, start, end, count, strand )
  }


  if( fileType == "unknown"){
    message("file not recognised")
    return(NULL)
  }

  if( drop_sex == TRUE){
    sorted <-
      filter(sorted, ! chr %in% c("chrX","chrY"))
  }

  if( drop_chrM == TRUE){
    sorted <-
      filter(sorted, chr != "chrM")
  }

  # find exact matches of both splice sites
  intersect_both <- sorted %>%
    dplyr::left_join(intron_db, by = c("chr","start", "end", "strand")) %>%
    dplyr::select( chr, start,end, count, strand, EnsemblID, gene_name, transcript  )

  # get just annotated junctions out
  annotated <- intersect_both %>%
    dplyr::filter( !is.na(transcript) ) %>%
    dplyr::mutate( type = "annotated") %>%
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript  )

  # by definition any junction that can't be found in the GENCODE intron table is novel
  novel_junctions <- dplyr::filter(intersect_both, is.na(transcript)) %>%
    dplyr::select( chr, start, end, strand, count)

  # annotate novel junctions to which genes they belong to
  intron_db_starts <-
    dplyr::select(intron_db, chr, start, strand, EnsemblID, gene_name) %>%
    dplyr::mutate(transcript = NA) %>%
    dplyr::distinct()

  intron_db_ends <-
    dplyr::select(intron_db, chr,end, strand, EnsemblID, gene_name) %>%
    dplyr::mutate(transcript = NA) %>%
    dplyr::distinct()
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
    dplyr::mutate( type = "skiptic") %>%
    left_join(intron_db_ends, by = c("chr","end", "strand")) %>%
    dplyr::distinct( chr, start,end,strand, .keep_all = TRUE) %>%
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript  )

  # SINGLY ANNOTATED JUNCTIONS

  # I can find left and right anchored junctions but I need to take strand into account
  # eg a left-anchored + junction is 5'-annotated and a right-anchored + junction is 3'-annotated
  # so I should categorise by this instead

  # annotate each junction with the gene that the anchored splice site belongs to

  # left-anchored
  anchored_start <- novel_junctions %>%
    dplyr::semi_join( intron_db, by = c("chr", "start", "strand") ) %>%
    dplyr::mutate( type = ifelse( strand == "+", "5'-anchored", "3'-anchored") ) %>%
    dplyr::left_join( intron_db_starts, by = c("chr","start", "strand") ) %>%
    dplyr::distinct( chr, start,end,strand, .keep_all = TRUE) %>% # remove any duplication weirdness
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript )

  # right anchored
  # + are 3'-anchored, - are 5'-anchored
  anchored_end <- novel_junctions %>%
    dplyr::semi_join( intron_db, by = c("chr", "end", "strand") ) %>%
    dplyr::mutate( type = ifelse( strand == "+", "3'-anchored", "5'-anchored") ) %>%
    dplyr::left_join( intron_db_ends, by = c("chr","end", "strand") ) %>%
    dplyr::distinct( chr, start,end,strand, .keep_all = TRUE) %>% # remove any duplication weirdness
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript )

  cryptic_anchored <- rbind( anchored_start, anchored_end) %>%
    dplyr::arrange( chr,start ) %>%
    dplyr::anti_join( skiptic, by = c("chr", "start", "end", "count", "strand") ) %>%
    dplyr::filter( !duplicated( paste(chr, start, end) ) )

  ## UNANCHORED JUNCTIONS

  cryptic_unanchored <- novel_junctions %>%
    dplyr::anti_join(skiptic, by = c("chr", "start", "end", "count", "strand")) %>%
    dplyr::anti_join(cryptic_anchored, by = c("chr", "start", "end", "count", "strand")) %>%
    dplyr::mutate( type = "cryptic_unanchored")

  # can I infer which genes an unanchored junction belongs to if it falls within the gene and is of the same strand?
  intron_db_gene_ranges <- intron_db %>%
    dplyr::group_by(EnsemblID) %>%
    dplyr::summarise( chr = unique(chr),
               start = min(start),
               end = max(end),
               gene_name = unique(gene_name),
               strand = unique(strand) ) %>%
    dplyr::arrange(chr, start)

  # full join and then filter on whether the start and end are within the range start and end
  cryptic_unanchored_annotated <-
    dplyr::full_join(cryptic_unanchored, intron_db_gene_ranges,
                     by = c("chr", "strand"),
                     suffix = c("", ".range")) %>%
    filter(start >= start.range & end <= end.range) %>%
    dplyr::distinct( chr, start,end,strand, .keep_all = TRUE) %>%
    mutate(transcript = NA) %>%
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript )

  cryptic_unanchored_novel <-
    dplyr::anti_join(cryptic_unanchored, cryptic_unanchored_annotated, by = c("chr", "start", "end", "count", "strand")) %>%
    mutate(EnsemblID = NA, gene_name = NA, transcript = NA) %>%
    dplyr::select( chr, start,end, strand, count, type, EnsemblID, gene_name, transcript )


  # bind all together
  all_junctions <- rbind( annotated, skiptic, cryptic_unanchored_novel, cryptic_unanchored_annotated, cryptic_anchored ) %>%
    dplyr::arrange(chr,start)

  # message(dim(all_junctions))
  # print(table(all_junctions$type))
  # print(head(intron_db))



  # create summary
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
createSimpleProportions <- function( data ){

  output <-
    map(data, ~{
  # remake counts from sample
    annotated_df <- .x
    annotated_df$counts <-
      group_by(annotated_df$all, type) %>%
      dplyr::summarise( n_unique = n(), sum_counts = sum(count)) %>%
      dplyr::mutate( prop_unique = n_unique / sum(n_unique),
                     prop_counts = sum_counts / sum(sum_counts) )

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
                                           propsum.skiptic = prop_sum_skiptic,
                                           stringsAsFactors = FALSE
    )
    return(annotated_df)
  })
  return(output)
}


#' Bootstrap estimate junction proportions
#' From list of junctions, remove a fixed proportion of a random set of genes and recalculate the proportions of each type of junction
#' Returns a dataframe summarising the proportions of each bootstrap
#' @param data
#' @param percent_genes
#' @param times
#'
#' @return
#' @export
#'
#' @examples
bootstrapProportions <- function(data, percent_genes = 0.9, times = 100){
  map( data, ~{
    res <- .x
    junctions <- res$all
    genes <- unique(junctions$gene_name)
    genes <- genes[!is.na(genes)]
    random_genes <-
      replicate(n = times,
                expr = sample(x = genes, size = ( length(genes) * percent_genes ) , replace = FALSE ),
                simplify = FALSE
      )
    res$bootstrap <-
      map_df( random_genes, ~{
        randomise <- res
        randomise$all <- dplyr::filter(randomise$all, gene_name %in% .x)
        #print(nrow(randomise$all))
        #return(randomise$all)
        random_props <- fidelio::createSimpleProportions(list(randomise))[[1]]

        return(random_props$proportions)
      })

    return(res)
  })
}


#' Filter sets of junctions
#'
#' @param data
#' @param min_count
#'
#' @return
#' @export
#'
#' @examples
filterJunctions <- function(data, min_count = 10){
  res <-
    map( data, ~{
      junctions <- .x$all
      filtered <- dplyr::filter(junctions, count >= min_count & !is.na(EnsemblID))
      .x$all <- filtered
      return(.x)
    })
  return(res)
}



#' Downsample set of junctions
#'
#' @param data
#' @param sampleSize
#'
#' @return
#' @export
#'
#' @examples
downSampleJunctions <- function(data){

  totals <-
    map_dbl(data, ~{
      sum(.x$all$count)
    })


  sampleSize <- min(totals)

  output <- map(data, ~{
    # downsample set of junctions to a number
    # weight sampling by counts
    junctions <- .x$all
    to_sample <-
      sample(1:nrow(junctions),
             prob = junctions$count,
             size = sampleSize,
             replace = TRUE )

    .x$all <-
      junctions[to_sample,] %>%
      group_by(chr,start,end,strand,type,EnsemblID, gene_name, transcript) %>%
      summarise( count = n() ) %>%
      select(chr,start,end,strand,count,type,EnsemblID, gene_name, transcript) %>%
      as_tibble()
    return(.x)
  })

	return(output)
}

