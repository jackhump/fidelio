

# for use with recount rse_jx objects

annotate_recount_junctions <- function(rse_jx, intron_db, drop_sex = TRUE, drop_chrM = TRUE, filter_threshold = NA){
  junc_ranges <- rowRanges(rse_jx)
  junc_counts <- assays(rse_jx)$counts

  # if requested, filter by minimum count
  if( !is.na(filter_threshold)){
    to_filter <- apply(junc_counts >= filter_threshold , MARGIN = 1, FUN = any)
    junc_counts <- junc_counts[to_filter,]
    junc_ranges <- junc_ranges[to_filter,]
  }

  # get columns from the ranges info - use ID as count for now
  sorted <- junc_ranges %>%
    as.data.frame() %>%
    select( "chr" = seqnames, start, end, count = junction_id, strand) %>%
    mutate( chr = as.character(chr), strand = as.character(strand)) %>%
    arrange( chr, start)

  # 1-based / 0-based fuckery
  sorted <- mutate(sorted, start = start - 1)


  # do the annotating
  all_junctions <- annotate_junctions(sorted, intron_db, drop_sex = drop_sex, drop_chrM = drop_chrM)

  # now get counts for each sample
  junc_ranges$annotation <- all_junctions$type[ match(junc_ranges$junction_id, all_junctions$count)]
  # for each column in junc_ranges
  count_recount_sample <- function(sample_id, annotation){
    sample <- junc_counts[,sample_id]
    split(sample, f = annotation) %>%
      map( sum)
  }
  sample_counts <-
    map_df( colnames(junc_counts), count_recount_sample, annotation = junc_ranges$annotation)


  sample_props <- sweep(sample_counts, MARGIN = 1, STATS = rowSums(sample_counts), FUN = "/" ) %>%
    mutate( sample = colnames(junc_counts) ) %>%
    select( sample, everything() )

  sample_counts <- sample_counts %>%
    mutate( sample = colnames(junc_counts) ) %>%
    select( sample, everything() )

  return( list(junctions = all_junctions, counts = sample_counts, props = sample_props))

}

