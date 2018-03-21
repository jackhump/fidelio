#setwd("~/SAN/HuRNASeq/GTEX/Disease_Brain")
library(purrr)
library(stringr)
library(dplyr)
library(data.table)

source("../splicing_fidelity.R")

options(echo=TRUE)
# human 
introns <- "/SAN/vyplab/HuRNASeq/GTEX/gencode_hg38_all_introns.bed.gz"
biomart_annotation <- "/SAN/vyplab/HuRNASeq/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab" 


intron_db <- prepareIntrons(intronList=introns)

# add biomart - get gene IDs for intron db
biomart <- fread(biomart_annotation, data.table = FALSE)
intron_db$gene_name <- biomart$external_gene_name[ 
  match( 
    str_split_fixed(intron_db$EnsemblID, "\\.", 2)[,1],
    biomart$EnsemblID
    )  ]




sample_list <- c("ion_ftd_junctions.txt","prudencio_junctions.txt")
datasets <- c("IoN_FTD", "Prudencio_ALS")

for( i in 1:2){
  file <- sample_list[i]
  data <- datasets[i]
  outFile <- paste0(datasets[i], "_splicing_fidelity.tab")

  print(file)

  support <- read.table(file, header=TRUE, stringsAsFactors=FALSE)


    # annotate junctions
  annotated <- support %>% 
    select(file_id = sample, file = junctions) %>% 
    purrr::pmap(annotateJunctions, intron_db = intron_db)

  # create proportions
  annotated <- purrr::map(annotated, createSimpleProportions)


  # write out proportions table
  proportions <- purrr::map_df( annotated, "proportions")

  # get total counts and add to proportions table

  all_counts <- purrr::map(annotated, "counts")

  getTotalCounts <- function( countObject ){
    sum(pull(countObject, sum_counts))
  }

  # create out table
  out <- proportions
  # %>% 
    #left_join( samples, by = c("ID" = "sample") ) %>%
  #  mutate( total_junction_counts = purrr::map_int(all_counts, getTotalCounts) )

  # write out table
  write.table(out, outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # save annotated as Rdata object
  outRData <- gsub(".tab$", ".Rdata", outFile)
  save(annotated, file = outRData)
}

# in Rstudio

# what does the residual of the 5'/3' regression really mean?
# the values change depending on whether I keep all data together or split up by region, for example
# they are also normally distributed around 0 so what does this mean?

prudencio <- read.table("Prudencio_ALS_splicing_fidelity.tab", header=TRUE, stringsAsFactors=FALSE)
metadata <- str_split_fixed(prudencio$ID, "_", 3)
prudencio$condition <- metadata[,1]
prudencio$region <- metadata[,2]
prudencio$sample <- paste( metadata[,1], metadata[,3] ) 

ggplot( prudencio, 
        aes( x = prop_unique_5anchored, y = prop_unique_3anchored )) + 
  geom_point( aes(colour = condition)) + 
  geom_line( aes(y = fitted(mod1), linetype = region))

mod1 <- lm( prop_unique_3anchored ~ prop_unique_5anchored + region, data = prudencio )
mod2 <- glm( prop_unique_3anchored ~ prop_unique_5anchored + region, family = "poisson", data = prudencio)

prudencio$direction_bias <- residuals(lm( prop_unique_3anchored ~ prop_unique_5anchored, data = prudencio))

regress <- function(dat){
  mod <- lm(prop_unique_3anchored ~ prop_unique_5anchored, data = dat)
  return(residuals(mod))
}

per_tissue_direction_bias <- prudencio %>%
  split( .$region) %>%
  purrr::map( regress) %>%
  unlist()

prudencio <- prudencio %>%
  arrange(region) %>%
  mutate( ptdb = per_tissue_direction_bias)

ggplot(prudencio, aes(y = direction_bias, x = ptdb )) + 
  geom_point(aes(colour= condition))


ggplot( prudencio, aes( y = direction_bias, x = paste(condition, region))) + 
  #geom_boxplot(fill = NA) + 
  geom_point(aes(colour = condition)) +
  geom_line( aes(group = sample))

## FTD brain

ftd <- read.table("IoN_FTD_splicing_fidelity.tab", header=TRUE, stringsAsFactors = FALSE)
ftd$direction_bias <- residuals(mod <- lm(prop_unique_3anchored ~ prop_unique_5anchored, data = ftd))
ftd$condition <- gsub("_[A-Z0-9]*$", "", ftd$ID)
ggplot(ftd, aes( y = direction_bias, x = condition)) + geom_point()

ftd$fittedPois <- fitted(mod2)

mod2 <- glm( prop_unique_3anchored ~ prop_unique_5anchored, family = "poisson", data = ftd )

ggplot(ftd, aes(y = prop_unique_3anchored, x = prop_unique_5anchored, colour = condition)) + 
  geom_point() + 
  geom_abline(intercept=0,slope=1) +
  geom_line( aes( y = fitted ), colour = "red" ) +
  geom_line( aes( y = fittedPois), colour = "blue" )

ftd %>%
  select( ID, condition ) %>%
  mutate( gaus = residuals(mod), pois = residuals(mod2)) %>%
  tidyr::gather(model, value, -ID,-condition) %>%
  ggplot( aes( colour = condition, group = ID, shape = model)) + 
  geom_point( aes( y = value, x = paste(condition, model) ) ) 


