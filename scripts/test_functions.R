library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(fidelio)

options(echo =TRUE)

# get intron database
hg38_introns <- "example/gencode_hg38_all_introns.bed.gz"
intron_db <- prepareIntrons(intronList=hg38_introns)

# add biomart - get gene IDs for intron db
biomart_annotation <- "example/biomart_annotations_human.tab"
biomart <- fread(biomart_annotation, data.table = FALSE)
intron_db$gene_name <- biomart$external_gene_name[ match( str_split_fixed(intron_db$EnsemblID, "\\.", 2)[,1], biomart$EnsemblID)  ]

files <- "example/test_junctions.txt"
files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

# add experiment column
files <- files %>%
  mutate(
    experiment_code = str_split_fixed(file_id, "_", 4)[,3],
    experiment_target = str_split_fixed(file_id, "_",4)[,1]
  )

rbp <- unique(files$experiment_target)

message(rbp)

# annotate junctions
annotated <- files[,1:2] %>%
  purrr::pmap(annotateJunctions, intron_db = intron_db)

# create proportions
annotated <- purrr::map(annotated, createSimpleProportions)

# create summary tables
proportion_table <-
  purrr::map( annotated, "proportions") %>%
  bind_rows()

p <- tidyr::gather(proportion_table, "key", "value", -ID) %>%
  tidyr::separate(key, into = c("method", "classification") ) %>%
  dplyr::filter( classification != "all") %>%
  ggplot( aes( x = classification, y = value, fill = ID) ) +
  geom_col(position = "dodge") +
  facet_wrap(~method, nrow = 2, scales = "free") +
  labs(title = "Test plot - U2AF2 or SNRNP70 knockdown") +
  theme_bw()

ggsave("example/example_plot.pdf")

readr::write_tsv( proportion_table, "example/example_table.txt" )

