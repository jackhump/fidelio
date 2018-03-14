# ENCODE RBPs
library(data.table)
library(dplyr)
library(purrr)
library(ggrepel)
#setwd("/Users/Jack/SAN/HuRNASeq/GTEX")
options(echo =TRUE)
source("splicing_fidelity.R")

hg38_introns <- "gencode_hg38_all_introns.bed.gz"

intron_db <- prepareIntrons(intronList=hg38_introns)

# add biomart - get gene IDs for intron db
biomart_annotation <- "../reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab" 
biomart <- fread(biomart_annotation, data.table = FALSE)
intron_db$gene_name <- biomart$external_gene_name[ match( str_split_fixed(intron_db$EnsemblID, "\\.", 2)[,1], biomart$EnsemblID)  ]

# are all the gene names in ENCODE recognised in Ensembl biomart or GENCODE?

#RBP_genes <- readLines("RBPs/all_RBP_genes.txt")
#RBP_genes[!RBP_genes %in% intron_db$gene_name]

# intronless genes - can't estimate differential expression from junction coverage, ha!
#"CSTF2T"  
#"DDX28"
#"HNRNPA0" 
#"PCBP1"
#"UTP3"

files <- "RBPs/all_ENCODE_junctions.txt"
files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

# add experiment column
files <- files %>% 
  mutate(
    experiment_code = str_split_fixed(file_id, "_", 4)[,3],
    experiment_target = str_split_fixed(file_id, "_",4)[,1] 
    )

# for each RBP 
RBPs <- unique(files$experiment_target)

# for finding strength of RBP knockdown in each dataset
knockdownStrength <- function(data, gene){
  ID <- data$ID
  libSize <-  sum(data$all$count)
  geneCounts <- filter(data$all, gene_name == gene) %>%
                .$count %>%
                sum()
  geneCPM <- (geneCounts / libSize ) * 1E6
  
  data$countData <- data.frame(ID, libSize, target = gene, geneCounts, geneCPM, stringsAsFactors = FALSE )
  return(data)
}

for( rbp in RBPs){

  message(rbp)
  
  # get metadata for RBP
  experiments <- filter(files, experiment_target == rbp)
  
  # annotate junctions
  annotated <- experiments[,1:2] %>%
    purrr::pmap(annotateJunctions, intron_db = intron_db)
  
  # create proportions
  annotated <- purrr::map(annotated, createSimpleProportions)
  
  # look for knockdown of target RBP
  annotated <- purrr::map(annotated,  knockdownStrength, gene = rbp )
  
  # write out proportions table
  proportions <- purrr::map_df( annotated, "proportions")
  knockdown <- purrr::map_df( annotated, "countData")
  # create out table
  out <- proportions %>% 
    left_join( experiments, by = c("ID" = "file_id") ) %>%
    left_join( knockdown )
  
  # write out table
  outFile <- paste0("RBPs/", rbp, "_fidelity.tab")
  write.table(out, outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # save annotated as Rdata object
  outFile <- gsub(".tab$", ".Rdata", outFile)
  save(annotated, file = outFile)
}

message( "all RBPs complete!")

exit()

# for transfer to Rstudio

#proportions2 <- read.table("RBPs/ALS_RBPs_fidelity.tab",header=TRUE, stringsAsFactors=FALSE) 

#proportions <- read.table("RBPs/all_RBPs_fidelity.tab",header=TRUE, stringsAsFactors=FALSE) 

# summarised counts are now in separate files for each RBP
# have to read each in an add to a giant list for plotting/exploration

readRBP <- function(rbp){
  message(rbp)
  fread(paste0("RBPs/", rbp, "_fidelity.tab"), data.table=FALSE, logical01 = FALSE)
}
  
proportions <- purrr::map_df( RBPs, readRBP)


# first test plot
# ggplot( proportions,
#         aes( x = prop_unique, y = prop_sum, colour = condition )) +
#   geom_text(aes(label = condition)) +
#   facet_wrap(~cell_type)

proportions <- mutate(proportions, experiment = str_split_fixed(ID, "_", 4)[,3])
  
# match in the specific control - is there a batch effect in the relative directions?
meta <- read.table("RBPs/all_ENCODE_metadata.txt",header=TRUE, sep = "\t", stringsAsFactors = FALSE)

proportions <- meta %>%
  select( Accession, Controls) %>%
  left_join(x = proportions, y = ., by = c( "experiment" = "Accession"))

save( proportions, meta, file = "RBPs/all_ENCODE_proportions.Rdata" )

#######################
# in Rstudio or whatever:

# look at spread of controls and variance between replicates
load("RBPs/all_ENCODE_proportions.Rdata")

proportions %>%
  filter( condition == "control") %>%
  #filter( cell_type == "K562") %>%
  mutate(ID = gsub("_[0-9]$", "", ID) ) %>%
  mutate( dataset = gsub("_[^_]+$", "", ID) ) %>%
  ggplot( aes(x = prop_sum_3anchored, y = prop_sum_5anchored, colour = condition, label = condition, group = ID)) + 
  geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  geom_text() + 
  facet_wrap( ~cell_type, scales = "free") +
  guides(colour=FALSE)


# give the mean value for the two replicates
proportions_clean <- proportions %>%
  mutate(ID = gsub("_[0-9]$", "", ID) ) %>%
  group_by( ID ) %>%
  dplyr::summarise(
             prop_unique = mean(prop_unique), 
             prop_sum = mean(prop_sum),
             prop_unique_5anchored = mean(prop_unique_5anchored),
             prop_unique_3anchored = mean(prop_unique_3anchored),
             prop_unique_unanchored = mean(prop_unique_unanchored),
             prop_unique_skiptic = mean(prop_unique_skiptic),
             prop_sum_5anchored = mean(prop_sum_5anchored),
             prop_sum_3anchored = mean(prop_sum_3anchored),
             prop_sum_unanchored = mean(prop_sum_unanchored),
             prop_sum_skiptic = mean(prop_sum_skiptic),
             condition = first(condition),
             cell_type = first(cell_type)) %>%
  mutate( dataset = gsub("_[^_]+$", "", ID) ) %>%
  arrange(desc(prop_sum), desc(prop_unique))



proportions_clean %>%
  #filter( condition == "control" ) %>%
  ggplot(
          aes( x = prop_unique, y = prop_sum, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 4) +
  theme_bw() +
  ylab("Proportion of total novel junctions") +
  xlab("Proportion of unique novel junctions") +
  labs(title = "Each knockdown linked to its control")
ggsave("RBPs/all_RBPs_scatter_fidelity.pdf", width = 15, height = 10)

# compare 5'anchored to 3'anchored splice sites
proportions_clean %>%
  #filter( condition == "control" ) %>%
  ggplot(
    aes( x = prop_sum_5anchored, y = prop_sum_3anchored, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 4) +
  theme_bw() +
 # ylab("Proportion of total novel junctions") +
 # xlab("Proportion of unique novel junctions") +
  labs(title = "Each knockdown linked to its control") +
  geom_abline(slope = 1, intercept =0, linetype=3)




# explore the plot with plotly!
library(plotly)
ggplotly()


# can I link each knockdown to its control and plot the relative change instead?

relative <- proportions_clean %>%
                 mutate(experiment = str_split_fixed(ID, "_", 4)[,3]) %>%
                 split( .$condition == "control" )    

relative <- left_join(relative[[1]],
                      relative[[2]], 
                      by = c("experiment","cell_type", "dataset"),
                      suffix = c("_kd", "_ctl")) %>%
            mutate( 
              rel_prop_unique = prop_unique_kd - prop_unique_ctl,
              rel_prop_sum = prop_sum_kd - prop_sum_ctl,
              rel_sum_5anchored = prop_sum_5anchored_kd - prop_sum_5anchored_ctl,
              rel_sum_3anchored = prop_sum_3anchored_kd - prop_sum_3anchored_ctl,
              rel_sum_unanchored = prop_sum_unanchored_kd - prop_sum_unanchored_ctl,
              rel_sum_skiptic = prop_sum_skiptic_kd - prop_sum_skiptic_ctl,
              rel_unique_5anchored = prop_unique_5anchored_kd - prop_unique_5anchored_ctl,
              rel_unique_3anchored = prop_unique_3anchored_kd - prop_unique_3anchored_ctl,
              rel_unique_unanchored = prop_unique_unanchored_kd - prop_unique_unanchored_ctl,
              rel_unique_skiptic = prop_unique_skiptic_kd - prop_unique_skiptic_ctl,
              condition = condition_kd)

relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_sum_5anchored, y = rel_sum_3anchored )) +
  geom_text( aes(label = condition_kd), size = 4) +
  #geom_point() +
  guides(colour = FALSE) +
  facet_wrap(~cell_type, scales = "free") + 
  theme_bw() +
  #ylab("Proportion of total novel junctions\nrelative to control") +
  #xlab("Proportion of unique novel junctions\nrelative to control") +
  labs(title = "All ENCODE RBPs") +
  geom_abline(slope = 1,intercept=0)

# create ratio of 5'-anchored to 3'-anchored
#relative$ratio <- relative$rel_sum_5anchored /  relative$rel_sum_3anchored
# doesn't work due to skew of 5'/3' - instead fit linear model and take residuals
mod <- lm( rel_sum_3anchored ~ rel_sum_5anchored, data = relative )
relative$residuals <- residuals(mod)
library(ggrepel)
# add jittered x
relative$jitter <- rnorm(nrow(relative), mean = 0, sd = 1)

relative %>% 
  ggplot(
    aes(y = residuals,  x = jitter, label = condition_kd)
  ) +
  geom_point( ) +
  geom_text_repel( data = tail( arrange(relative, abs(residuals)), 30 ) ) +
  facet_wrap(~cell_type, scales = "free") +
  theme_bw()

# kitty has counts of junctions by gene
getCounts <- function(rbp){
  message(rbp)
  load( paste0("../..//IoN_RNAseq/Kitty/Encode/RBPs/", rbp, "_counts.Rdata"))
  purrr::map_df(counts, "counts")
}

genes <- purrr::map_df( RBPs, getCounts )

"/SAN/vyplab/IoN_RNAseq/Kitty/Encode/RBPs/"




relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_sum_unanchored, y = rel_prop_sum )) +
  geom_text( aes(label = condition_kd), size = 3) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type, scales = "free") + 
  theme_bw() +
  #ylab("Proportion of total novel junctions\nrelative to control") +
  #xlab("Proportion of unique novel junctions\nrelative to control") +
  labs(title = "All ENCODE RBPs")



ggsave("RBPs/all_RBPs_scatter_fidelity_relative.pdf")

# for one-slide meeting

relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_anchored, y = rel_sum_anchored )) +
  geom_text( aes(label = condition_kd), size = 2.5) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type, scales = "free") + 
  theme_bw() +
  ylab("Proportion of total singly novel junctions\nrelative to control") +
  xlab("Proportion of unique singly novel junctions\nrelative to control") +
  labs(title = "All ENCODE RBPs")
ggsave("~/Documents/Misc/encode_rbps_scatter.pdf", width = 12, height = 8)

ggplotly()

zscore <- function(x){
  # return vector
  z <- (x - mean(x) ) / sd(x)
}

relative$prop_sum_Z <- zscore(relative$rel_prop_sum)
relative$prop_unique_Z <- zscore(relative$rel_prop_unique)

select(relative, condition, cell_type, prop_sum_Z, prop_unique_Z) %>%
  mutate(sum_Z = prop_sum_Z + prop_unique_Z) %>%
  arrange( desc(sum_Z) ) %>%
  ggplot( aes(y = sum_Z, x = sum_Z, colour = cell_type, label = condition )) + geom_text()


relative_controls <- meta %>%
          select( Accession, Controls) %>%
          left_join(x = relative, y = ., by = c( "experiment" = "Accession"))



# colour by control
# relative_controls %>%
#   #filter( condition != "HNRNPC" ) %>%
#   ggplot(
#     aes( x = rel_prop_unique, y = rel_prop_sum, colour = Controls )) +
#   #geom_text(aes(label = condition)) +
#   #geom_line( aes(group = dataset), colour = "black", linetype = 2) +
#   geom_point( aes(alpha = 0.2) ) +
#   guides(colour = FALSE) +
#   facet_wrap(~cell_type) + 
#   theme_bw()

top_controls <- group_by(relative_controls, Controls) %>%
                summarise( n = n() ) %>%
                arrange(desc(n)) %>%
                head(20)

relative_controls %>%
  filter( Controls %in% top_controls$Controls) %>%
  select(condition, cell_type, prop_sum_Z, prop_unique_Z, Controls) %>%
  mutate(sum_Z = prop_sum_Z + prop_unique_Z) %>%
  arrange( (sum_Z) ) %>%
  View()
  ggplot( aes(y = sum_Z, x = sum_Z, colour = Controls, label = condition )) + geom_text()



relative_controls %>%
  filter( Controls %in% top_controls$Controls) %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_prop_unique, y = rel_prop_sum, colour = Controls )) +
  #geom_text(aes(label = condition)) +
  #geom_line( aes(group = dataset), colour = "black", linetype = 2) +
  #geom_point( aes(alpha = 0.2) ) +
  geom_text( aes( label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total novel junctions\nrelative to control") +
  xlab("Proportion of unique novel junctions\nrelative to control") +
  labs(title = "Coloured by batch of control samples")

ggsave("RBPs/all_RBPs_scatter_fidelity_relative_batch.pdf")

ggplotly()

# separate by different classes of novel junction

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_anchored, y = prop_sum_anchored, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total anchored novel junctions") +
  xlab("Proportion of anchored novel junctions") +
  labs(title = "Anchored novel")

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_unanchored, y = prop_sum_unanchored, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total novel junctions") +
  xlab("Proportion of unanchored novel junctions") +
  labs(title = "Unanchored novel")

proportions_clean %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = prop_unique_skiptic, y = prop_unique_sum, colour = condition )) +
  #geom_text(aes(label = condition)) +
  geom_line( aes(group = dataset), colour = "gray", linetype = 2) +
  #geom_point( ) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) +
  geom_text( aes(label = condition), size = 2) +
  theme_bw() +
  ylab("Proportion of total skiptic junctions") +
  xlab("Proportion of skiptic junctions") +
  labs(title = "Skiptic")

 # and the relative plots:

rel_anchored <-relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_anchored, y = rel_sum_anchored, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total anchored novel junctions\nrelative to control") +
  xlab("Proportion of unique anchored novel junctions\nrelative to control") +
  labs(title = "Relative anchored novel")

rel_unanchored <- relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_unanchored, y = rel_sum_unanchored, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total unanchored novel junctions\nrelative to control") +
  xlab("Proportion of unique unanchored novel junctions\nrelative to control") +
  labs(title = "Relative unanchored ")

rel_skiptic <- relative %>%
  #filter( condition != "HNRNPC" ) %>%
  ggplot(
    aes( x = rel_unique_skiptic, y = rel_sum_skiptic, colour = condition )) +
  geom_text( aes(label = condition_kd), size = 2) +
  guides(colour = FALSE) +
  facet_wrap(~cell_type) + 
  theme_bw() +
  ylab("Proportion of total skiptic junctions\nrelative to control") +
  xlab("Proportion of unique skiptic junctions\nrelative to control") +
  labs(title = "Relative Skiptic")

ggplotly()

pdf("RBPs/all_RBPs_scatter_fidelity_relative_by_type.pdf")
print(rel_anchored)
print(rel_unanchored)
print(rel_skiptic)
dev.off()
  