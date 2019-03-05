library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(fidelio)
library(patchwork)

options(echo =TRUE)

# get intron database - latest GENCODE
hg38_introns <- "example/gencode.hg38.v29_all_introns.bed.gz"
#hg38_introns <- "example/gencode_hg38_all_introns.bed.gz"
#biomart_annotation <- "data/biomart_annotations_human.tab"

intron_db <- fidelio::prepareIntrons(intronList=hg38_introns)

# add biomart - get gene IDs for intron db


#files <- "example/test_junctions.txt"

#files <- "example/SNRNP70_junctions.txt"
files <- "example/junctions/DDX5_K562/DDX5_K562_junctions.txt"
files <- read.table(files, header=TRUE, stringsAsFactors=FALSE)

# add experiment column
files <- files %>%
  mutate(
    experiment_code = str_split_fixed(file_id, "_", 4)[,3],
    experiment_target = str_split_fixed(file_id, "_",4)[,1],
    replicate = str_split_fixed(file_id, "_",5)[,5] # watch out for this!
  )
file.exists(files$file)

rbp <- unique(files$experiment_target)

print(rbp)

# annotate junctions
annotated <- files[,1:2] %>%
  purrr::pmap(annotateJunctions, intron_db = intron_db)

# create proportions
annotated <- fidelio::createSimpleProportions(annotated)


# Compare filtering -> bootstrapping to filtering -> downsampling -> bootstrapping

bootstrap_plot <- function(data, mytitle ="bootstrap plot", mysubtitle = "data"){
  props <- map_df(data, "bootstrap")
  props %>%
  tidyr::gather( "key", "value", -ID) %>%
  tidyr::separate(key, into = c("method", "classification") ) %>%
  dplyr::filter( classification != "all") %>%
  filter(method == "propsum") %>%
  left_join(files, by = c("ID" = "file_id") ) %>%
  ggplot( aes( x = classification, y = value, colour = condition, group = ID) ) +
  #geom_point(position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = paste(classification,ID) ), fill = NA, notch = TRUE) +
  #scale_colour_manual( values = c("control" = "black", "SNRNP70" = "red",
  #                                "U2AF1" = "blue")) +
  #facet_wrap(~method, nrow = 2, scales = "free") +
  labs(title = mytitle,
       subtitle = mysubtitle) +
  theme_bw()
}

raw <- annotated
raw <- fidelio::bootstrapProportions(raw, times = 100, percent_genes = 0.1)
raw_50 <- fidelio::bootstrapProportions(raw, times = 100, percent_genes = 0.5)
raw_90 <- fidelio::bootstrapProportions(raw, times = 100, percent_genes = 0.9)

bootstrap_plot(raw, mysubtitle = "10% sample") +
  bootstrap_plot(raw_50, mysubtitle = "50% sample") +
  bootstrap_plot(raw_90, mysubtitle = "90% sample")

# remove intergenic junctions and any junction with count < 10
filtered <-
  annotated %>%
  fidelio::filterJunctions(min_count = 10) %>%
  fidelio::createSimpleProportions() %>%
  fidelio::bootstrapProportions(percent_genes = 0.9, times = 100)

filtered_5 <-
  annotated %>%
  fidelio::filterJunctions(min_count = 5) %>%
  fidelio::createSimpleProportions() %>%
  fidelio::bootstrapProportions(percent_genes = 0.9, times = 100)


bootstrap_plot(raw_90, mysubtitle = "raw - 90% sample") +
  bootstrap_plot(filtered_5, mysubtitle = "filtered >= 5 and bootstrapped 90%") +
  bootstrap_plot(filtered, mysubtitle = "filtered >= 10 and bootstrapped 90%")


# downsample junctions to lowest number and redo bootstrapping

downsampled <-
  annotated %>%
  fidelio::filterJunctions(min_count = 10) %>%
  fidelio::downSampleJunctions() %>%
  fidelio::bootstrapProportions(percent_genes = 0.9, times = 100)


#%>%
  #fidelio::createSimpleProportions() %>%

bootstrap_plot(raw_90, mysubtitle = "no filtering") + #ylim(0,0.006) +

bootstrap_plot(filtered, mysubtitle = "junction and intergenic filtering") + #ylim(0,0.006) +

bootstrap_plot(downsampled, mysubtitle = "filtering plus downsampling to equal depth") #+ ylim(0,0.006)

junctionANOVA <- function(data){
  bootstraps <-
    map_df(data, "bootstrap") %>%
    tidyr::gather( "key", "value", -ID) %>%
    tidyr::separate(key, into = c("method", "classification") ) %>%
    dplyr::filter( classification != "all") %>%
    filter(method == "propsum") %>%
    left_join(files, by = c("ID" = "file_id")) #%>%
    #mutate(replicate = as.factor(experiment_code))

    # fit model per junction type
    anova_res <-
      split(bootstraps, bootstraps$classification) %>%
      purrr::map_df( ~{
        res <- aov( value ~ condition + replicate, data = .x )
        res <- broom::tidy(res)
        res
      }, .id = "classification") %>%
      dplyr::filter(term == "condition")

    # calculate effect size for each junction type
    effect_sizes <-
      bootstraps %>%
      group_by(condition, classification) %>%
      summarise(mean = mean(value)) %>%
      tidyr::spread(key = condition, value = mean)
    names(effect_sizes)[3] <- "knockdown"
    effect_sizes$deltaProp <- effect_sizes[[3]] - effect_sizes[[2]]
    effect_sizes$log2Prop <- log2(effect_sizes[[3]]/effect_sizes[[2]])

    full_res <- left_join(anova_res, effect_sizes, by = "classification")
    return(full_res)
   }

junctionANOVA(raw)
junctionANOVA(filtered)
junctionANOVA(downsampled)


# bootstrap_props <- map_df(annotated,
#                                    fidelio::bootstrapProportions,
#                                    times = 100)

#bootstrap_plot <-

#bootstrap_plot_filtered

#ggsave(plot =bootstrap_plot_filtered, filename = "example/bootstrap_plot.pdf")

#filtered_bootstrap_props

#readr::write_tsv( proportion_table, "example/example_table.txt" )

# statistical testing

# test on SNRNP70

# U2AF1



# create summary tables
proportion_table <-
  purrr::map( annotated, "proportions") %>%
  bind_rows()

proportion_table %>%
  tidyr::gather( "key", "value", -ID) %>%
  tidyr::separate(key, into = c("method", "classification") ) %>%
  dplyr::filter( classification != "all") %>%
  filter(method == "propsum") %>%
  left_join(files, by = c("ID" = "file_id")) %>%
  mutate(replicate = as.factor(experiment_code)) %>%
  ggplot( aes(x  = classification, y = value)) +
  geom_point( aes(colour =condition, group = condition), position = position_dodge(width = 0.5) )

