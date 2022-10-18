library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Prepare files for downsampling high coverage individuals

options(scipen = 999)

# Read in output of samtools coverage

data_path <- "data/out/0.1_bwa/cov/HS"
files <- dir(data_path, pattern = "*cov.gz")

cov <- tibble(filename = c(files, files)) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents))

# extract chr length scaffs

cov <- cov %>%
  filter(
    `#rname` == "NW_024070185.1" | `#rname` == "NW_024070186.1" |
      `#rname` == "NW_024070187.1" |
      `#rname` == "NW_024070188.1" |
      `#rname` == "NW_024070189.1" |
      `#rname` == "NW_024070190.1" |
      `#rname` == "NW_024070191.1" |
      `#rname` == "NW_024070192.1" |
      `#rname` == "NW_024070193.1" |
      `#rname` == "NW_024070194.1" |
      `#rname` == "NW_024070195.1" |
      `#rname` == "NW_024070196.1" |
      `#rname` == "NW_024070197.1" |
      `#rname` == "NW_024070198.1" |
      `#rname` == "NW_024070199.1" |
      `#rname` == "NW_024070200.1" |
      `#rname` == "NW_024070201.1" |
      `#rname` == "NW_024070202.1" |
      `#rname` == "NW_024070203.1" |
      `#rname` == "NW_024070204.1" |
      `#rname` == "NW_024070205.1" |
      `#rname` == "NW_024070206.1" |
      `#rname` == "NW_024070207.1" |
      `#rname` == "NW_024070208.1" |
      `#rname` == "NW_024070209.1" |
      `#rname` == "NW_024070210.1" |
      `#rname` == "NW_024070211.1" |
      `#rname` == "NW_024070212.1" |
      `#rname` == "NW_024070213.1"
  ) %>%
  mutate(bam = case_when(grepl("rmdup", `filename`) ~ "rmdup",
                         TRUE ~ "RG"))


hist(cov$meandepth) # depth
hist(cov$coverage) # breadth


meta <- read.csv("data/meta/SHO_WGS_IDs_metadata_clean.csv")

cov <- cov %>%
  mutate(filename = gsub("_mapped_sorted_RG_rmdup_cov.gz", "", filename)) %>%
  mutate(filename = gsub("_mapped_sorted_RG_cov.gz", "", filename)) %>%
  mutate(filename = gsub("_merged_rmdup_cov.gz", "", filename)) %>%
  mutate(filename = gsub("_merged_cov.gz", "", filename)) %>%
  right_join(meta, by = c("filename" = "Sample_ID"))

# Mean depth per sample

cov %>%
  dplyr::group_by(filename, WGS_run) %>%
  dplyr::summarise(depth = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth))

# Mean depth across wgs run and rmdup

cov %>%
  dplyr::group_by(filename, WGS_run, bam) %>%
  dplyr::summarise(mean = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth)) %>%
  ungroup() %>%
  dplyr::group_by(WGS_run, bam) %>%
  dplyr::summarise(meandepth = mean(mean),
                   min = min(min),
                   max = max(max))


#~~~ Calculating proportion to subsample

# Average low coverage is ~ 6X post rmdup

subsample <- cov %>%
  dplyr::group_by(filename, WGS_run) %>%
  dplyr::summarise(depth = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth)) %>%
  filter(WGS_run == 1) %>%
  mutate(prop_sub = 7/depth) %>%
  mutate(prop_sub = round(prop_sub, 2)) %>%
  mutate(prop_sub = gsub("0.", "1.", prop_sub)) # replacing 0 with 1 for samtools view -s syntax

write.table(subsample[c(1,6)], "data/out/0.1_bwa/cov/HS/proportion_to_downsample.txt",
            quote = F, col.names = F, row.names = F)

# move onto eddie

# scp data/out/0.1_bwa/cov/proportion_to_downsample.txt ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/file_lists/ 
  
  
