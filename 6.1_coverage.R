library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
source("scripts/theme_emily.R")

# Sequencing depth ~ inbreeding
# samtools coverage files transferred from datastore onto local

DS_path <- "data/out/0.1_bwa/cov/DS"
DS_files <- dir(DS_path, pattern = "*cov.gz")

DS_cov <- tibble(filename = c(DS_files, DS_files)) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(DS_path, .)))) %>%
  unnest(cols = c(file_contents))

NS_path <- "data/out/0.1_bwa/cov/NS"
NS_files <- dir(NS_path, pattern = "*cov.gz")

NS_cov <- tibble(filename = c(NS_files, NS_files)) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(NS_path, .)))) %>%
  unnest(cols = c(file_contents))

cov <- rbind(NS_cov, DS_cov) %>%
  filter(
    `#rname` == "NW_024070185.1" | 
      `#rname` == "NW_024070186.1" |
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
                         TRUE ~ "RG")) %>%
  dplyr::group_by(filename) %>%
  dplyr::summarise(meandepth = mean(meandepth),
                   min = min(meandepth),
                   max = max(meandepth)) %>%
  separate(filename, c("FID"))


# Froh


#~~ Read in metadata

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv") %>%
  mutate(Origin = case_when(Origin == "USA (SSP)" | Origin == "USA (Ranch)" ~ "USA",
                            Origin == "EEP" ~ "EEP",
                            Origin == "EAD A" ~ "EAD A",
                            Origin == "EAD B" ~ "EAD B")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B"))) %>%
  mutate(manage = as.factor(case_when(Origin == "EEP" | Origin == "USA" ~ "Managed",
                                      TRUE ~ "Unmanaged"))) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA","EAD A", "EAD B")))


#~~ Chr lengths

chr_lengths <- fread("data/meta/GCF_014754425.2_SCBI_Odam_1.1_genomic_chrom_sizes.txt") %>%
  filter(V1 != "NW_024070207.1") # sex chr

# Total length of autosomal scaffolds
auto_length <- sum(chr_lengths$V2) / 1000


#~~ PLINK out

file_path <- "data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf.hom"

roh_lengths <- fread(file_path)
length(unique(roh_lengths$FID))

roh_lengths <- roh_lengths %>%
  dplyr::filter(CHR != "NW_024070207.1") %>% # remove X chr
  left_join(meta, by = c("IID" = "Sample_ID")) %>%
  mutate(WGS_run = as.factor(WGS_run))

froh <- roh_lengths %>%
  dplyr::group_by(FID, Origin, WGS_run, manage) %>%
  dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
  mutate(percent_genome = (KBSUM/auto_length) * 100) %>%
  mutate(FROH = KBSUM/auto_length)

# Plot

df <- left_join(froh, cov, by = "FID")

mod <- lm(meandepth ~ FROH, data=df)
r2 <- summary(mod)$r.squared

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

seqdepth_froh <- ggplot(df, aes(meandepth, FROH, col = Origin)) +
  geom_point(size = 2) +
  scale_color_manual(values = col_palette, name = "Population") +
  theme_emily() +
  #theme(legend.position = "none") +
  # scale_x_reordered() +
  xlab("Mean sequencing depth") + ylab("% genome in ROH") +
  facet_wrap(~Origin)

seqdepth_froh

ggsave("figs/froh_depth.png", seqdepth_froh, height = 4, width = 6)
 