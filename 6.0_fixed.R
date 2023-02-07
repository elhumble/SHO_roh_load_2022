library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(tidyr)
library(purrr)
source("theme_emily.R")

# Script to calculate the proportion of homozygous derived genotypes due to alleles fixed in each population

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/snpeff/freqs/*afreq data/out/6_load/DS_NS/snpeff/freqs/

#~~ Read in files
  
data_path <- "data/out/6_load/DS_NS/snpeff/freqs/"
files <- dir(data_path, pattern = "*\\.afreq")

# Visualise freqs for each locus across populations

freqs <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  mutate(filename = gsub("ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_", 
                         "", filename)) %>%
  separate(filename, c("snp_class", "Origin"), sep = "_snps_") %>%
  mutate(Origin = gsub(".afreq", "", Origin))

#~~ Alt freq distribution

ggplot(freqs, aes(ALT_FREQS)) +
  geom_histogram(bins = 15) + 
  facet_wrap(Origin~snp_class, scales = "free") +
  theme_emily()

ggplot(filter(freqs, ALT_FREQS > 0), aes(ALT_FREQS)) +
  geom_histogram(bins = 15) + 
  facet_wrap(Origin~snp_class, scales = "free") +
  theme_emily()

#~~ Read in genotype data for SNP classes

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv") %>%
  mutate(Origin = case_when(Origin == "USA (SSP)" | Origin == "USA (Ranch)" ~ "USA",
                            Origin == "EEP" ~ "EEP",
                            Origin == "EAD A" ~ "EAD A",
                            Origin == "EAD B" ~ "EAD B")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B"))) %>%
  mutate(manage = as.factor(case_when(Origin == "EEP" | Origin == "USA" ~ "Managed",
                                      TRUE ~ "Unmanaged"))) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA","EAD A", "EAD B")))


data_path <- "data/out/6_load/DS_NS/snpeff/"
files <- dir(data_path, pattern = "*\\.traw")

# Number of 0,1,2 genotypes per individual for lof and missense

df <- tibble(filename = files[c(3:4)]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  filter(name != "MSH682" & name != "MSH238" & name != "MSH250") # related and admixed

# Missense

EAD_A_missense <- filter(df, Origin == "EAD A") %>%
  filter(value == 0) %>%
  filter(filename == 
           "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, Origin == "EAD_A" & snp_class == "missense"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Variable"))

EAD_A_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EAD_B_missense <- filter(df, Origin == "EAD B") %>%
  filter(value == 0) %>%
  filter(filename == 
           "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, Origin == "EAD_B" & snp_class == "missense"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_B_missense %>% 
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)


USA_missense <- filter(df, Origin == "USA") %>%
  filter(value == 0) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, Origin == "USA" & snp_class == "missense"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

USA_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EEP_missense <- filter(df, Origin == "EEP") %>%
  filter(value == 0) %>%
  filter(filename == 
           "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, Origin == "EEP" & snp_class == "missense"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EEP_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

# LoF

EAD_A_lof<- filter(df, Origin == "EAD A") %>%
  filter(value == 0) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, Origin == "EAD_A" & snp_class == "lof"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_A_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EAD_B_lof<- filter(df, Origin == "EAD B") %>%
  filter(value == 0) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, Origin == "EAD_B" & snp_class == "lof"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_B_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)


USA_lof<- filter(df, Origin == "USA") %>%
  filter(value == 0) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, Origin == "USA" & snp_class == "lof"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

USA_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)


EEP_lof <- filter(df, Origin == "EEP") %>%
  filter(value == 0) %>%
  filter(filename == 
           "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, Origin == "EEP" & snp_class == "lof"), by = c("SNP" = "ID")) %>%
  mutate(fixed = case_when(ALT_FREQS == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EEP_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

