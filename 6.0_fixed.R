library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
source("scripts/theme_emily.R")

# Script to count fixed putative deleterious alleles per pop

#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/snpeff/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.frq.strat data/out/6_load/DS_NS/snpeff/
#scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/snpeff/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.frq.strat data/out/6_load/DS_NS/snpeff/
  
col_palette <- c(wes_palette("Cavalcanti1")[2],
                   wes_palette("Cavalcanti1")[3],
                   wes_palette("Cavalcanti1")[1],
                   wes_palette("GrandBudapest1")[4])

missense_frq <- fread("data/out/6_load/DS_NS/snpeff/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.frq.strat")

hist(missense_frq$MAF)

missense_frq %>%
  group_by(CLST) %>%
  filter(MAF >0.5) %>%
  summarise(sum = n())

missense_frq %>%
  filter(MAF != 0) %>%
  group_by(CLST) %>%
  summarise(mean = mean(MAF))

ggplot(filter(missense_frq, MAF != 0), aes(MAF, fill = CLST)) +
  geom_histogram(bins = 15) + facet_wrap(~CLST) +
  theme_emily() + scale_fill_manual(values =  col_palette) +
  ggtitle("Missense")

lof_frq <- fread("data/out/6_load/DS_NS/snpeff/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.frq.strat")

hist(lof_frq$MAF)

lof_frq %>%
  group_by(CLST) %>%
  filter(MAF == 1) %>%
  summarise(sum = n())

lof_frq %>%
  filter(MAF != 0) %>%
  group_by(CLST) %>%
  summarise(mean = mean(MAF))

ggplot(filter(lof_frq, MAF != 0), aes(MAF, fill = CLST)) +
  geom_histogram(bins = 15) + facet_wrap(~CLST) +
  theme_emily() + scale_fill_manual(values =  col_palette) +
  ggtitle("LoF")

# join freq dfs

freqs <- rbind(missense_frq, lof_frq)


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

# Number of 0,1,2 genotypes per individual across SNP classes

df <- tibble(filename = files[c(2:3)]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("_.+", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  filter(name != "MSH682" & name != "MSH238" & name != "MSH250") # related and admixed

EAD_A_missense <- filter(df, Origin == "EAD A") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(missense_frq, CLST=="EAD"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Variable"))

EAD_A_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n))

EAD_B_missense <- filter(df, Origin == "EAD B") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, CLST=="Hilwa"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_B_missense %>% 
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)


USA_missense <- filter(df, Origin == "USA") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, CLST=="USA"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

USA_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EEP_missense <- filter(df, Origin == "EEP") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_missense_snps.traw") %>%
  left_join(filter(freqs, CLST=="EEP"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EEP_missense %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

# LoF

EEP_lof<- filter(df, Origin == "EEP") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, CLST=="EEP"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EEP_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

USA_lof<- filter(df, Origin == "USA") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, CLST=="USA"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

USA_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EAD_A_lof<- filter(df, Origin == "EAD A") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, CLST=="EAD"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_A_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)

EAD_B_lof<- filter(df, Origin == "EAD B") %>%
  filter(value == 2) %>%
  filter(filename == "ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_polarised_annotated_filtered_lof_snps.traw") %>%
  left_join(filter(freqs, CLST=="Hilwa"), by = "SNP") %>%
  mutate(fixed = case_when(MAF == 1 ~ "Fixed",
                           TRUE ~ "Varible"))

EAD_B_lof %>%
  group_by(fixed) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n)*100)
