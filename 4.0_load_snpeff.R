# Genetic load plotting and analysis (SNPeff)

source("theme_emily.R")
library(dplyr)
library(data.table)
library(readxl)
library(tidyr)
library(ggplot2)
library(gghalves)
library(forcats)
library(purrr)
library(stringr)
library(wesanderson)
library(tidytext)
library(patchwork)
library(lme4)
library(broom.mixed)
library(performance)
library(sjPlot)
#devtools::install_github("clauswilke/relayer")
#remotes::install_github("clauswilke/relayer")
library(relayer)

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/snpeff/*traw data/out/6_load/DS_NS/snpeff/

#~~ metadata

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv") %>%
  mutate(Origin = case_when(Origin == "USA (SSP)" | Origin == "USA (Ranch)" ~ "USA",
                            Origin == "EEP" ~ "EEP",
                            Origin == "EAD A" ~ "EAD A",
                            Origin == "EAD B" ~ "EAD B")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B"))) %>%
  mutate(manage = as.factor(case_when(Origin == "EEP" | Origin == "USA" ~ "Managed",
                                      TRUE ~ "Unmanaged"))) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA","EAD A", "EAD B")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Load genotype files      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load impact SNP plink raw 012 files
# missense, synon, lof

data_path <- "data/out/6_load/DS_NS/snpeff/"
files <- dir(data_path, pattern = "*\\.traw")

# Number of 0,1,2 genotypes per individual across SNP classes

df <- tibble(filename = files[2:5]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  dplyr::count(name, filename, WGS_run, Origin, manage, value) %>% # count unique values of one or more vars
  group_by(name, filename) %>%
  mutate(n = as.numeric(n)) %>%
  filter(!is.na(value)) %>%
  filter(name != "MSH682" & name != "MSH238" & name != "MSH250") # related and admixed

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Load ROH data        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# SHO genome size

chr_lengths <- fread("data/meta/GCF_014754425.2_SCBI_Odam_1.1_genomic_chrom_sizes.txt") %>%
  filter(V1 != "NW_024070207.1")

# Total length of autosomal scaffolds
sho_genome_kb <- sum(chr_lengths$V2) / 1000

# roh

file_path <- "data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf.hom"

roh_all <- fread(file_path)

roh <- roh_all %>% 
  group_by(IID) %>% 
  dplyr::summarise(froh = sum(KB)/sho_genome_kb) %>% 
  dplyr::rename(iid = IID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Homozygous mutation load : Number of alternative homozygote genotypes     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Number of alternative homozygote genotypes in a given SNP class 

hom_load <- df %>%
  group_by(filename, name, WGS_run, Origin, manage) %>%
  mutate(total_alt_hom = case_when(value == 0 ~ n, 
                                   TRUE ~ 0)) %>%
  summarise(total_geno = sum(n),
            total_alt_hom = sum(total_alt_hom)) %>%
  ungroup() %>%
  group_by(name, WGS_run, Origin, manage) %>%
  mutate(sum_alt_hom = sum(total_alt_hom)) %>%
  ungroup() %>%
  group_by(filename, name, WGS_run, Origin, manage) %>%
  mutate(prop_alt_hom = total_alt_hom / sum_alt_hom) %>%
  mutate(WGS_run = case_when(WGS_run == 1 ~ "high_cov",
                             TRUE ~ "low_cov")) %>%
  left_join(roh, by = c("name" = "iid")) %>% 
  mutate(WGS_run = as.factor(WGS_run)) %>% 
  mutate(snp_class = str_split(filename, "_")) %>% 
  unnest(cols = snp_class) %>% 
  filter(snp_class %in% c("lof", "missense", "synonymous", "intergenic")) %>%
  mutate(snp_class = factor(snp_class, levels = c("missense",
                                                  "lof",
                                                  "synonymous",
                                                  "intergenic"))) %>%
  mutate(snp_class = case_when(snp_class == "lof" ~ "LoF",
                               snp_class == "missense" ~ "Missense",
                               snp_class == "synonymous" ~ "Synonymous",
                               snp_class == "intergenic" ~ "Intergenic")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B")))


# Figure

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

# Absolute number of hom alternative genotypes for each SNP class

ggplot(hom_load, aes(reorder_within(Origin, total_alt_hom, snp_class), 
                                 total_alt_hom, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  theme(legend.position = "none") +
  scale_x_reordered() +
  facet_wrap( ~ snp_class, scales = "free") +
  xlab("Population") + ylab("Number of derived homozgyote genotypes")


#~~ Models for homozygous mutation load

fit_lof <- hom_load %>% filter(snp_class == "LoF") %>%
  lm(total_alt_hom ~ manage, data = .)

tidy(fit_lof, conf.int = TRUE)
check_model(fit_lof)
plot_model(fit_lof)
summary.lm(fit_lof)

fit_missense <- hom_load %>% filter(snp_class == "Missense") %>%
  lm(total_alt_hom ~ manage, data = .)

tidy(fit_missense, conf.int = TRUE)
check_model(fit_missense)
plot_model(fit_missense)
summary.lm(fit_missense)

#~~ Correlation with inbreeding

froh_load <- ggplot(filter(hom_load, snp_class != "Intergenic" & snp_class != "Synonymous"), 
       aes(froh, total_alt_hom, col = Origin)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_emily() +
  theme(strip.text = element_text(face = "plain")) +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") + 
  geom_smooth(method = "lm", se = F) + 
  scale_color_manual(values = col_palette, name = "Population") +
  ylab("Number of alternative homozygotes") + xlab(expression(italic(F["ROH"])))

froh_load

ggsave("figs/froh_hom_load.png", froh_load, width = 7, height = 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Heterozygous mutation load : Number of heteroyzgous genotypes      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

het_load <- df %>%
  # filter(!grepl("intergenic", filename)) %>%
  group_by(filename, name, WGS_run, Origin, manage) %>%
  mutate(total_het = case_when(value == 1 ~ n,
                                  TRUE ~ 0)) %>%
  summarise(total_het = sum(total_het)) %>%
  mutate(WGS_run = case_when(WGS_run == 1 ~ "high_cov",
                             TRUE ~ "low_cov")) %>%
  mutate(WGS_run = as.factor(WGS_run)) %>% 
  mutate(snp_class = str_split(filename, "_")) %>% 
  unnest(cols = snp_class) %>% 
  filter(snp_class %in% c("lof", "missense", "synonymous", "intergenic")) %>%
  mutate(snp_class = factor(snp_class, levels = c("missense",
                                                  "lof",
                                                  "synonymous",
                                                  "intergenic"))) %>%
  mutate(snp_class = case_when(snp_class == "lof" ~ "LoF",
                               snp_class == "missense" ~ "Missense",
                               snp_class == "synonymous" ~ "Synonymous",
                               snp_class == "intergenic" ~ "Intergenic")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B"))) %>%
  mutate(snp_class = factor(snp_class, levels = c("LoF", "Missense", "Synonymous", "Intergenic")))

# Figure

ggplot(het_load, aes(reorder_within(Origin, total_het, snp_class), 
                    total_het, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =4) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_x_reordered() +
  facet_wrap(~ snp_class, scales = "free") +
  theme(legend.position="none") +
  xlab("Population") + ylab("Number of heterozygotes") +
  scale_fill_manual(values = col_palette, name = "Population")


#~~ Models

fit_lof_het <- lm(total_het ~ manage, data = het_load %>% filter(snp_class == "LoF"))
tidy(fit_lof_het, conf.int = TRUE)
check_model(fit_lof_het)
plot_model(fit_lof_het)
summary.lm(fit_lof_het)

fit_missense_het <- lm(total_het ~ manage, data = het_load %>% filter(snp_class == "Missense"))
tidy(fit_missense_het, conf.int = TRUE)
check_model(fit_missense_het)
plot_model(fit_missense_het)
summary.lm(fit_missense_het)


#~~ Correlation with inbreeding

froh_het_load <- het_load %>%
  left_join(roh, by = c("name" = "iid")) %>%
  filter(snp_class != "Intergenic" & snp_class != "Synonymous") %>%
  ggplot(aes(froh, total_het, col = Origin)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_emily() +
  theme(strip.text = element_text(face = "plain")) +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") + 
  geom_smooth(method = "lm", se = F) + 
  scale_color_manual(values = col_palette, name = "Population") +
  ylab("Number of heterozygotes") + xlab(expression(italic(F["ROH"])))

froh_het_load

ggsave("figs/froh_het_load.png", froh_het_load, width = 7, height = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Mutation load : Number of alternative deleterious alleles      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

genotyped <- df %>%
  group_by(name, Origin) %>%
  summarise(n = sum(n))

mutation_load <- df %>%
  # filter(!grepl("intergenic", filename)) %>%
  group_by(filename, name, WGS_run, Origin, manage) %>%
  mutate(total_alt = case_when(value == 0 ~ n*2,
                               value == 1 ~ n,
                               value == 2 ~ 0)) %>%
  summarise(total_alt = sum(total_alt)) %>%
  mutate(snp_class = str_split(filename, "_")) %>% 
  unnest(cols = snp_class) %>% 
  filter(snp_class %in% c("lof", "missense", "synonymous", "intergenic")) %>%
  mutate(snp_class = factor(snp_class, levels = c("missense",
                                                  "lof",
                                                  "synonymous",
                                                  "intergenic"))) %>%
  mutate(snp_class = case_when(snp_class == "lof" ~ "LoF",
                               snp_class == "missense" ~ "Missense",
                               snp_class == "synonymous" ~ "Synonymous",
                               snp_class == "intergenic" ~ "Intergenic")) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA", "EAD A", "EAD B"))) %>%
  mutate(snp_class = factor(snp_class, levels = c("LoF", "Missense", "Synonymous", "Intergenic"))) %>%
  left_join(genotyped, by = c("name", "Origin")) %>%
  mutate(prop = total_alt/n)

ggplot(mutation_load, aes(reorder_within(Origin, total_alt, snp_class), 
                          total_alt, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  facet_wrap(~ snp_class, scales = "free") +
  scale_x_reordered() +
  theme(legend.position="none") +
  xlab("Population") + ylab("Total number of alternative alleles") +
  scale_fill_manual(values = col_palette, name = "Population")


#~~ Models

fit_missense_mut <- lm(total_alt ~ manage, data = mutation_load %>% filter(snp_class == "Missense"))

tidy(fit_missense_mut, conf.int = TRUE)
check_model(fit_missense_mut)
plot_model(fit_missense_mut)
summary.lm(fit_missense_mut)

fit_lof_mut <- lm(total_alt ~ manage, data = mutation_load %>% filter(snp_class == "LoF"))

tidy(fit_lof_mut, conf.int = TRUE)
check_model(fit_lof_mut)
plot_model(fit_lof_mut)
summary.lm(fit_lof_mut)

fit_inter_mut <- lm(total_alt ~ manage, data = mutation_load %>% filter(snp_class == "Intergenic"))

tidy(fit_inter_mut, conf.int = TRUE)
check_model(fit_inter_mut)
plot_model(fit_inter_mut)
summary.lm(fit_inter_mut)


#~~ Correlation with inbreeding

froh_mutation_load <- mutation_load %>%
  left_join(roh, by = c("name" = "iid")) %>%
  filter(snp_class != "Intergenic" & snp_class != "Synonymous") %>%
  ggplot(aes(froh, total_alt, col = Origin)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_emily() +
  theme(strip.text = element_text(face = "plain")) +
  facet_wrap( ~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") + 
  geom_smooth(method = "lm", se = F) + 
  scale_color_manual(values = col_palette, name = "Population") +
  ylab("Number of derived alleles") + xlab(expression(italic(F["ROH"])))

froh_mutation_load

ggsave("figs/froh_mutation_load.png", froh_mutation_load, width = 7, height = 3)


#~~~~~~~~~~~~~#
#     Rxy     #
#~~~~~~~~~~~~~#

data_path <- "data/out/6_load/DS_NS/snpeff/"
files <- dir(data_path, pattern = "*\\.traw")

# intergenic variants for standardising

rxy_intergenic_total <- tibble(filename = files[2]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  summarise(total_alleles = n()*2)

rxy_intergenic_total_alt <- tibble(filename = files[2]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  mutate(total_alt = case_when(value == 0 ~ value*2,
                               value == 1 ~ value,
                               value == 2 ~ 0)) %>%
  summarise(total_alt = sum(total_alt))

rxy_intergenic <- left_join(rxy_intergenic_total, rxy_intergenic_total_alt, by = c("SNP", "manage")) %>%
  mutate(freq = total_alt / total_alleles) %>%
  ungroup() %>%
  select(SNP, manage, freq) %>%
  pivot_wider(names_from = manage, values_from = freq)


rxy_intergenic_ab <- rxy_intergenic %>%
  mutate(manage_unmanage = Managed * (1 - Unmanaged),
         unmanage_manage = Unmanaged * (1 - Managed)) %>%
  select(-c(SNP, Managed, Unmanaged))


rxy_intergenic_ba <- rxy_intergenic %>%
  mutate(unmanage_manage = Unmanaged * (1 - Managed),
         manage_unmanage = Managed * (1 - Unmanaged)) %>%
  select(-c(SNP, Managed, Unmanaged))


# LoF variants

rxy_lof_total <- tibble(filename = files[3]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  summarise(total_alleles = n()*2)

rxy_lof_total_alt <- tibble(filename = files[3]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  mutate(total_alt = case_when(value == 0 ~ value*2,
                               value == 1 ~ value,
                               value == 2 ~ 0)) %>%
  summarise(total_alt = sum(total_alt))

rxy_lof <- left_join(rxy_lof_total, rxy_lof_total_alt, by = c("SNP", "manage")) %>%
  mutate(freq = total_alt / total_alleles) %>%
  ungroup() %>%
  select(SNP, manage, freq) %>%
  pivot_wider(names_from = manage, values_from = freq)


rxy_lof_ab <- rxy_lof %>%
  mutate(manage_unmanage = Managed * (1 - Unmanaged),
         unmanage_manage = Unmanaged * (1 - Managed)) %>%
  select(-c(SNP, Managed, Unmanaged))

rxy_lof_ba <- rxy_lof %>%
  mutate(unmanage_manage = Unmanaged * (1 - Managed),
         manage_unmanage = Managed * (1 - Unmanaged)) %>%
  select(-c(SNP, Managed, Unmanaged))

# Missense variants

rxy_missense_total <- tibble(filename = files[4]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  summarise(total_alleles = n()*2)

rxy_missense_total_alt <- tibble(filename = files[4]) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .)))) %>%
  unnest(cols = c(file_contents)) %>%
  select(-c(CHR, `(C)M`, POS, COUNTED, ALT)) %>%
  pivot_longer(cols = -c(filename, SNP)) %>%
  mutate(name = gsub("0_", "", name)) %>%
  left_join(meta, by = c("name" = "Sample_ID")) %>%
  # right_join(ids) %>%
  group_by(SNP, manage) %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  mutate(total_alt = case_when(value == 0 ~ value*2,
                               value == 1 ~ value,
                               value == 2 ~ 0)) %>%
  summarise(total_alt = sum(total_alt))

rxy_missense <- left_join(rxy_missense_total, rxy_missense_total_alt, by = c("SNP", "manage")) %>%
  mutate(freq = total_alt / total_alleles) %>%
  ungroup() %>%
  select(SNP, manage, freq) %>%
  pivot_wider(names_from = manage, values_from = freq)

rxy_missense_ab <- rxy_missense %>%
  mutate(manage_unmanage = Managed * (1 - Unmanaged),
         unmanage_manage = Unmanaged * (1 - Managed)) %>%
  select(-c(SNP, Managed, Unmanaged))

rxy_missense_ba <- rxy_missense %>%
  mutate(unmanage_manage = Unmanaged * (1 - Managed),
         manage_unmanage = Managed * (1 - Unmanaged)) %>%
  select(-c(SNP, Managed, Unmanaged))

#~~ Create matrix

rxy_intergenic_ab_mat <- as.matrix(rxy_intergenic_ab)
rxy_intergenic_ba_mat <- as.matrix(rxy_intergenic_ba)

rxy_lof_ab_mat <- as.matrix(rxy_lof_ab)
rxy_lof_ba_mat <- as.matrix(rxy_lof_ba)

rxy_missense_ab_mat <- as.matrix(rxy_missense_ab)
rxy_missense_ba_mat <- as.matrix(rxy_missense_ba)

#~~ Run the calculation

lab_lof <- colSums(rxy_lof_ab_mat) / colSums(rxy_intergenic_ab_mat)
lba_lof <- colSums(rxy_lof_ba_mat) / colSums(rxy_intergenic_ba_mat)

lab_missense <- colSums(rxy_missense_ab_mat) / colSums(rxy_intergenic_ab_mat)
lba_missense <- colSums(rxy_missense_ba_mat) / colSums(rxy_intergenic_ba_mat)

rxy_lof <- as.matrix(c(lab_lof / lba_lof)) %>%
  as_tibble(rownames = "comparison") %>%
  dplyr::rename(Rxy = V1)

rxy_missense <- as.matrix(lab_missense / lba_missense) %>%
  as_tibble(rownames = "comparison") %>%
  dplyr::rename(Rxy = V1)

#~~ Subsampling

# subsample here 80:20, 100 times / frac = fraction of snps 

subsample <- function(iter, frac) {
  
  snps_lof <- sample(1:nrow(rxy_lof_ab_mat), round(frac * nrow(rxy_lof_ab_mat), 0), replace=TRUE)
  snps_intergenic <- sample(1:nrow(rxy_intergenic_ab_mat), round(frac * nrow(rxy_intergenic_ab_mat), 0), replace=TRUE)
  snps_missense <- sample(1:nrow(rxy_missense_ab_mat), round(frac * nrow(rxy_missense_ab_mat), 0), replace=TRUE)           
  
  # don't subsample intergenic snps [snps_intergenic, ]
  lab_lof_sub <- colSums(rxy_lof_ab_mat[snps_lof, ]) / colSums(rxy_intergenic_ab_mat)
  lba_lof_sub <- colSums(rxy_lof_ba_mat[snps_lof, ]) / colSums(rxy_intergenic_ba_mat)
  
  lab_missense_sub <- colSums(rxy_missense_ab_mat[snps_missense, ]) / colSums(rxy_intergenic_ab_mat)
  lba_missense_sub <- colSums(rxy_missense_ba_mat[snps_missense, ]) / colSums(rxy_intergenic_ba_mat)
  
  rxy_lof_sub <- as.matrix(c(lab_lof_sub / lba_lof_sub)) %>%
    as_tibble(rownames = "comparison") %>%
    dplyr::rename(Rxy = V1)
  
  rxy_missense_sub <- as.matrix(lab_missense_sub / lba_missense_sub) %>%
    as_tibble(rownames = "comparison") %>%
    dplyr::rename(Rxy = V1)
  
  out <- list(rxy_lof_sub = rxy_lof_sub, rxy_missense_sub = rxy_missense_sub)
  
}

rxy_samples <- map(1:100, subsample, 1)


#~~ Plot

rxy_samples_missense <- rxy_samples %>% 
  map("rxy_missense_sub") %>% 
  bind_rows(.id = "iter") %>% 
  group_by(comparison) %>% 
  summarise(lower = quantile(Rxy, c(0.025)), upper = quantile(Rxy, 0.975)) %>%
  dplyr::left_join(rxy_missense) %>%
  pivot_longer(cols = c(lower, upper, Rxy)) %>%
  mutate(snp_class = "Missense")

rxy_samples_lof <- rxy_samples %>% 
  map("rxy_lof_sub") %>% 
  bind_rows(.id = "iter") %>% 
  group_by(comparison) %>% 
  summarise(lower = quantile(Rxy, c(0.025)), upper = quantile(Rxy, 0.975)) %>%
  left_join(rxy_lof) %>%
  pivot_longer(cols = c(lower, upper, Rxy)) %>%
  mutate(snp_class = "LoF")

rxy_df <- rbind(rxy_samples_missense, rxy_samples_lof) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(snp_class = as.factor(snp_class))


#~~ Unmanaged versus Managed

rxy <- filter(rxy_df, comparison == "unmanage_manage") %>%
  ggplot(aes(Rxy, snp_class, col = snp_class, fill = snp_class)) +
  geom_pointrange(aes(xmax = upper, xmin = lower), 
                  size = 0.5, 
                  position=position_dodge(0.2)) +
  geom_vline(xintercept = 1) +
  theme_emily() +
  theme(legend.position = "none",
        plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 12)) +
  scale_color_manual(values = c("black", "black")) +
  xlim(0.85, 1.15) + ylab("") +
  #xlim(0.7, 1.3) + ylab("") +
  ggtitle("Unmanaged versus managed")

rxy

#~~~~~~~~~~~~~~~~~~~~#
#      Full plot     #
#~~~~~~~~~~~~~~~~~~~~#

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

hom_load_fig <- ggplot(filter(hom_load, snp_class == "LoF" | 
                                    snp_class == "Missense"),
  aes(factor(manage, levels = c("Managed", "Unmanaged")), 
      total_alt_hom, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  theme(axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = -2))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") +
  xlab("") + ylab("Number of \n derived homozygotes")

hom_load_fig


loc_levs <- levels(het_load$Origin)
cut.values <- setNames(col_palette[c(1,2,3,4)], loc_levs)

het_load_fig <- ggplot(filter(het_load, snp_class == "LoF" | snp_class == "Missense"), 
                     aes(reorder_within(manage, -total_het, snp_class),
                         total_het)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1),
                  aes(fill = Origin)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8, aes(fill = Origin)) +
  theme(legend.position="bottom") +
  geom_half_boxplot(side = "r", outlier.color = NA,
                   width = 0.6, lwd = 0.3, color = "black",
                   alpha = 0.8, aes(fill2 = Origin)) %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
  guides(fill = guide_legend(order = 1)) +
  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = loc_levs[1:2], name = "Managed:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:2)], name = "Unmanaged:") +
  theme_emily() +
  theme(strip.text = element_text(face = "plain"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text(margin = margin(r = -2))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") +
  xlab("Population") + ylab("Number of \n heterozygotes")

het_load_fig


# total load figure main text

total_load_fig <- ggplot(filter(mutation_load, snp_class == "LoF" | snp_class == "Missense"), 
                         aes(reorder_within(manage, total_alt, snp_class),
                             total_alt)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1),
                  aes(fill = Origin)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8, aes(fill = Origin)) +
  theme(legend.position="bottom") +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8, aes(fill2 = Origin)) %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
  guides(fill = guide_legend(order = 1)) +
  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = loc_levs[1:2], name = "Managed:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:2)], name = "Unmanaged:") +
  theme_emily() +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text(margin = margin(r = -2))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") +
  xlab("Population") + ylab("Number of \n derived alleles")

total_load_fig


rxy_fig <- rxy + plot_spacer() + plot_layout(widths = c(1.6,1))

#het_load_fig + hom_load_fig + rxy_fig + plot_layout(guides = "collect", nrow = 3,
                                                heights = c(1.6,1.6, 0.8))

# ggsave("figs/Figure_2.png", het_load_fig + hom_load_fig + rxy_fig + 
#          plot_layout(guides = "collect", nrow = 3,
#                      heights = c(1.6,1.6, 1)),
#        height = 8, width = 7)

het_load_fig + hom_load_fig + total_load_fig + plot_layout(guides = "collect", nrow = 3,
                                                    heights = c(1.6,1.6, 0.8))

ggsave("figs/Figure_2.png", het_load_fig + hom_load_fig + total_load_fig + 
         plot_layout(guides = "collect", nrow = 3,
                     heights = c(1.6,1.6, 1.6)),
       height = 8, width = 7)

# Rxy plot

ggsave("figs/Rxy_SNPeff.png", rxy_fig, height = 3, width = 6)

# Total load derived neutral alleles

derived_neutral <- ggplot(filter(mutation_load, snp_class == "Intergenic"), 
                         aes(reorder_within(manage, -total_alt, snp_class),
                             total_alt)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1),
                  aes(fill = Origin)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8, aes(fill = Origin)) +
  theme(legend.position="bottom") +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8, aes(fill2 = Origin)) %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
  guides(fill = guide_legend(order = 1)) +
  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = loc_levs[1:2], name = "Managed:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:2)], name = "Unmanaged:") +
  theme_emily() +
  theme(strip.text = element_text(face = "plain"),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text(margin = margin(r = 15))) +
  scale_x_reordered() +
  xlab("Population") + ylab("Number of \n derived neutral alleles")

derived_neutral

ggsave("figs/derived_alleles.png", derived_neutral,
       height = 4, width = 5)
  

