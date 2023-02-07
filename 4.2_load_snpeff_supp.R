# Genetic load plotting

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
library(relayer)

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/snpeff/*traw data/out/6_load/DS_NS/snpeff/

#~~ metadata

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv") %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (SSP)", "USA (Ranch)", "EAD A", "EAD B"))) %>%
  mutate(manage = as.factor(case_when(Origin == "EEP" | Origin == "USA (SSP)" | Origin == "USA (Ranch)" ~ "Managed",
                                      TRUE ~ "Unmanaged"))) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (SSP)", "USA (Ranch)", 
                                            "EAD A", "EAD B")))

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

# sho genome size

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
  #left_join(ml, by = c("name" = "ID")) %>%
  # dplyr::slice(3) %>% # alt hom == 2
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
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (SSP)", "USA (Ranch)",
                                            "EAD A", "EAD B")))


# Figure

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("GrandBudapest2")[4],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

# Absolute number of hom alternative genotypes for each SNP class

ggplot(hom_load, aes(reorder_within(Origin, total_alt_hom, snp_class), 
                     total_alt_hom, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  theme(legend.position = "none") +
  scale_x_reordered() +
  facet_wrap( ~ snp_class, scales = "free") +
  # ggtitle("Number of homozygote genotypes in each SNP class") +
  xlab("Population") + ylab("Number of homozgyote genotypes")


#~~ Correlation with inbreeding

ggplot(filter(hom_load, snp_class != "Intergenic" & snp_class != "Synonymous"), 
                    aes(froh, total_alt_hom, col = Origin)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_emily() +
  theme(strip.text = element_text(face = "plain")) +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") + 
  geom_smooth(method = "lm", se = F) + 
  scale_color_manual(values = col_palette, name = "Population") +
  ylab("Number of alternative homozygotes") + xlab(expression(italic(F["ROH"])))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Heterozygous mutation load : Number of alternative heteroyzgous genotypes      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

het_load <- df %>%
  # filter(!grepl("intergenic", filename)) %>%
  group_by(filename, name, WGS_run, Origin, manage) %>%
  mutate(total_het = case_when(value == 1 ~ n,
                               TRUE ~ 0)) %>%
  summarise(total_het = sum(total_het)) %>%
  mutate(WGS_run = case_when(WGS_run == 1 ~ "high_cov",
                             TRUE ~ "low_cov")) %>%
  #left_join(roh, by = c("name" = "iid")) %>% 
  #left_join(ml, by = c("name" = "ID")) %>%
  # dplyr::slice(3) %>% # alt hom == 2
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
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (SSP)", "USA (Ranch)",
                                            "EAD A", "EAD B"))) %>%
  mutate(snp_class = factor(snp_class, levels = c("LoF", "Missense", "Synonymous", "Intergenic")))

# Figure

ggplot(het_load, aes(reorder_within(Origin, total_het, snp_class), 
                     total_het, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_emily() +
  scale_x_reordered() +
  facet_wrap(~ snp_class, scales = "free") +
  theme(legend.position="none") +
  xlab("Population") + ylab("Number of heterozygotes") +
  scale_fill_manual(values = col_palette, name = "Population")

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
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (SSP)", "USA (Ranch)",
                                            "EAD A", "EAD B"))) %>%
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
  scale_x_reordered() +
  facet_wrap(~ snp_class, scales = "free") +
  theme(legend.position="none") +
  xlab("Population") + ylab("Total number of alternative alleles") +
  scale_fill_manual(values = col_palette, name = "Population")

#~~~~~~~~~~~~~~~~~~~~#
#      Full plot     #
#~~~~~~~~~~~~~~~~~~~~#

hom_load_fig <- ggplot(filter(hom_load, snp_class == "LoF" | 
                                 snp_class == "Missense"), 
                        # aes(reorder_within(manage, total_alt_hom, snp_class), 
                        #     total_alt_hom, fill = Location_C)) +
                        aes(factor(manage, levels = c("Managed", "Unmanaged")), 
                            total_alt_hom, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.y = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 5))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") +
  #coord_flip() +
  #ggtitle("Number of homozygote genotypes in each SNP class") +
  xlab("") + ylab("Number of \n derived homozygotes")

hom_load_fig

loc_levs <- levels(het_load$Origin)
cut.values <- setNames(col_palette[c(1,2,3,4,5)], loc_levs)

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
                    breaks = loc_levs[1:3], name = "Managed:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:3)], name = "Unmanaged:") +
  theme_emily() +
  theme(strip.text = element_text(face = "plain"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_text(margin = margin(r = 5))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF")), scales = "free") +
  xlab("Population") + ylab("Number of \n heterozygotes")

het_load_fig


total_load_fig <- ggplot(filter(mutation_load, snp_class == "LoF" | snp_class == "Missense"), 
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
                    breaks = loc_levs[1:3], name = "Managed:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:3)], name = "Unmanaged:") +
  theme_emily() +
  theme(strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = -2))) +
  scale_x_reordered() +
  facet_wrap(~ factor(snp_class, levels = c("Missense", "LoF", "Intergenic")), scales = "free") +
  xlab("Population") + ylab("Number of \n derived alleles")

total_load_fig

het_load_fig + hom_load_fig + total_load_fig + plot_layout(guides = "collect", nrow = 3,
                                          heights = c(1.6,1.6,1.6))


ggsave("figs/load_USA_split.png", het_load_fig + hom_load_fig + total_load_fig +
         plot_layout(guides = "collect", nrow = 3,
                     heights = c(1.6,1.6,1.6)),
       height = 10, width = 7)

