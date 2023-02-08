# Inbreeding landscape, Ranch / AZA split

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gghalves)
library(patchwork)
library(wesanderson)
library(forcats)
library(purrr)
library(lme4)
library(broom.mixed)
library(performance)
library(sjPlot)
source("theme_emily.R")
options(scipen=999)
#install.packages("remotes")
#remotes::install_github("clauswilke/relayer")
library(relayer)
library(ggnewscale)

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf.hom data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf.hom
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_focal_bcftools_roh.txt data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_focal_bcftools_roh.txt

#~~ Read in metadata

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv") %>%
  mutate(Origin = case_when(Origin == "USA (SSP)" ~ "USA (AZA)",
                            TRUE ~ Origin)) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (AZA)", "USA (Ranch)", "EAD A", "EAD B"))) %>%
  mutate(manage = as.factor(case_when(Origin == "EEP" | Origin == "USA (AZA)" | Origin == "USA (Ranch)" ~ "Managed",
                                      TRUE ~ "Unmanaged"))) %>%
  mutate(Origin = factor(Origin, levels = c("EEP", "USA (AZA)", "USA (Ranch)", 
                                            "EAD A", "EAD B")))


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

#~~ Froh plot

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("GrandBudapest2")[4],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

froh_plot <- ggplot(froh, aes(fct_reorder(manage, percent_genome), 
                              percent_genome, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =4,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  #theme(legend.position = "none") +
  # scale_x_reordered() +
  xlab("") + ylab("% genome in ROH")  

froh_plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Specify length distribution      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

len_dis <- roh_lengths %>%
  filter(CHR != "NW_024070207.1") %>%
  mutate(length_Mb = KB/1000) %>%
  mutate(class = case_when(length_Mb >= 25.000000000 ~ 2,
                           length_Mb < 25.000000000 & length_Mb >= 12.500000000 ~ 4,
                           length_Mb < 12.500000000 & length_Mb >= 6.250000000 ~ 8,
                           length_Mb < 6.250000000 & length_Mb >= 3.125000000 ~ 16,
                           length_Mb < 3.125000000 & length_Mb >= 1.562500000 ~ 32,
                           length_Mb < 1.562500000 & length_Mb >= 0.781250000 ~ 64,
                           length_Mb < 0.781250000 & length_Mb >= 0.390625000 ~ 128)) %>%
  mutate(length_class = case_when(class == 2 ~">25.0",
                                  class == 4 ~ "12.5–25",
                                  class == 8 ~ "6.25–12.5",
                                  class == 16 ~ "3.12–6.25",
                                  class == 32 ~ "1.56–3.12",
                                  class == 64 ~ "0.78–1.56",
                                  class == 128 ~ "0.39–0.78")) %>%
  mutate(length_class = fct_reorder(length_class, class)) %>% 
  mutate(class = as.factor(class)) %>%
  select(-Origin, class) %>%
  group_by(FID, length_class, .drop = F) %>%
  dplyr::summarise(prop_IBD = sum((length_Mb / 2439) * 100)) %>%
  left_join(meta, by = c("FID" = "Sample_ID"))

hist(len_dis$prop_IBD)

# Mean IBD for each length class

len_dis %>%
  group_by(length_class) %>%
  summarise(mean = mean(prop_IBD))

# Mean IBD for each length class by Population

len_dis %>% group_by(Origin, length_class) %>%
  summarise(mean = mean(prop_IBD),
            min = min(prop_IBD),
            max = max(prop_IBD))

# Numbers of ROH in each category across population

inds <- meta %>%
  group_by(Origin) %>%
  count(name = "inds")

len_dis %>% 
  group_by(Origin, length_class) %>%
  count() %>%
  left_join(inds, by = "Origin") %>%
  mutate(prop = (n/inds) * 100)


#~~ Figure

recent_demog <- ggplot(len_dis, aes(length_class, prop_IBD, fill = Origin)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size = 2,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.8, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  scale_fill_manual(values = col_palette, name = "Population") +
  theme_emily() +
  theme(axis.ticks.x = element_line(colour = "#cccccc", size = 0.3),
        legend.position = "none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("ROH length (Mb)") + ylab("% genome in ROH")


recent_demog


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Inbreeding Ne         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


inb_ne <- function(mean_froh, class){
  
  Ne <- 1 / (2*(1 - ((1-mean_froh) ^ (1/class))))
  
}


df <- len_dis %>%
  mutate(prop_IBD = prop_IBD / 100) %>%
  ungroup() %>%
  select(-WGS_run) %>%
  mutate(length_class = as.character(length_class)) %>%
  mutate(class = case_when(length_class == "12.5–25" ~ 4,
                           length_class == "6.25–12.5" ~ 8,
                           length_class == "3.12–6.25" ~ 16,
                           length_class == "1.56–3.12" ~ 32,
                           length_class == "0.78–1.56" ~ 64,
                           length_class == "0.39–0.78" ~ 128)) %>%
  select(-length_class) %>%
  pivot_wider(names_from = class, values_from = prop_IBD) %>%
  pivot_longer(-c(FID, Origin, manage), names_to = "class", values_to = "prop_IBD") %>%
  mutate(prop_IBD = replace_na(prop_IBD, 0)) %>% # Get values for all combinations
  mutate(class = as.numeric(class))


#~~ Uncertainty

bootstrap_Ne <- function(iter, df, frac) {
  
  x <- df %>%
    group_by(class, Origin) %>%
    slice_sample(prop = frac, replace = TRUE) %>% 
    summarise(mean_froh = mean(prop_IBD)) %>%
    mutate(Ne = inb_ne(mean_froh, class))
}

uncertainty <- map(1:100, bootstrap_Ne, df, 1) %>%
  bind_rows(.id = "iter") %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>%
  drop_na(Ne) %>%
  group_by(class, Origin) %>% 
  summarise(lower = quantile(Ne, c(0.025)), upper = quantile(Ne, 0.975))


ne_plot <-  df %>% 
  group_by(class, Origin) %>%
  summarise(mean_froh = mean(prop_IBD)) %>%
  mutate(Ne = inb_ne(mean_froh, class)) %>% 
  left_join(uncertainty) %>%
  mutate(length_class = case_when(class == 2 ~ ">2",
                                  class == 4 ~ "2-4",
                                  class == 8 ~ "4-8",
                                  class == 16 ~ "8-16",
                                  class == 32 ~ "16-32",
                                  class == 64 ~ "32-64",
                                  class == 128 ~ "64-128")) %>%
  mutate(length_class = as.factor(length_class)) %>%
  #filter(class != 128) %>%
  ggplot(aes(length_class, Ne, col = Origin)) +
  #geom_point(size = 2, position = position_dodge(width=0.5)) + 
  #geom_path() +
  scale_color_manual(values = col_palette, name = "Population") +
  geom_pointrange(aes(ymax = upper, ymin = lower), size = 0.4,
                  position = position_dodge(width = 0.5)) +
  theme_emily() +
  theme(legend.position = "none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_log10(breaks = c(100,300,1000,3000)) + xlab("Generations back in time") +
  ylab(expression(italic(N[e])))


ne_plot


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Plot runs of ROH       #       
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Identify individuals with intermediate froh within each group

froh %>% arrange(FROH) %>%
  group_by(Origin) %>%
  mutate(counter = row_number())

aza_roh <- ungroup(froh) %>%
  filter(Origin == "USA (AZA)") %>%
  arrange(FROH) %>%
  slice(4:5) # middle values

ranch_roh <- ungroup(froh) %>%
  filter(Origin == "USA (Ranch)") %>%
  arrange(FROH) %>%
  slice(5:6)

eep_roh <- ungroup(froh) %>%
  filter(Origin == "EEP") %>%
  arrange(FROH) %>%
  slice(3:4)

ead_a_roh <- ungroup(froh) %>%
  filter(Origin == "EAD A") %>%
  arrange(FROH) %>%
  slice(4:5)

ead_b_roh <- ungroup(froh) %>%
  filter(Origin == "EAD B") %>%
  arrange(FROH) %>%
  slice(7:8)

rep_inds <- rbind(aza_roh, ranch_roh, eep_roh, ead_a_roh, ead_b_roh) %>%
  mutate(FID = as.factor(FID))


runs <- roh_lengths %>%
  filter(CHR != "NW_024070207.1") %>%
  mutate(
    SCAFF = case_when(
      CHR == "NW_024070185.1" ~ 1,
      CHR == "NW_024070186.1" ~ 2,
      CHR == "NW_024070187.1" ~ 3,
      CHR == "NW_024070188.1" ~ 4,
      CHR == "NW_024070189.1" ~ 5,
      CHR == "NW_024070190.1" ~ 6,
      CHR == "NW_024070191.1" ~ 7,
      CHR == "NW_024070192.1" ~ 8,
      CHR == "NW_024070193.1" ~ 9,
      CHR == "NW_024070194.1" ~ 10,
      CHR == "NW_024070195.1" ~ 11,
      CHR == "NW_024070196.1" ~ 12,
      CHR == "NW_024070197.1" ~ 13,
      CHR == "NW_024070198.1" ~ 14,
      CHR == "NW_024070199.1" ~ 15,
      CHR == "NW_024070200.1" ~ 16,
      CHR == "NW_024070201.1" ~ 17,
      CHR == "NW_024070202.1" ~ 18,
      CHR == "NW_024070203.1" ~ 19,
      CHR == "NW_024070204.1" ~ 20,
      CHR == "NW_024070205.1" ~ 21,
      CHR == "NW_024070206.1" ~ 22,
      CHR == "NW_024070207.1" ~ 23,
      CHR == "NW_024070208.1" ~ 24,
      CHR == "NW_024070209.1" ~ 25,
      CHR == "NW_024070210.1" ~ 26,
      CHR == "NW_024070211.1" ~ 27,
      CHR == "NW_024070212.1" ~ 28,
      CHR == "NW_024070213.1" ~ 29
    )
  ) %>%
  mutate(POS1 = POS1 / 1e+6,
         POS2 = POS2 / 1e+6,
         MB = KB / 1000) %>%
  dplyr::filter(FID %in% rep_inds$FID) %>%
  # filter(SCAFF %in% chr_names) %>%
  mutate(FID = factor(FID, levels = rep_inds$FID))

yax <- data.frame(FID = fct_inorder(levels(runs$FID))) %>%
  mutate(yax = seq(from = 2,
                   to = 2*length(unique(runs$FID)),
                   by = 2))

runs <- left_join(runs, yax, by = "FID")


shade <- runs %>%
  group_by(SCAFF) %>%
  summarise(min = min(POS1), max = max(POS2)) %>%
  mutate(min = case_when(SCAFF == 2 | SCAFF == 4 | SCAFF == 6 | SCAFF == 8 | SCAFF == 10 |
                           SCAFF == 12 | SCAFF == 14 | SCAFF == 16 | SCAFF == 18 | SCAFF == 20 |
                           SCAFF == 22 | SCAFF == 25 | SCAFF == 27 | SCAFF == 29 ~ 0,
                         TRUE ~ min)) %>%
  mutate(max = case_when(SCAFF == 2 | SCAFF == 4 | SCAFF == 6 | SCAFF == 8 | SCAFF == 10 |
                           SCAFF == 12 | SCAFF == 14 | SCAFF == 16 | SCAFF == 18 | SCAFF == 20 |
                           SCAFF == 22 | SCAFF == 25 | SCAFF == 27 | SCAFF == 29 ~ 0,
                         TRUE ~ max))


col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4],
                 wes_palette("GrandBudapest2")[4])

# Hacky legend

runs$Location_C <- factor(runs$Origin, levels = c("EAD B", "EAD A","EEP", "USA (AZA)", "USA (Ranch)"))

chr_names <- as.character(1:29)
names(chr_names) <- as.character(1:29)
chr_names[c(7,9,11, 13, 14, 15, 17, 19, 21, 22, 23, 25, 27, 29)] <- ""


loc_levs <- levels(runs$Location_C)
cut.values <- setNames(col_palette[c(4,3,1,2,5)], loc_levs)


roh_runs <- runs %>%
  filter(MB > 1.56) %>%
  # filter(SCAFF == 2) %>%
  ggplot() +
  geom_rect(data = shade, aes(xmin = min, xmax = max, ymin = 0, ymax = max(yax$yax)), 
            fill ='grey90', alpha=0.2) +
  geom_hline(data = yax, aes(yintercept = yax), color = "grey90", size = 0.4) +
  geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.5, fill = Location_C)) +
  geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.5, fill2 = Location_C)) %>% rename_geom_aes(new_aes = c(fill = "fill2")) +
  guides(fill=guide_legend(order = 1)) +  ##
  theme(legend.position="bottom") +
  scale_fill_manual(aesthetics = "fill", values = cut.values,
                    breaks = loc_levs[1:2], name = "Unmanaged:") +
  scale_fill_manual(aesthetics = "fill2", values = cut.values,
                    breaks = loc_levs[-(1:2)], name = "Managed:") +
  theme_emily() +
  facet_grid(~SCAFF, scales = 'free_x', space = 'free_x', switch = 'x',
             labeller = as_labeller(chr_names)) +
  theme(#strip.placement = 'outside',
    strip.text.x = element_text(face = "plain", 
                                margin = margin(b = 0, t = -1), size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0, "lines"),
    #plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
    axis.line.x = element_blank(),
    # legend.position="bottom",
    axis.title.x = element_text(margin=margin(t=-10)),
    axis.title.y = element_text(margin=margin(r=0)),
    axis.text.y = element_blank(),
    axis.line.y = element_blank()) +
  coord_cartesian(clip = 'off') +
  xlab("Chromosome") + ylab("Individuals")

roh_runs

#~~~~~~~~~~~~~~~~~~~~~~~#
#       Figure 1        #
#~~~~~~~~~~~~~~~~~~~~~~~#

aa <- froh_plot + roh_runs + 
  plot_layout(widths = c(1, 1.6))

bb <- (recent_demog + theme(legend.position = "none")) + 
  ne_plot + 
  plot_layout(widths = c(1.6, 1))


aa / bb + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


ggsave("figs/ROH_Ranch_SSP.png", aa / bb + plot_layout(guides = "collect") + 
         plot_annotation(tag_levels = "A"), height = 6, width = 9)


