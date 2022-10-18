library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
source("scripts/theme_emily.R")
library(wesanderson)
library(stringr)
library(patchwork)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Relatedness            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in genome and NGSrelate output

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.genome data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.genome
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.res data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.res

gen <- fread("data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.genome", header = T)

summary(gen$PI_HAT)
hist(gen$PI_HAT)

#~~ Specify inference criteria

# https://academic.oup.com/bioinformatics/article/26/22/2867/228512

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1-(1/(2^(3/2))) ~ "Second-degree",
                              kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1 -(1/2^(5/2)) ~ "Third-degree",
                              kinship < 1/2^(9/2) & Z0 > 1-(1/2^(5/2)) ~ "Unrelated",
                              TRUE ~ "Unknown"))


#~~ Combine with ngsrelate output

ngsrel <- fread("data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr_mdepth_miss_maf_ld.res")

# individual order is the same
gen$R1 <- ngsrel$R1 
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING


#~~ Plotting

# Z scores

pal <- c(wes_palette("Cavalcanti1")[3],
         wes_palette("Cavalcanti1")[1],
         wes_palette("Cavalcanti1")[2],
         wes_palette("GrandBudapest1")[4],
         wes_palette("GrandBudapest2")[4])

ggplot(gen, aes(Z0, Z1)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=1)") +
  scale_colour_manual(values = pal) +
  theme_emily() +
  
  ggplot(gen, aes(Z0, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=2)") +
  ylim(c(0,1)) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  
  ggplot(gen, aes(Z1, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=1)") + ylab("Pr(IBD=2)") +
  ylim(c(0,1)) +
  scale_colour_manual(values = pal) +
  theme_emily()


# R1, R0, KING

ggplot(gen, aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.8,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Full-sibling",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.position = "none")

ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.8,
             aes(colour = factor(criteria, 
                                 levels = c("Parent-offspring",
                                            "Full-sibling",
                                            "Second-degree",
                                            "Third-degree",
                                            "Unrelated",
                                            "Unknown")))) +
  scale_colour_manual(values = pal) +
  theme_emily() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.9)) +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch"))

relate <- ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.3, shape = 21,
             aes(fill = ifelse(KING > 1/2^(5/2), "red", "black"),
             col = ifelse(KING > 1/2^(5/2), "red", "black"))) +
  scale_colour_manual(values=c("black", "red")) + 
  scale_fill_manual(values=c("black", "red")) + 
  theme_emily() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(alpha=1),
                               keyheight = 0.1,
                               default.unit = "inch")) +
  geom_hline(yintercept = 1/2^(5/2), linetype = "dashed", alpha = 0.3)

relate


#~~~ Individuals to remove

filter(gen, R1 > 0.9)

filter(gen, R1 < 0.9) %>%
  mean(PI_HAT, na.rm = T)

unrelated <- filter(gen, KING < 1/2^(5/2))
mean(unrelated$PI_HAT) 

# MSH648 (AZA) /MSH682 (Private ranch)
# Remove ranch animal MSH682 to balance populations

ggsave("figs/R1_KING.png", relate, height = 4, width = 5)

