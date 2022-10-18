# Process PCAngsd output

library(readr)
library(dplyr)
library(readxl)
library(tidyr)
library(data.table)
library(ggplot2)
library(wesanderson)
source("scripts/theme_emily.R")
library(patchwork)

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/5_angsd/DS_NS/pcangsd/ORYX_GL_DS_NS_PCAngsd_covmat.cov data/out/5_angsd/DS_NS/pcangsd/

# meta

meta <- fread("data/meta/SHO_WGS_IDs_metadata_clean.csv")

# load filenames

ids <- fread("data/out/5_angsd/DS_NS/DS_NS_rmdup_bams.txt", header = F) %>%
  mutate(V1 = gsub("/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/0.1_bwa/",
                   "",
                   V1)) %>%
  mutate(V1 = gsub("_mapped_sorted_RG","_merged", V1)) %>%
  mutate(V1 = gsub("_downsample","", V1)) %>%
  separate(V1, c("Run", "ID")) %>%
  left_join(meta, by = c("ID" = "Sample_ID")) 

head(ids)


# Load the covariance matrix

cov <- as.matrix(fread("data/out/5_angsd/DS_NS/pcangsd/ORYX_GL_DS_NS_PCAngsd_covmat.cov", header = F))
  
# Add a column with population assignments
pop <- ids$Origin

# Run PCA
mme.pca <- eigen(cov)

eigenvectors <- mme.pca$vectors #extract eigenvectors 
pca.vectors <- as.data.frame(cbind(pop, eigenvectors)) #combine with our population assignments
df <- type_convert(pca.vectors)

df$pop <- ids$Origin
df$ID <- ids$ID

# Variance explained

pca.eigenval.sum <- sum(mme.pca$values) # sum of eigenvalues
(mme.pca$values[1]/pca.eigenval.sum)*100 # Variance explained by PC1
(mme.pca$values[2]/pca.eigenval.sum)*100 # Variance explained by PC2
(mme.pca$values[3]/pca.eigenval.sum)*100 # Variance explained by PC3
(mme.pca$values[4]/pca.eigenval.sum)*100 # Variance explained by PC4


# Individuals to remove from EAD A due to high EAD B ancestry

filter(df, pop == "EAD A" & V2 < 0)

# MSH238
# MSH250

#~~ Plot PCAs and remove intermediate EAD animals

col_palette <- c(wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4],
                 wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("GrandBudapest2")[4])

pc1_2 <- ggplot(filter(df, ID != "MSH238" & ID != "MSH250"), aes(V2, V3)) + 
  geom_point(aes(colour = factor(pop)), size = 2) +
  theme_emily() +
  scale_color_manual(values = col_palette) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain"),
        legend.position = "none") +
  ylab(paste0("PC2 (",round((mme.pca$values[2]/pca.eigenval.sum)*100, digits = 2),"%)")) +
  xlab(paste0("PC1 (",round((mme.pca$values[1]/pca.eigenval.sum)*100, digits = 2),"%)")) +
  ggtitle("A") + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain"))

pc1_2

pc1_3 <- ggplot(filter(df, ID != "MSH238" & ID != "MSH250"), aes(V2, V4)) + 
  geom_point(aes(colour = factor(pop)), size = 2) +
  theme_emily() +
  scale_color_manual(values = col_palette, name = "Population") +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ylab(paste0("PC3 (",round((mme.pca$values[3]/pca.eigenval.sum)*100, digits = 2),"%)")) +
  xlab(paste0("PC1 (",round((mme.pca$values[1]/pca.eigenval.sum)*100, digits = 2),"%)")) +
  ggtitle("B") + 
  theme(plot.title=element_text(hjust=0, size = 20, face = "plain"))
 

pc1_2 + pc1_3

 
ggsave("figs/PCAngsd.png", pc1_2 + pc1_3, width = 9, height = 4)

