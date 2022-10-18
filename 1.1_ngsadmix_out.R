# Process NGSadmix output

library(tidyverse)
library(purrr)
library(forcats)
library(readxl)
library(data.table)
library(wesanderson)
library(patchwork)
source("scripts/theme_emily.R")

# meta

meta <- read.csv("data/meta/SHO_WGS_IDs_metadata_clean.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Ln and Delta K      #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/5_angsd/DS_NS/beagle/*log data/out/5_angsd/DS_NS/ngsadmix

data_path <- "data/out/5_angsd/DS_NS/ngsadmix"
files <- dir(data_path, pattern = "*.log")

df <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             ~ readLines(file.path(data_path, .))))

get_best_like <- function(df) {
  like <- df[grep(pattern = "best like=",df)]
  like <- substr(like,regexpr(pattern="=",like)[1]+1,regexpr(pattern="=",like)[1]+14)
}

lnl <- df %>%
  mutate(lnl = map(df$file_contents, get_best_like)) %>%
  select(-file_contents) %>%
  unnest(cols = c(lnl)) %>%
  #mutate(filename = gsub("ORYX_GL_HS_NS_NGSadmix_", "", filename)) %>%
  mutate(filename = gsub("ORYX_GL_DS_NS_NGSadmix_", "", filename)) %>%
  mutate(filename = gsub("_out.log", "", filename)) %>%
  separate(filename, c("K", "run"), sep = "_") %>%
  mutate(run = gsub("run", "", run)) %>%
  mutate(run = as.numeric(run)) %>%
  mutate(K = as.factor(K)) %>%
  mutate(lnl = as.numeric(lnl)) %>%
  group_by(K) %>%
  summarise(mean = mean(lnl),
            minlnl = min(lnl),
            maxlnl = max(lnl),
            sd = sd(lnl)) %>%
  mutate(K = as.numeric(gsub("K", "", K)))

lnl_plot <- ggplot(lnl, aes(K, mean)) +
  geom_point() +
  geom_pointrange(aes(ymin = minlnl, ymax = maxlnl),
                  size = 0.5) +
  theme_emily() +
  ylab("Likelihood") +
  scale_y_continuous(labels = scales::scientific) +
  ggtitle("B")

lnl_plot

#~~ Delta K

deltaK <- lnl %>%
  mutate(LprimeK = c(NA,mean[-1]-mean[-length(mean)]),
         LdblprimeK = c(NA,LprimeK[-c(1,length(mean))]-(LprimeK)[-(1:2)],NA),
         delta = LdblprimeK/sd) 
  
deltak_plot <- ggplot(deltaK, aes(K, delta)) +
  geom_line() +
  geom_point() +
  theme_emily() +
  scale_x_continuous(limits=c(2, 5)) +
  ylab("Delta K") +
  ggtitle("A") +
  scale_y_continuous(labels = scales::scientific)


lnl_plot + deltak_plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Plotting         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/5_angsd/DS_NS/beagle/*qopt data/out/5_angsd/DS_NS/ngsadmix

data_path <- "data/out/5_angsd/DS_NS/ngsadmix/"
files <- dir(data_path, pattern = "*.qopt")

df <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .))))

# add individual names from order in which they are run in angsd

ids <- fread("data/out/5_angsd/DS_NS/DS_NS_rmdup_bams.txt", header = F) %>%
  mutate(V1 = gsub("/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/0.1_bwa/",
                   "",
                   V1)) %>%
  separate(V1, c("Run", "ID")) %>%
  left_join(meta, by = c("ID" = "Sample_ID"))

head(ids)

# Unnest data

df <- unnest(df, cols = c(file_contents)) %>%
  mutate(ID = rep(ids$ID, 60)) %>%
  pivot_longer(cols = c(V1, V2, V3, V4, V5, V6)) %>%
  drop_na() %>%
  left_join(meta, by = c("ID" = "Sample_ID"))

# Split by K

K1 <- filter(df, grepl("K1_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K2 <- filter(df, grepl("K2_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K3 <- filter(df, grepl("K3_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K4 <- filter(df, grepl("K4_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K5 <- filter(df, grepl("K5_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))

K6 <- filter(df, grepl("K6_run1_out.qopt$", filename)) %>%
  mutate(ID = as.factor(ID),
         name = as.factor(name))


# Plot

# https://luisdva.github.io/rstats/model-cluster-plots/

col_palette <- c( wes_palette("Cavalcanti1")[1],
                  wes_palette("GrandBudapest1")[4],
                  wes_palette("Cavalcanti1")[2],
                  wes_palette("Cavalcanti1")[3],
                  wes_palette("GrandBudapest2")[4],
                  wes_palette("Cavalcanti1")[4])


k2_plot <- ggplot(K2, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Origin), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=2")

k2_plot

k3_plot <- ggplot(K3, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Origin), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=3")

k3_plot

k4_plot <- ggplot(K4, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Origin), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=4")

k4_plot

k5_plot <- ggplot(K5, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Origin), switch = "x", scales = "free", space = "free") +
  theme_minimal() + 
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  ylab("K=5")

k5_plot

k6_plot <- ggplot(K6, aes(factor(ID), value, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Origin), switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  scale_fill_manual(values = col_palette) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        #axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 1, vjust = -60),
        legend.position = "none") +
  ylab("K=6")

k6_plot

k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1)

# Admixed individuals from EAD to remove

filter(K2, Origin == "EAD A" & name == "V2" & value > 0.5)

# MSH238
# MSH250

ggsave("figs/NGSadmix.png",
       k2_plot + k3_plot + k4_plot + k5_plot + k6_plot + plot_layout(ncol = 1),
       width = 7, height = 9)

ggsave("figs/deltaK.png",
       deltak_plot + lnl_plot + plot_layout(ncol = 2),
       width = 8, height = 4)
 
