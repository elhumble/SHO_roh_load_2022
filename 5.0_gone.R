library(data.table)
library(dplyr)
library(ggplot2)
library(wesanderson)
source("theme_emily.R")

# Contemporary Ne using Gone

#scp ehumble@eddie.ecdf.ed.ac.uk://exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/8_gone/DS_NS/EAD_A/outfileLD_Ne_estimates data/out/8_gone/outfileLD_Ne_estimates_EAD_A
#scp ehumble@eddie.ecdf.ed.ac.uk://exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/8_gone/DS_NS/EAD_B/outfileLD_Ne_estimates data/out/8_gone/outfileLD_Ne_estimates_EAD_B
#scp ehumble@eddie.ecdf.ed.ac.uk://exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/8_gone/DS_NS/EEP/outfileLD_Ne_estimates data/out/8_gone/outfileLD_Ne_estimates_EEP
#scp ehumble@eddie.ecdf.ed.ac.uk://exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/8_gone/DS_NS/USA/outfileLD_Ne_estimates data/out/8_gone/outfileLD_Ne_estimates_USA

#~~ Output

col_palette <- c(wes_palette("Cavalcanti1")[2],
                 wes_palette("Cavalcanti1")[3],
                 wes_palette("Cavalcanti1")[1],
                 wes_palette("GrandBudapest1")[4])

gone_ne_EAD_A <- fread("data/out/8_gone/outfileLD_Ne_estimates_EAD_A") %>%
  mutate(Origin = "EAD_A")

gone_ne_EAD_B <- fread("data/out/8_gone/outfileLD_Ne_estimates_EAD_B") %>%
  mutate(Origin = "EAD_B")

gone_ne_EEP <- fread("data/out/8_gone/outfileLD_Ne_estimates_EEP" )%>%
  mutate(Origin = "EEP")

gone_ne_USA <- fread("data/out/8_gone/outfileLD_Ne_estimates_USA" )%>%
  mutate(Origin = "USA")

gone_ne <- rbind(gone_ne_EAD_A, gone_ne_EAD_B, gone_ne_EEP, gone_ne_USA) %>%
  mutate(Origin = case_when(Origin == "EAD_A" ~ "EAD A",
                             Origin == "EAD_B" ~ "EAD B",
                             TRUE ~ Origin))


gone_ne_plot <-ggplot(gone_ne, aes(Generation, Geometric_mean, col = Origin)) +
  geom_line(size = 0.9) +
  theme_emily() +
  scale_color_manual(values = col_palette[c(3,4,1,2)]) + 
  scale_x_log10(limits= c(1, 128), breaks = c(1,2,4,8,16,32,64,128)) +
  scale_y_log10() +
  xlab("Generations back in time") +
  ylab(expression(paste("Estimated ", italic(N[e]))))

gone_ne_plot

ggsave("figs/GONe.png", gone_ne_plot, height = 4, width = 6)

