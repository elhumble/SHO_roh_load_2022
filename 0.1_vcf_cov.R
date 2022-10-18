library(data.table)
library(dplyr)
library(tidyr)
source("scripts/theme_emily.R")

# (future.globals.maxSize = 4000 * 1024^5)

# VCFtools depth output on vcf files

#~~~~~~~~~~~~~~~~~~~~~#
#       DS & NS       #
#~~~~~~~~~~~~~~~~~~~~~#

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr.gdepth data/out/4_vcftools/DS_NS/ORYX_geno_DS_NS_biallelic_chr.gdepth
# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_pop_gen_2020/data/out/4_IBD_IBC/DS_NS/ORYX_geno_DS_NS_biallelic_chr.idepth data/out/4_vcftools/DS_NS/ORYX_geno_DS_NS_biallelic_chr.idepth

#~~ Depth

idepth <- fread("data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr.idepth", colClasses = "numeric")
summary(idepth$MEAN_DEPTH)


depth <- fread("data/out/4_ibd_ibc/DS_NS/ORYX_geno_DS_NS_biallelic_chr.gdepth", colClasses = "numeric") %>%
  unite(SNP, CHROM, POS)

View(depth[1:100,1:21])

depth <- select(depth, -SNP)

depth <- depth %>%
  mutate(mean = rowMeans(.), 
         sum = rowSums(.), 
         logmeandepth = log10(mean))

summary(depth$mean)

# DS_NS
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    5.215    5.646    5.792    6.123 3032.692

CI <- 0.95
CI_meanDP <- stats::quantile(depth$mean, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_sumDP <- stats::quantile(depth$sum, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
CI_meanDP
CI_sumDP

deptha <- filter(depth, mean < 40)
hist(deptha$mean, breaks = 100)

ggplot(deptha, aes(mean)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 5)



