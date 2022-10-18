library(dplyr)
library(tidyr)
library(purrr)
library(data.table)

# Check coverage of each outgroup species alignment

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/out/6_load/DS_NS/polarise/wildebeest_topi_hartebeest/bam/cov/*txt data/out/6_load/DS_NS/polarise/wildebeest_topi_hartebeest/

data_path <- "data/out/6_load/DS_NS/polarise/wildebeest_topi_hartebeest/"
files <- dir(data_path, pattern = "*.txt")

meandepth <- tibble(filename = c(files)) %>%
  mutate(file_contents = map(filename,
                             ~ fread(file.path(data_path, .), nrows = 29, header = T))) %>%
  unnest(cols = c(file_contents))

meandepth %>%
  group_by(filename) %>%
  summarise(mean = mean(meandepth),
            min = min(meandepth),
            max = max(meandepth))
