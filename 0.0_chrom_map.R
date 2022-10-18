# Make chrom map for plink editing

chrom<- fread("data/meta/GCF_014754425.2_SCBI_Odam_1.1_genomic_NAMES.txt",
              header = F)

file <- chrom %>%
  separate(V1, c("CHROM", "GENUS", "SP",
                 "ISOLATE", "CODE", "UNPLACED", "GENOMIC", "SCAFFOLD"), 
           sep = " ") %>%
  mutate(CHROM = gsub(">", "", CHROM)) %>%
  mutate(NEWCHROM = 1:nrow(.)) %>%
  mutate(NEWCHROM2 = gsub("NW_", "", CHROM)) %>%
  select(CHROM, NEWCHROM, NEWCHROM2)

write.table(file[,c(1,3)], "data/meta/chrom_map.txt", quote = F,
            row.names = F, col.names = F)


