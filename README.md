**Chromosomal-level genome assembly of the scimitar-horned oryx: insights into diversity and demography of a species extinct in the wild**

**Summary**
-------------
This repository contains the scripts used for downstream analyses of NGDadmix, PCAngsd, ROH and deleterious mutations. Please see github repository [SHO_reseq_2022](https://github.com/elhumble/SHO_reseq_2022) for prior steps including alignment, SNP calling and variant annotation.

**Code structure**
-------------

*Scripts 0.0-0.1: Pre-processing*  
**0.0_chrom_map.R:** Creates a chrom map file for use with PLINK 
**0.1_vcf_cov.R:** Explores SNP coverage in vcf file 

*Script 1.0-1.2: Relatedness and population structure*  
**1.0_IBD.R:** Relatedness analysis 
**1.1_ngsadmix_out.R:** NGSadmix analysis and visualisation 
**1.2_pcangsd_out.R:** PCAnsgd analysis and visualisation 

*Script 2.0-2.2: ROH analysis*  
**2.0_ROH_plink.R:** ROH analysis of PLINK output  
**2.1_ROH_bcftools.R:** ROH analysis of PLINK output  
**2.2_ROH_plink_supp:** ROH analysis of PLINK output for supplementary 

*Script 3.0: Outgroup depth* 
**3.0_outgroup_depth.R:** Depth of coverage of outgroup alignments for polarisation 

*Script 4.0-4.3: Mutation load analysis*  
**4.0_load_snpeff.R:** Mutation load analysis using SNPeff annotations 
**4.1_load_vep.R:** Mutation load analysis using VEP annotations  
**4.2_load_snpeff_supp.R:** Mutation load analysis using SNPeff annotations for supplementary 

*Script 5.0: Ne* 
**5.0_gone.R:** Visualisation of GONe Ne estimates 

*Script 3.0: Fixed loci and coverage* 
**6.0_fixed.R:** Assessment of fixed deleterious loci 
**6.1_coverage.R:** Relationship between inbreeding and depth of coverage 

**Data**
-------------
Data files can be generated using scripts from the [first part](https://github.com/elhumble/SHO_reseq_2022) of the pipeline.  

###### Please feel free to [get in touch](mailto:emily.humble@ed.ac.uk) if you have any questions.
