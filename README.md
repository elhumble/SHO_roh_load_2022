**Analysis code for:**
-------------
Humble E, Stoffel MA, Dicks K, Ball AD, Gooley RM, Chuven J, Pusey R, Al Remeithi M, Koepfli KP, Pukazhenthi B, Senn H, Ogden R: **Conservation management strategy impacts inbreeding and mutation load in scimitar-horned oryx.**  *PNAS* **120**, 18 (2023). https://doi.org/10.1073/pnas.22107561.  

**Code structure**
-------------

*Relatedness and population structure*  
`1.0_IBD.R`  
`1.1_ngsadmix_out.R`  
`1.2_pcangsd_out.R`  

*ROH analysis*   
`2.0_ROH_plink.R`  
`2.1_ROH_bcftools.R`  
`2.2_ROH_plink_supp`  

*Mutation load analysis*    
`4.0_load_snpeff.R`  
`4.1_load_vep.R`  
`4.2_load_snpeff_supp.R`  

*Ne and fixed loci*  
`5.0_gone.R`  
`6.0_fixed.R`  

*Coverage assessment*   
`0.1_vcf_cov.R`  
`3.0_outgroup_depth.R`  
`6.1_coverage.R`  

*File creation for PLINK*  
`0.0_chrom_map.R`  

**Data**
-------------
Data files can be generated using scripts from the [first part](https://github.com/elhumble/SHO_reseq_2022) of the pipeline.  

###### Please feel free to [get in touch](mailto:emily.humble@ed.ac.uk) if you have any questions.
