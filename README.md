# Comorbidity_pleiotropy
## **Overview:**
These codes were used in the paper "Genome-wide Pleiotropy Analyses Identified Shared Polygenic Architecture between Somatic Conditions and Brain Disorders"

* Competing risk models: R package cmprsk v 2.2-11 ([https://cran.r-project.org/web/packages/cmprsk/index.html](https://cran.r-project.org/web/packages/cmprsk/index.html)). 
* Two sample MR analysis: R package TwoSampleMR ([https://mrcieu.github.io/TwoSampleMR/](https://mrcieu.github.io/TwoSampleMR/)) and RadialMR ([https://github.com/WSpiller/RadialMR/](https://github.com/WSpiller/RadialMR/)). 
* Bayesian colocalization analysis: R package coloc v 5.2.2 ([https://cran.r-project.org/web/packages/coloc/](https://cran.r-project.org/web/packages/coloc/)). 
* LDSC: ([https://github.com/bulik/ldsc/](https://github.com/bulik/ldsc/))
* PLACO ([https://github.com/RayDebashree/PLACO](https://github.com/RayDebashree/PLACO))
* FUMA ([https://fuma.ctglab.nl/](https://fuma.ctglab.nl/))
* MAGMA ([https://ctg.cncr.nl/software/magma/](https://ctg.cncr.nl/software/magma/), also implemented in FUMA), 
* Custom scripts for the analyses are listed as below.
* We provided the data example in the link: [[https://drive.google.com/drive/my-drive](https://drive.google.com/drive/folders/1EYpq8FZ8N5lkFfsB0z3sR-4yfQ4-1Hlb?usp=sharing)](https://drive.google.com/drive/folders/1EYpq8FZ8N5lkFfsB0z3sR-4yfQ4-1Hlb?usp=sharing)

### Analysis 1: Cross-trait LDSC

We performed the cross-trait LDSC analyses on a total of 135 (9 × 15) pairs of diseases. LDSC requires the two single-trait GWAS summary statistics as input and estimates the pairwise genetic correlations (rgs) between pairs of traits. 

See install instruction from [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)

**Citation:**
Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. *Nature Genetics*, 2015.

* Prepare GWAS summary statistics and LD score
    * See *Methods-Data Sources for GWAS section* for sources and quality-control steps in GWAS data preparations
* munge_sumstats.py
    > python /Users/yujiezhao/Desktop/COMMO/LDSC/ldsc-master/munge_sumstats.py --sumstats /Users/yujiezhao/Desktop/COMMO/LDSC/data/QC_GWAS_new/D1_gwas_qc.txt --N-col N --chunksize 500000 --out /Users/yujiezhao/Desktop/COMMO/LDSC/data/sum_gwas/D1 --merge-alleles /Users/yujiezhao/Desktop/COMMO/LDSC/w_hm3.noMHC.snplist

* ldsc.py
    > python ldsc.py \
  --rg /Users/yujiezhao/Desktop/COMMO/LDSC/data/sum_gwas/D1.sumstats.gz,/Users/yujiezhao/Desktop/COMMO/LDSC/data/sum_gwas/D2.sumstats.gz \
  --ref-ld-chr /Users/yujiezhao/Desktop/COMMO/LDSC/eur_w_ld_chr/ \
  --w-ld-chr /Users/yujiezhao/Desktop/COMMO/LDSC/eur_w_ld_chr/  \
  --out /Users/yujiezhao/Desktop/COMMO/LDSC/data/cross_genecorr/D1_D2
  
**LDSC_Bubble_plot.R**: Fig.2 plot

### Analysis 2: PLACO 
We applied the pleiotropic analysis under the composite null hypothesis (PLACO) to conduct a genome-wide search for pleiotropic SNPs influencing risks of all pairwise diseases with significant genetic correlation. 

See install instruction from [https://github.com/RayDebashree/PLACO](https://github.com/RayDebashree/PLACO)

**Citation:**
Ray, D., Chatterjee, N. (2020) "A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between Type 2 Diabetes and Prostate Cancer". *PLoS Genetics* 16(12): e1009218.

**cal_placo.R**: used to calculate Placo.p.value for each SNP shared in the disease pairs with significant genetic correlation and draw qq-plots. (we provided T2D.sumstats.gz and ALD.sumstats.gz for data example)

### Analysis 3: FUMA 
We applied functional mapping and annotation of genetic associations (FUMA)(SNP2GENE function v1.5.4) to single-trait GWAS and significant pleiotropic variants from PLACO for **<u>pleiotropic loci characterization</u>** and **<u>functional annotation</u>**. 

See instruction from [https://fuma.ctglab.nl/snp2gene](https://fuma.ctglab.nl/snp2gene)

**Citation:**
Watanabe, K., Taskesen, E., Van Bochoven, A. & Posthuma, D. Functional mapping and annotation of genetic associations with FUMA. *Nat Commun* 8, 1826 (2017).

**fuma_pleiotropy_loci.R**: used to select pleiotropy loci for disease pairs and collect functional annotation results.(take T2D and ALD for example)

**eqtlmapping_collection.R**: used to collect eQTL mapping results in these pleiotropy loci. eQTL mapping provides significant SNP-gene-tissue pairs in brain tissue types from GTEx/v8 project. No additional variant filtering by functional annotation was applied in the eQTL mapping.

### Analysis 4: Bayesian Colocalization Analysis (COLOC)

To check for possible bias due to LD, we utilized the Bayesian colocalization method implemented in the R package coloc (coloc.abf function) to test the probability of colocalization for each genomic locus from FUMA. 

See install instruction from [https://cran.r-project.org/web/packages/coloc/](https://cran.r-project.org/web/packages/coloc/)

**Citation:**
Giambartolomei, C. et al. Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics. *PLoS Genet* 10, e1004383 (2014).

**coloc_manu_loci.R**: used to identify colocalized loci and best causal variants in the pleiotropy loci. (take T2D and ALD for example)

### Analysis 5: MAGMA
We conducted a genome-wide gene-based association analysis in Multi-Marker Analysis of GenoMic Annotation (MAGMA) implemented in FUMA.

See instruction from [https://ctg.cncr.nl/software/magma](https://ctg.cncr.nl/software/magma)

**Citation:**
De Leeuw, C. A., Mooij, J. M., Heskes, T. & Posthuma, D. MAGMA: Generalized Gene-Set Analysis of GWAS Data. *PLoS Comput Biol* 11, e1004219 (2015).

**magma_gene_based_association.R**: 
For each pleiotropic loci of disease pairs, we identified the genes significantly associated with the pleiotropic SNPs based on PLACO outputs, after Bonferroni correction for the number of genes tested. We further subset the genes significantly associated with pleiotropic SNPs based on both two GWAS results of pairwise diseases. Those genes with significance declared at loci-specific Bonferroni corrected 2-side P < 0.05 based on PLACO results, and 2-side P < 0.05 based on two GWAS were identified as candidate pleiotropic genes.(take T2D and ALD for example)

**manhattanplot.R**: Fig 3 plot

### Analysis 6: AHBA expression map

Brain transcriptional data was from the Allen Human Brain Atlas (AHBA) ([http://human.brain- map.org/](http://human.brain-%20map.org/)).

AHBA gene expression data was preprocessed according to a recommended pipeline ([https://github.com/BMHLab/AHBAprocessing](https://github.com/BMHLab/AHBAprocessing)). 

**Citations**:
Hawrylycz, M. J. et al. An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature* 489, 391–399 (2012).
Arnatkevic̆iūtė, A., Fulcher, B. D. & Fornito, A. A practical guide to linking brain-wide gene expression and neuroimaging data. *NeuroImage* 189, 353–367 (2019).

**gem_ahba.m**: used to examine the mean normalized gene expression map and plot Fig.4A.


### Analysis 7: Enrichment
See instruction of Mouse Genome Informatics platform from [http://www.informatics.jax.org/](http://www.informatics.jax.org/)

See instruction of FUMA GENE2FUNC module from [https://fuma.ctglab.nl/gene2func](https://fuma.ctglab.nl/gene2func)

**Citations:**
Blake, J. A. et al. Mouse Genome Database (MGD): Knowledgebase for mouse–human comparative biology. *Nucleic Acids Research* 49, D981–D987 (2021).

Watanabe, K., Taskesen, E., Van Bochoven, A. & Posthuma, D. Functional mapping and annotation of genetic associations with FUMA. *Nat Commun* 8, 1826 (2017).

**MGI_phenotype_enrichment.R**: used to perform phenotype-specific analysis and plot Fig.4B

**FUMA_tissue_pathway_specificity.R**: used to extract results of tissue-specific analysis and gene-set enrichment analysis obtained from FUMA (GENE2FUNC module).

### Analysis 8: Mendelian Randomization(MR) analysis

**cal_MR_pairdiseases.R**: used to perform MR analyses in pairwise diseases to test the bidirectional causality.

### Analysis 9: Comorbid Association Analysis (chi-squared test)
**chisq_test_comorbid.R**: used to examine the comorbid association between two diseases in UKB by comparing the difference in the proportion of participants in the four groups whether or not they were diagnosed with these two types of diseases. 

### Analysis 10: Competing Risk Models 
**predict_PRS_cmprsk.R**: used in survival analyses to assess the utility of PRS in the prediction of developing comorbid disorders.
