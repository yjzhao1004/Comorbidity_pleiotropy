# Comorbidity_pleiotropy
## **Overview:**
These codes were used in the paper "Genome-wide Pleiotropy Analyses Identified Shared Polygenic Architecture between Somatic Conditions and Psychiatric Disorders"

* Cox models: R package survival ([https://github.com/therneau/survival](https://github.com/therneau/survival)).  
* Bayesian colocalization analysis: R package coloc v 5.2.2 ([https://cran.r-project.org/web/packages/coloc/](https://cran.r-project.org/web/packages/coloc/)).
* PLINK2: ([https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/)).
* LDSC: ([https://github.com/bulik/ldsc/](https://github.com/bulik/ldsc/)).
* FUMA ([https://fuma.ctglab.nl/](https://fuma.ctglab.nl/)).
* MAGMA ([https://ctg.cncr.nl/software/magma/](https://ctg.cncr.nl/software/magma/), also implemented in FUMA).
* CondFDR and ConjFDR ([https://github.com/precimed/pleiofdr](https://github.com/precimed/pleiofdr)).
* BOLT-LMM ([https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)).
* PRS-CS ([https://github.com/getian107/PRScs](https://github.com/getian107/PRScs))
* Custom scripts for the analyses are listed as below.

### Analysis 1: CondFDR and ConjFDR - Genome-wide Pleiotropy Search for Pleiotropic SNPs.
1. We investigated the cross-trait pleiotropy enrichment based on overlapping SNP associations between 8 somatic conditions and 10 psychiatric disorders at condFDR < 0.01. 
2. We use conjFDR to identify pleiotropic SNP in pairwise disease at conjFDR < 0.05.
   
See instruction from [https://github.com/precimed/pleiofdr](https://github.com/precimed/pleiofdr)

**Citation:**
Andreassen, O. A. et al. Improved Detection of Common Variants Associated with Schizophrenia and Bipolar Disorder Using Pleiotropy-Informed Conditional False Discovery Rate. _PLoS Genet_ 9, e1003455 (2013)

### Analysis 2: Pleiotropic Loci and Functional Annotation - FUMA 
We applied functional mapping and annotation of genetic associations (FUMA)(SNP2GENE function v1.5.4) to single-trait GWAS and significant pleiotropic variants.

See instruction from [https://fuma.ctglab.nl/snp2gene](https://fuma.ctglab.nl/snp2gene)

**Citation:**
Watanabe, K., Taskesen, E., Van Bochoven, A. & Posthuma, D. Functional mapping and annotation of genetic associations with FUMA. *Nat Commun* 8, 1826 (2017).

### Analysis 3: Bayesian Colocalization Analysis (COLOC)

To check for possible bias due to LD, we utilized the Bayesian colocalization method implemented in the R package coloc (coloc.abf function) to test the probability of colocalization for each genomic locus from FUMA. 

See install instruction from [https://cran.r-project.org/web/packages/coloc/](https://cran.r-project.org/web/packages/coloc/)

**Citation:**
Giambartolomei, C. et al. Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics. *PLoS Genet* 10, e1004383 (2014).

**coloc_manu_loci.R**: used to identify colocalized loci and best causal variants in the pleiotropy loci. 

### Analysis 4: MAGMA
We conducted a genome-wide gene-based association analysis in Multi-Marker Analysis of GenoMic Annotation (MAGMA) implemented in FUMA.

See instruction from [https://ctg.cncr.nl/software/magma](https://ctg.cncr.nl/software/magma)

**Citation:**
De Leeuw, C. A., Mooij, J. M., Heskes, T. & Posthuma, D. MAGMA: Generalized Gene-Set Analysis of GWAS Data. *PLoS Comput Biol* 11, e1004219 (2015).

**magma_gene_based_association.R**: 
We mapped GWAS SNPs from each disease to their respective protein-coding genes using MAGMA. We then mapped the boundaries of each gene to the identified pleiotropic loci to identify overlapping genes located within these pleiotropic loci.

### Analysis 5: Enrichment
See instruction of Mouse Genome Informatics platform from [http://www.informatics.jax.org/](http://www.informatics.jax.org/)

See instruction of FUMA GENE2FUNC module from [https://fuma.ctglab.nl/gene2func](https://fuma.ctglab.nl/gene2func)

**drug_enrichment.R**: Drug perturbation analysis.

**Citations:**
Blake, J. A. et al. Mouse Genome Database (MGD): Knowledgebase for mouse–human comparative biology. *Nucleic Acids Research* 49, D981–D987 (2021).

Watanabe, K., Taskesen, E., Van Bochoven, A. & Posthuma, D. Functional mapping and annotation of genetic associations with FUMA. *Nat Commun* 8, 1826 (2017).

Napolitano, F. et al. gene2drug: a computational tool for pathway-based rational drug repositioning. *Bioinformatics* 34, 1498–1505 (2018).

### Analysis 6: Survival analyses 
**cox_survival_analysis.R**: used in survival analyses to assess the utility of PRS of one disease in the prediction of the risk of developing the other disease in the pair.
