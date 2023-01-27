# mdd-impute
Code to reproduce analyses in "Phenotype integration improves power and preserves specificity in biobank-based genetic studies of MDD":
https://www.biorxiv.org/content/10.1101/2022.08.15.503980v2


*************
softImpute-based code is all in run_softImpute directory:

setup_data.R and setup_misc.R -- simple scripts for QC and printing out some data summaries (output in Rout/)

run.R -- estimates imputation correlation and MSE using copy-masking or random masking
run_saveimputes.R -- saves imputed values for downstream plotting (eg SuppFig1)
run_MTAG_variables.R -- for SuppFig 2A
run_sex_strat.R -- for SuppFig 2B,C
run_wbmi.R -- for SuppFig 2D

predstrength.R -- applies softImpute to two halves of the real phenotype matrix

final.R -- applies softImpute to the real phenotype matrix, generating the imputed phenotypes used in GWAS/PRS/etc; also generates Fig3a+c

my_softImpute.R -- softImpute wrapper. Based closely on script from Dahl et al 2016 Nature Genetics, which is based closely on a Trevor Hastie vignette: http://web.stanford.edu/~hastie/swData/softImpute/vignette.html


*************
PRS Pleiotropy code is all in PRS_pleio directory:

reshape_data.R -- transforms matrix of results for each (pheno,PRS,fold) triple into simpler matrices of size #PRSx#pheno

qc.R -- removes uninteresting secondary phenos, eg assessment center, genotype PCs

Fig6.R -- creates Figure 6 and Extended Data Figure 4
