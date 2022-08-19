# mdd-impute
Code to reproduce analyses in MDD phenotype imputation paper: "softImpute" by Dahl, ..., Cai


*************
softImpute-based code is all in run_softImpute directory:

setup_data.R and setup_misc.R -- simple scripts for QC and printing out some data summaries (output in Rout/)

estimate_imp_r2.R -- estimates imputation correlation and MSE using copy-masking or random masking
plot_imp_r2.R -- generates Extended Data Fig 1 summarizing softImpute performance

predstrength.R -- applies softImpute to two halves of the real phenotype matrix
plot_predstrength.R -- calculates and plots the r2 between factors estimated from separate halves of the phenotype matrix

final.R -- applies softImpute to the real phenotype matrix, generating the imputed phenotypes used in GWAS/PRS/etc; also generates Fig3a+c

my_softImpute.R -- softImpute wrapper. Based closely on script from Dahl et al 2016 Nature Genetics, which is based closely on a Trevor Hastie vignette: http://web.stanford.edu/~hastie/swData/softImpute/vignette.html

*************
PRS Pleiotropy code is all in PRS_pleio directory:

reshape_data.R -- transforms matrix of results for each (pheno,PRS,fold) triple into simpler matrices of size #PRSx#pheno

qc.R -- removes uninteresting secondary phenos, eg assessment center, genotype PCs

Fig6.R and Fig6_resid -- creates Figure 6 and Extended Data Figure 4
