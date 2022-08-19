rm( list=ls() )
load('Rdata/setup_misc.Rdata')
goodnames	    <- as.character( read.table( 'namefile.txt', head=F )[,2] )
names(goodnames)<- as.character( read.table( 'namefile.txt', head=F )[,1] )

niter	<- 10 
miss	<- colMeans(is.na(dat))[sel]*100 
phens	<- colnames(dat)[sel]
np		<- length(phens)

errs		<- array( NA, c( 2, 2, nf, np, niter ), dimnames=list( c( 'copy', 'car' ), c( 'mse', 'cor' ), as.character(fs), phens, 1:niter )  )
for( it in 1:niter )
	for( k in 1:2 )
{
	load( paste0( 'Rdata/', it, '_', k, '.Rdata') )
	errs[k,'cor',,,it]	<- allerrs['cor',,]^2 * sign(allerrs['cor',,]) ## signed R2 reduces bias
	errs[k,'mse',,,it]	<- allerrs['mse',,]
	rm( allerrs ) 
}
age10_phens		<- phens[ grep( 'atage10'			, phens ) ] 
icd10_phens		<- phens[ grep( 'ICD10'				, phens ) ] 
edu_phens		<- phens[ grep( 'educ'				, phens ) ] 

loc		<- c( 'LifetimeMDD', 'SelfRepDep', 'DepAll', 'Psypsy', 'GPpsy', 'neuroticismscore.baseline' )
nloc	<- length(loc)
cols    <- c('#a53606', '#b32db5', '#0e288e', 'green1', 'green4', 'cyan', 'pink1', '#164c64', 'red' )

r2s		<- rowMeans( errs['copy','cor',1,,] )
r2s_car <- rowMeans( errs['car' ,'cor',1,,] )

pdf( 'figs/ExtFig1a.pdf', width=4.2, height=4.5, useDingbats=FALSE ) 
par( mar=c(4.3,4.3,1,1) ) 
plot( miss, r2s, main='', ylab='Imputation R2 (Copy-Masking)', xlab='Phenotype Missingness (%)', pch=16, type='n' ) 
box(); axis(1); axis(2)
points( miss        , r2s, pch=16, col='grey', cex=.7 )
points( miss[loc]   , r2s[loc], pch=16, col=cols, cex=1.3 )
legend( 'topright', fill=cols, leg=goodnames[loc], bty='n', cex=.88 ) 
dev.off() 

pdf( 'figs/ExtFig1b.pdf', width=5, height=5, useDingbats=FALSE ) 
par( mar=c(4.3,4.3,1,1) ) 
plot( r2s_car, r2s, main='', ylab='Imputation R2 with Random-Masking', xlab='Imputation R2 with Copy-Masking', type='n' ) 
points( r2s_car, r2s, pch=16, col='grey', cex=.8 )
abline(a=0,b=1,col='grey')
points( r2s_car[loc			], r2s[loc			], pch=16, col=cols )
points( r2s_car[age10_phens	], r2s[age10_phens  ], pch=16, col=cols[4+3] )
points( r2s_car[icd10_phens ], r2s[icd10_phens  ], pch=16, col=cols[5+3] )
points( r2s_car[edu_phens	], r2s[edu_phens    ], pch=16, col=cols[6+3] )
legend( -.03, .87, fill=cols, leg=c( goodnames[loc], 'Age 10 Height/Weight', 'ICD10 Dep', 'Edu Attain' ), bty='n', cex=.85 )
dev.off()

################## ESS
ess_gain<- rowMeans( errs['copy','cor',1,,] ) * ( N * miss/100 )
ess		<- ess_gain / ( N * (1-miss/100) ) * 100
mytab	<- cbind( phens, round(miss,1), round(r2s,3), round(ess_gain), round(ess,1) ) 
colnames(mytab)	<- c( 'Phenotype', '% Missing', 'Imputation R2', 'Imputed N_eff', '% Gain in N_eff' )
write.table( mytab, file='figs/imputation_table.txt', row.names=F, quote=F )

ess['nsamesexpartner.baseline']	<- NA   ## gives ridiculous results due to noise because highly missing (98%, next most missing is 93%), could solve noise with many niter at high computational cost
pdf( 'figs/ExtFig1c.pdf', width=3.2, height=3.2, useDingbats=FALSE ) 
par( mar=c(4.2,4.2,1,1) ) 
hist( ess, breaks=61, main='', xlab='Effective Sample Size Gain (%)', ylab='# Phenotypes' ) 
points( ess[loc], rep(-.1,nloc), pch=16, col=cols, cex=1.2 )
legend( 'topright', fill=cols[-4], leg=goodnames[loc], bty='n', cex=.9 ) 
dev.off() 
