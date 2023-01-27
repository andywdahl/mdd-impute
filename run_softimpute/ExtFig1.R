rm( list=ls() )
load('Rdata/setup_misc.Rdata')
goodnames	<- as.character( read.table( 'data/namefile.txt', head=F )[,2] )
names(goodnames)	<- as.character( read.table( 'data/namefile.txt', head=F )[,1] )

niter	<- 10 
miss	<- colMeans(is.na(dat))[sel]*100 
phens	<- colnames(dat)[sel]
np		<- length(phens)

errs	<- array( NA, c( 2, 2, nf, np, niter ), dimnames=list( c( 'copy', 'car' ), c( 'mse', 'cor' ), as.character(fs), phens, 1:niter )  )
for( it in 1:niter )
	for( k in 1:2 )
{
	load( paste0( 'Rdata/', it, '_', k, '.Rdata') )
	errs[k,'cor',,,it]	<- allerrs['cor',,]^2 * sign(allerrs['cor',,])
	errs[k,'mse',,,it]	<- allerrs['mse',,]
	rm( allerrs ) 
}
r2s			<- rowMeans( errs['copy','cor',1,,] )
r2s_car <- rowMeans( errs['car','cor',1,,] )

################## convert to ESS
N*(1-miss['LifetimeMDD']/100) #observed
# 67164 
N*(1-miss['LifetimeMDD']/100) + N*miss['LifetimeMDD']/100*r2s['LifetimeMDD'] #effective
# 173875.5 
ess_gain<- rowMeans( errs['copy','cor',1,,] ) * ( N * miss/100 )
ess			<- ess_gain / ( N * (1-miss/100) ) * 100

################## output table
mytab	<- cbind( phens, round(miss,1), round(r2s,3), round(ess_gain), round(ess,1) ) 
colnames(mytab)	<- c( 'Phenotype', '% Missing', 'Imputation R2', 'Imputed N_eff', '% Gain in N_eff' )
write.table( mytab, file='figs/imputation_table.txt', row.names=F, quote=F, sep='\t' )

################## load heldout values 
savefile	<- 'Rdata/impvec.Rdata' 
if( ! file.exists( savefile ) ){
	#load('Rdata/setup_data.Rdata')
	impmats	<- list()
	for( it in 1:niter ){
		impmats[[it]]			<- dat + NA
		load( paste0( 'Rdata/saveimp_', it, '.Rdata')  )
		print( ls() )
		imp	<- imp$Y
		imp[ is.na(dat) ]			<- NA ### remove unobserved values
		imp[ which(imp==dat) ]<- NA ### remove observed + unmasked values 
		impmats[[it]]	<- imp
		rm(imp) 
	}
	impvec	<- sapply( impmats, function(impmat) impmat[,'LifetimeMDD'] )
	mdd			<- dat[,'LifetimeMDD']
	save( mdd, impvec, file=savefile )
} else {
	load( savefile )
}

################## Figure
r2s			<- r2s*100
r2s_car	<- r2s_car*100 
impvec	<- ( impvec - min(mdd,na.rm=T) ) / ( max(mdd,na.rm=T) - min(mdd,na.rm=T) ) ### rescale to 0-1
mdd			<- ( mdd    - min(mdd,na.rm=T) ) / ( max(mdd,na.rm=T) - min(mdd,na.rm=T) )  ### rescale to 0-1 

loc		<- c( 'LifetimeMDD', 'SelfRepDep', 'DepAll', 'Psypsy', 'GPpsy', 'neuroticismscore.baseline' ) 
nloc	<- length(loc)
cols	<- c('#a53606', '#b32db5', '#0e288e', 'green1', 'green4', 'cyan', 'pink1', '#164c64', 'red' ) 

pdf( 'figs/ExtFig1.pdf', width=7.7, height=6, useDingbats=FALSE ) 
layout( cbind( c(1,3), c(2,4), 5 ), width=c(3,3,1.7) )
par( mar=c(4.3,4.3,1,1) ) 

#### A: imputation R2 vs missingness
plot( miss, r2s, main='', ylab='Imputation R2 with Copy-Masking (%)', xlab='Phenotype Missingness (%)', pch=16, type='n' ) 
box(); axis(1); axis(2)
points( miss, r2s, pch=16, col='grey', cex=.8 )
points( miss[loc], r2s[loc], pch=16, col=cols, cex=1.5 )

#### C: copy-mask vs CAR mask
plot( r2s_car, r2s, main='', xlab='Imputation R2 with Random-Masking (%)', ylab='Imputation R2 with Copy-Masking (%)', type='n' ) 
points( r2s_car, r2s, pch=16, col='grey', cex=.8 )
abline(a=0,b=1,col='grey')
points( r2s_car[loc					], r2s[loc				], pch=16, cex=1.5, col=cols )
points( r2s_car[age10_phens	], r2s[age10_phens], pch=16, cex=1.5, col=cols[4+3] )
points( r2s_car[icd10_phens ], r2s[icd10_phens], pch=16, cex=1.5, col=cols[5+3] ) 
points( r2s_car[edu_phens		], r2s[edu_phens	], pch=16, cex=1.5, col=cols[6+3] ) 

#### B: Dist'n of gain in ESS across traits
ess1		<- ess; ess1['nsamesexpartner.baseline']	<- NA
hist( ess1, breaks=61, main='', xlab='Effective Sample Size Gain (%)', ylab='# Phenotypes' ) 
points( ess1[loc], rep(-.1,nloc), pch=16, col=cols, cex=1.5 )

#### D: Dist'n of imputed LifetimeMDD for held-out cases and controls
breaks<- hist( c(impvec,mdd), plot=F, breaks=21 )$breaks
mids	<- hist( c(impvec,mdd), plot=F, breaks=21 )$mids
d1		<- hist( impvec[which(mdd==min(mdd,na.rm=T)),], breaks=breaks, plot=F )$density
d2		<- hist( impvec[which(mdd==max(mdd,na.rm=T)),], breaks=breaks, plot=F )$density

plot( range(breaks), range(c(d1,d2)), type='n', xlab='Imputed Values of Masked LifetimeMDD', ylab='Density' )
lines( mids, d1, col=1, lwd=1.4 )
lines( mids, d2, col=2, lwd=1.4 )
legend( 'topright', fill=1:2, col=1:2, leg=c( 'Controls', 'Cases' ), bg='white' ) 

par( mar=rep(0,4) )
plot.new()
legend( 'center', fill=cols, leg=c( goodnames[loc[-length(loc)]], 'ICD10 Dep', goodnames[loc[length(loc)]], 'Age 10 Height/Weight', 'Edu Years' ), bty='n', cex=1.1 ) 

dev.off()

##### threshold to binary 0/1 label for review response 
mat	<- matrix( 
c(
sum( 1-mdd[ impvec < .5 ] , na.rm=T ),
sum( 1-mdd[ impvec > .5 ] , na.rm=T ),
sum(   mdd[ impvec < .5 ] , na.rm=T ),
sum(   mdd[ impvec > .5 ] , na.rm=T )
)
,2,2)
colnames(mat)	<- c('Obs Ctrl', 'Obs Case')
rownames(mat)	<- c('Imp Ctrl', 'Imp Case')
mat
#              Obs Ctrl  Obs Case
# Imp Ctrl     2202      611
# Imp Case       10      121
