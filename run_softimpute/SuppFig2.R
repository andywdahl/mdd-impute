rm( list=ls() )
load('Rdata/setup_misc.Rdata')
goodnames	<- as.character( read.table( 'data/namefile.txt', head=F )[,2] )
names(goodnames)	<- as.character( read.table( 'data/namefile.txt', head=F )[,1] )

cols	<- c( '#e6beff','#fabebe', '#008080',  '#808080','#800000', '#911eb4', '#9a6324', '#ffd8b1', '#f58231', 1, '#4363d8', '#aaffc3',  '#46f0f0', '#bcf60c', '#ffe119',  '#e6194b', '#f032e6','#3cb44b','#808000','#000075', '#66ccff' )

niter	<- 10
miss	<- colMeans(is.na(dat))[sel]*100 
phens	<- colnames(dat)[sel]
np		<- length(phens)

imp_methods	<- c( 'base', 'MTAG_variables_', 'sex_strat_f', 'sex_strat_m', 'wbmi_' )

allr2s	<- array( NA, c( length(imp_methods), nf, np, niter ), dimnames=list( imp_methods, as.character(fs), phens, 1:niter )  )
for( it in 1:niter )
	for( meth in 1:5 )
{ 
	if( meth == 1 ){
		load( paste0( 'Rdata/'								, it, '_', 1, '.Rdata') )
	} else if( meth == 2 ){
		load( paste0( 'Rdata/MTAG_variables_'	, it, '_', 1, '_1', '.Rdata') )
	} else if( meth == 3 ){
		load( paste0( 'Rdata/sexstrat_'				, it, '_', 1, '_1', '.Rdata') )
	} else if( meth == 4 ){
		load( paste0( 'Rdata/sexstrat_'				, it, '_', 1, '_0', '.Rdata') )
	} else if( meth == 5 ){
		load( paste0( 'Rdata/wbmi_'						, it, '_', 1, '.Rdata') )
	}
	sub	<- intersect( dimnames( allerrs )[[3]], phens ) 
	allr2s[meth,,sub,it]	<- allerrs['cor',,sub]^2 * sign(allerrs['cor',,sub])
	rm( allerrs ) 
} 
r2s			<- 100*rowMeans( allr2s['base'						,1,,] )
r2s_MTAG<- 100*rowMeans( allr2s['MTAG_variables_'	,1,,] )
r2s_wbmi<- 100*rowMeans( allr2s['wbmi_'						,1,,] )
r2s_m		<- 100*rowMeans( allr2s['sex_strat_m'			,1,,] )
r2s_f		<- 100*rowMeans( allr2s['sex_strat_f'			,1,,] ) 

pv_fxn <- function( x, y, ... )
sapply( phens, function(p){
	tmp	<- NA
	try( tmp	<- t.test( x[p,],y[p,], ... )$p.value , silent=T )
	tmp
}) 
pvs_MTAG	<- pv_fxn( allr2s['base',1,,],allr2s['MTAG_variables_'	,1,,] )
pvs_f			<- pv_fxn( allr2s['base',1,,],allr2s['sex_strat_f'			,1,,] )
pvs_m			<- pv_fxn( allr2s['base',1,,],allr2s['sex_strat_m'			,1,,] )
pvs_wbmi	<- pv_fxn( allr2s['base',1,,],allr2s['wbmi_'						,1,,], paired=T ) ### paired bc the folds only align here; others are misaligned bc they change high-missing columns or rows

png( 'figs/SuppFig2.png', width=7.5, height=3, units='in', res=1200, pointsize=6 )#, useDingbats=FALSE ) 
layout( cbind( matrix( 1:8, 2, 4, byrow=T ), c(9,9) ) )
par( mar=c(4.3,4.3,1,1) ) 

myplot	<- function(x,y,xlab='Baseline Imputation R2 (%)',ylab='',pvs,let){ 
	plot(   x, y, main='', xlab=xlab, ylab=ylab, type='n', xlim=c(0,100), ylim=c(0,100), cex.lab=1.2 ) 
	abline(a=0,b=1,col='grey')
	bhs	<- p.adjust(pvs,'BH')
	points( x, y														, pch=16, cex=1+.9*(bhs							<.05), col='grey'		 )
	points( x[main_phens	], y[main_phens	]	, pch=16, cex=1+.9*(bhs[main_phens]	<.05), col=cols			 )
	points( x[age10_phens	], y[age10_phens]	, pch=16, cex=1+.9*(bhs[age10_phens]<.05), col=cols[4+3] )
	points( x[icd10_phens ], y[icd10_phens]	, pch=16, cex=1+.9*(bhs[icd10_phens]<.05), col=cols[5+3] ) 
	points( x[edu_phens		], y[edu_phens	]	, pch=16, cex=1+.9*(bhs[edu_phens	]	<.05), col=cols[6+3] ) 
	legend( 'topleft', bty='n', leg=let, cex=1.8, adj=1.7 )
}

myplot( r2s, r2s_MTAG, pvs=pvs_MTAG, let='A', ylab='Imputation R2 with MTAG.All Traits (%)' )
myplot( r2s, r2s_m   , pvs=pvs_m   , let='B', ylab='Imputation R2 with Males Only (%)' ) 
myplot( r2s, r2s_f   , pvs=pvs_f   , let='C', ylab='Imputation R2 with Females Only (%)' ) 
myplot( r2s, r2s_wbmi, pvs=pvs_wbmi, let='D', ylab='Imputation R2 with Baseline + BMI (%)' )

hist( pvs_MTAG, xlim=0:1, breaks=31, main='', cex.lab=1.2, xlab='P-value for: Baseline vs MTAG.All-based' ); legend( 'topright', bty='n', leg='E', cex=1.8 )
hist( pvs_m		, xlim=0:1, breaks=31, main='', cex.lab=1.2, xlab='P-value for: Baseline vs Male-only'			); legend( 'topright', bty='n', leg='F', cex=1.8 ) 
hist( pvs_f		, xlim=0:1, breaks=31, main='', cex.lab=1.2, xlab='P-value for: Baseline vs Female-only'		); legend( 'topright', bty='n', leg='G', cex=1.8 ) 
hist( pvs_wbmi, xlim=0:1, breaks=31, main='', cex.lab=1.2, xlab='P-value for: Baseline vs +BMI'						); legend( 'topright', bty='n', leg='H', cex=1.8 ) 

par( mar=c(10,0,10,0) )
plot.new()
legend( 'top', col=1, leg=c( 'BH FDR < .05', 'N.S.' ), pt.cex=c(1.9,1.0), pch=16, bty='n', cex=1.6 ) 
legend( 'bottom', fill=cols, leg=c( goodnames[main_phens], 'Age 10 Height/Weight', 'ICD10 Dep', 'Edu Attain' ), bty='n', cex=1.6 ) 

dev.off() 

## summarize diffs between joint and sex-specific for review response: 
sub	<- phens[which( abs( r2s - r2s_m ) > 5 )]
round(cbind( r2s, r2s_m, r2s_m-r2s )[sub,])

sub	<- phens[which( abs( r2s - r2s_f ) > 5 )]
round(cbind( r2s, r2s_f, r2s_f-r2s )[sub,])

''
r2s  ['LifetimeMDD']
r2s_m['LifetimeMDD']
r2s_f['LifetimeMDD']

''
r2s     ['LifetimeMDD']
r2s_MTAG['LifetimeMDD']
pvs_MTAG['LifetimeMDD']

mean(   r2s + 0*r2s_MTAG, na.rm=T )
mean( r2s_MTAG, na.rm=T )
t.test( r2s + 0*r2s_MTAG, r2s_MTAG)$p.value 

cor( r2s, r2s_wbmi, use='complete.obs' )^2

####
mytab	<- cbind( phens, round(r2s,4), round(r2s_MTAG,4), round(r2s_m,4), round(r2s_f,4), round(r2s_wbmi,4) )
colnames(mytab)	<- c( 'Phenotype', 'Baseline Imputation R2', 'MTAG.All Traits Imputation R2', 'Male-only Imputation R2', 'Female-only Imputation R2', 'Appending BMI Imputation R2' )
write.table( mytab, file='figs/SuppTab2.txt', row.names=F, quote=F, sep='\t' )
