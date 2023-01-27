rm( list=ls() ) 

savefile	<- 'Rdata/variances.Rdata'
if( !file.exists( savefile ) ){

	load( 'Rdata/final.Rdata' )
	imp	<- out$Y; rm( out )

	load('Rdata/setup_misc.Rdata')
	print( summary( lm( c(imp) ~ c( dat ) ) ) ) # sanity check
	imponly		<- imp
	imponly[ !is.na(dat) ]	<- NA

	miss	<- apply( is.na(dat), 2, mean )*100
	x			<- which( colnames(dat) == 'LifetimeMDD' ) 

	imp_vars	<- apply( imp			, 2, sd )
	obs_vars	<- apply( dat			, 2, sd, na.rm=T )
	impo_vars	<- apply( imponly	, 2, sd, na.rm=T ) 

	imp_corrs	<- cor( imp			, use='pairwise.complete.obs' )
	obs_corrs	<- cor( dat			, use='pairwise.complete.obs' )
	impo_corrs<- cor( imponly	, use='pairwise.complete.obs' )

	diag(imp_corrs	) <- NA
	diag(obs_corrs	) <- NA
	diag(impo_corrs	) <- NA

	save( imp_vars, obs_vars, impo_vars, imp_corrs, obs_corrs, impo_corrs, miss, x, file=savefile )
} else {
	load( savefile ) 
}

quantile( impo_vars, c(.1,.5,.9), na.rm=T )
#10%        50%        90% 
#0.08240287 0.20325341 0.35976844 
impo_vars[x] 
#  0.3196333 

quantile( (impo_corrs/obs_corrs)[x,], c(.1,.5,.9), na.rm=T )
#10%       50%       90% 
#1.351972  3.955227 12.326621 

pdf( 'figs/SuppFig1.pdf', width=5.2*2, height=5.2, useDingbats=FALSE ) 
par( mfrow=c(1,2) )
par( mar=c(4.7,4.7,.5,.5) )

plot( miss, impo_vars, xlab='Phenotype Missingness (%)', ylab='Relative Std Deviation (Imputed-Only/Observed)', type='n' )
points( miss    , impo_vars			, pch=16, cex=.7  )
points( miss[x]	, impo_vars[x]	, pch=16, cex=1.1, col=2 )
legend( 'topright', pch=16, col=2, leg=c( 'LifetimeMDD' ), bty='n' )

plot( c(-1,1), c(-1,1), xlab='Observed Correlations', ylab='Imputed-Only Correlations', type='n' )#, cex.axis=1.5, cex.lab=2.5 )
abline( a=0, b=1, col=2 )
abline( h=0, col='grey' )
abline( v=0, col='grey' )
points( obs_corrs    , impo_corrs    , pch=16, cex=.7  )
points( obs_corrs[x,], impo_corrs[x,], pch=16, cex=1.1, col=2 )
legend( 'bottomright', pch=16, col=2, leg=c( 'LifetimeMDD' ), bty='n' )#, cex=1.5 )

dev.off() 
