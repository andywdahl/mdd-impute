rm( list=ls() )
load('Rdata/setup_misc.Rdata')

if( !file.exists( outfile	<- 'Rdata/crossvalidated_udv.Rdata' ) ){
    ### calculate the r2 across halves of the data
	xmax	<- 190 # p=217, but only ~200 are actually returned/nonzero
	Vs		<- array( NA, dim=c( 2, ncol(dat), xmax ) )
	Ds		<- array( NA, dim=c( 2,            xmax ) )
	for( it in 1:2 ){
		load( paste0( 'Rdata/crossvalidated_2fold_', it, '.Rdata' ) ) 
		Vs[it,,]	<- out$v[,1:xmax]
		Ds[it, ]	<- out$d[ 1:xmax] 
		dimnames(Vs)[[2]]	<- colnames( out$Y ) 
		rm(out)
	} 
	save( xmax, Vs, Ds, Vs1, Ds1, file=outfile )
} else {
	load( outfile )
}

xmax	<- 50
accs	<- sapply( 1:xmax, function(x) cor( Vs[1,,x], Vs[2,,x] )^2 )

pdf( 'figs/Figure3b.pdf', width=6, height=6, useDingbats=FALSE )
par( mar=c(5,5,1,1) ) 

plot( accs, type='n', ylab='Prediction Strength (R2)', xlab='softImpute Factor (Index)' ) #, main='Across-fold R2 for V[,Index]'
abline( h=0, col='grey', lty=1 )
abline( h=1, col='grey', lty=1 )

points( accs, col='grey', pch=16, cex=.7 ) 

sub		<- c(1:5,10,20,30)
cols <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#0072B2", "#000000", "#000000", "#000000")
points( sub, accs[sub], col=cols, pch=16, cex=1.1 )

dev.off() 
