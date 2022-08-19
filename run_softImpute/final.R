rm( list=ls() ) 
load('Rdata/setup_misc.Rdata')

################# run softImpute on real data
savefile	<- 'Rdata/final.Rdata'
if( !file.exists( savefile ) ){
	sink( 'Rout/final.Rout' ) 
	source( 'my_softImpute.R' )
	set.seed(1234) 
	out	<- my_softImpute( dat, lambda.len=100, nfolds=10, return_udv=TRUE )
	save( out, file=savefile ) 
	write.table( cbind( IIDs, FIDs, out$u )	, col.names=T, row.names=F, sep='\t', 'output/ukb_md_u.txt' )
	write.table( out$v						, col.names=T, row.names=F, sep='\t', 'output/ukb_md_v.txt' )
	write.table( cbind( IIDs, FIDs, out$Y )	, col.names=T, row.names=F, sep='\t', 'output/ukb_md_imputed.txt' )
	rm( out, dat ) 
	sink()
}
load( savefile ) 

################# plot factor loadings: v is #phenos x #factors
v	<- out$v
colnames(v)	<- paste0( 'Index ', 1:ncol(v) ) 

goodnames	<- as.character( read.table( 'namefile.txt', head=F )[,2] ) ### tidier names for UKB phenos
rownames(v)	<- goodnames

v	<- v * matrix( sign(as.numeric(v['Psypsy',])), nrow(v), ncol(v), byrow=T )  ## polarize so always has + load on Psypsy (arbitrary choice)

factornames	<- c( 'Neuroticism', 'Age', 'SES/EA', 'Cohabit', 'Sex/Gender' ) ## our interpretation of top factors

factorplot  <- function(x,legtext,nplot=16,...){
	x <- x[ sort.list(abs(x),dec=T) ] [1:nplot]
	cols	<- 1+(sign(x)+1)/2
	x <- x^2*sign(x)*100
	plot( c(1,nplot), c(0,max(c(x,16))), type='n', ylab='% of Factor Explained', cex.lab=1.8, axes=F, xlab='', ... )
	axis(2); box()
	abline( h=0, col='grey' )
	points( abs(x), col=cols, pch=16 )
	x[1] <- sign(x[1]) * min( abs(x[1])-.2, 14.5 )
	text( abs(x) + 1.7, names(x), srt=90, col=cols, cex=1.1 )
	legend( 'topright', leg=legtext, bty='n', cex=2.4 )
}

pdf( 'figs/Figure3c.pdf', width=22, height=6, useDingbats=FALSE )
par( mar=c(.5,5,.2,.2), mfrow=c(1,5) )
for( j in 1:5 )
	factorplot( v[,j], legtext=c( paste0( 'Factor ', j ), factornames[j] ) )
dev.off()

pdf( 'figs/analyze_factors_MD_UKB_6-15.pdf', width=22, height=6*2, useDingbats=FALSE )
par( mar=c(.5,5,.2,.2), mfrow=c(2,5) )
for( j in 6:15 )
	factorplot( v[,j], legtext=paste0( 'Factor ', j ) ) 
dev.off()

################# plot factor strengths (eigenvalues)
pdf( 'figs/Figure3a.pdf', width=5, height=5, useDingbats=FALSE ) 
par(mar=c(5,5,.5,.5))
y	<- out$d
y	<- y/sum(y)*100
plot( y, xlab='softImpute Factor (Index)', ylab='Latent Factor Strength (%)', axes=F, type='n', cex.lab=1.2 ) 
abline( h=1, col=2, lwd=1 ) 
points( y, col='grey', pch=16 ) 
axis(1); axis(2) 

sub		<- c(1:5,10,20,30)
cols    <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#0072B2", "#000000", "#000000", "#000000") ## Okabe_Ito 
points( sub, y[sub], col=cols, pch=16 )

x	 <- sub+15
x[1] <- x[1] + 15
x[2] <- x[2] + 2
x[3] <- x[3] + 8
x[4] <- x[4] + 8
x[5] <- x[5] + 15
text( x, y[sub]+c(0,0,0,.05,.015,.015,.04,.04), col=cols, c( factornames, '#10', '#20', '#30' ), cex=1.0 )

dev.off()
