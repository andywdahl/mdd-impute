rm( list=ls() ) 
load('Rdata/setup_misc.Rdata')

savefile	<- 'Rdata/final.Rdata'
if( !file.exists( savefile ) ){
	sink( 'Rout/final.Rout' ) 
	source( 'my_softImpute.R' )
	set.seed(1234) 
	out	<- my_softImpute( dat, lambda.len=100, nfolds=10, return_udv=TRUE )
	save( out, file=savefile ) 
	rm( out, dat ) 
	sink()
} else {
	load( savefile ) 
}

# ran once, then I manually created namefile.txt
#allnames	<- colnames(out$Y)
#write.table( cbind( allnames, allnames ), 'data/namefile_raw.txt', col.names=F, row.names=F, quote=F ) 
goodnames	<- as.character( read.table( 'data/namefile.txt', head=F )[,2] )
head( cbind( goodnames, colnames( out$Y ) ) )

if( !file.exists( 'output/ukb_md_imputed.txt' ) ){
	write.table( cbind( IIDs, FIDs, out$u )	, col.names=T, row.names=F, sep='\t', 'output/ukb_md_u.txt' )
	write.table( out$v											, col.names=T, row.names=F, sep='\t', 'output/ukb_md_v.txt' )
	write.table( cbind( IIDs, FIDs, out$Y )	, col.names=T, row.names=F, sep='\t', 'output/ukb_md_imputed.txt' )
} 

v	<- out$v
rownames(v)	<- goodnames
colnames(v)	<- paste0( 'Index ', 1:ncol(v) ) 
v	<- v * matrix( sign(as.numeric(v['Psypsy',])), nrow(v), ncol(v), byrow=T ) 
#apply( v, 2, function(x) sum(x^2) ) # all 1s

factornames	<- c( 'Neuroticism', 'Age', 'SES/EA', 'Cohabit', 'Sex/Gender' )

########################
factorplot  <- function(x,legtext,nplot=16,...){
	x <- x[ sort.list(abs(x),dec=T) ] [1:nplot]
	cols	<- 1+(sign(x)+1)/2
	x <- x^2*sign(x)*100
	plot( c(1,nplot), c(0,max(c(x,16))), type='n', ylab='% of Factor Explained', cex.lab=1.8, axes=F, xlab='', ... )
	axis(2)
	abline( h=0, col='grey' )
	box()
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
for( j in 1:10+5 )
	factorplot( v[,j], legtext=paste0( 'Factor ', j ) ) 
dev.off()
########################

########################
pdf( 'figs/Figure3a.pdf', width=5, height=5, useDingbats=FALSE ) 
par(mar=c(5,5,.5,.5))
y	<- out$d
y	<- y/sum(y)*100
plot( y, xlab='softImpute Factor (Index)', ylab='Latent Factor Strength (%)', axes=F, type='n', cex.lab=1.2 ) 
abline( h=1, col=2, lwd=1 ) 
points( y, col='grey', pch=16 ) 
axis(1); axis(2) 

sub		<- c(1:5,10,20,30)
#cols	<- rep( c( 3, 4, 1, 'purple' ), 2 )
#cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#F0E442", "#CC79A7", "#000000") ## Okabe_Ito 
cols <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#0072B2", "#000000", "#000000", "#000000") ## Okabe_Ito 
points( sub, y[sub], col=cols, pch=16 )

x	<- sub+15
x[1] <- x[1] + 15
x[2] <- x[2] + 2
x[3] <- x[3] + 8
x[4] <- x[4] + 8
x[5] <- x[5] + 15
text( x, y[sub]+c(0,0,0,.05,.015,.015,.04,.04), col=cols, c( factornames, '#10', '#20', '#30' ), cex=1.0 )

dev.off()
########################

#np	<- length(y)
#mp	<- RMTstat:::qmp( 1:np/(1+np), ndf=370*1e3, pdim=np  )
#mp	<- mp / sum(mp) * 100
#points( np:1, mp, col='grey', pch=16 )

#round( out$d, 2 )
#names(out) 
#quantile(apply( out$v, 2, mean ) )
#quantile(apply( out$v^2, 2, sum ) )
#quantile(apply( out$u, 2, mean ) )
#quantile(apply( out$u^2, 2, sum ) )
