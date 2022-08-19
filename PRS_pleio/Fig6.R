rm( list=ls() )
load( file='data/allmaplesyrup_qcd.Rdata' )
rownames(y)[9] <- 'MTAG.All'

js	<- c( 4, 3, 1, 2, 10, 11 )#, 3 )#, 5 )
minps   <- apply( pvs, 2, min )
sub     <- which( minps < .05/ncol(pvs) )
y       <- y[,sub]
nphen   <- length(sub)
phens   <- phens[sub]

y	<- y*100
js0 <- c(4,3,12)

pdf( 'figs/Fig6.pdf', width=12, height=4, useDingbats=FALSE )
par( mfrow=c(1,3) )
par( mar=c(4.2,4.2,.7,.7) )

###### Panel A
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2 )
for( i in js0 )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )

sel         <- phens[ sort.list( y['GPpsy',], dec=T )[1:5] ]
goodnames   <- sel
names(goodnames)   <- sel
goodnames[4]    <- 'Neurot'
y.off   <- c( 0, 0, 5, -3, 0 )
x0s     <- c( 10, 15, 20, 20, 30 )
for( p in sel ){
    x0  <- x0s[which(sel==p)]+10
    y0  <- y['GPpsy',p]+y.off[which(sel==p)]
    text( x0 + ifelse( p == 'SelfRepDep', 2, 0 ), y0+2, lab=goodnames[p], srt=15, cex=1.2, col=1 )
    for( j in js0 ){
        x1  <- which( phens[sort.list(y[j,],dec=T)] == p )
        y1  <- y[j,p]
        lines( c(x0-6,x1), c(y0,y1), col=cols[j] )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}
legend( -6, 106, leg='A', bty='n', cex=1.9 )
legend( 'topright', fill=cols[js0], leg=rownames(y)[js0], bty='n', cex=1.3 )


###### Panel B
js_mtag	<- c( 4, 3, 7, 9, 8 )# 6=AllDepEnvs excluded bc almost identical to 9=Multitraits
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2  )
for( i in js_mtag )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )
legend( -6, 106, leg='B', bty='n', cex=1.9 )
legend( 'topright', fill=cols[js_mtag], leg=rownames(y)[js_mtag], bty='n', cex=1.3 )

base        <- matrix( y['LifetimeMDD',], nmeth, nphen, byrow=T )
diffs	    <- ( y - base ) / base * 100
sel         <- phens[ sort.list( diffs['MTAG.Envs',], dec=T )[1:5] ]
goodnames   <- sel
names(goodnames)  <- sel
goodnames[1:5]    <- c( 'College' ,'Alevels' ,'OtherQuals' ,'CurrSmoke' ,'MaternSmokeAtBirth' )

y.off   <- c( 0, 0, 5, -3, -15 )+30
x0s     <- c( 22, NA, NA, 39, 52 )
for( p in sel[c(1,4,5)] ){
    x0  <- x0s[which(sel==p)]
    y0  <- y['MTAG.Envs',p]+y.off[which(sel==p)]
    text( x0 + ifelse( p == sel[5], -1, -5 ), y0+6, lab=goodnames[p], srt=15, cex=1.2, col=1 )
    for( j in js_mtag[c(1,2,3,5)] ){
        x1  <- which( phens[sort.list(y[j,],dec=T)] == p )
        y1  <- y[j,p]
        lines( c(x0-6,x1), c(y0,y1), col=cols[j] )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

###### Panel C
js		<- c( 4, 3, 1, 2, 10, 11 )
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2  )
for( i in js )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )
legend( -6, 106, leg='C', bty='n', cex=1.9 )
legend( 'topright', fill=cols[js], leg=rownames(y)[js], bty='n', cex=1.3 )

y.off   <- c( 0, 0, 5, -3, -15 )+30
x0s     <- c( 22, NA, NA, 39, 52 )
for( p in sel[c(1,4,5)] ){
    x0  <- x0s[which(sel==p)]
    y0  <- y['MTAG.Envs',p]+y.off[which(sel==p)]
    text( x0 + ifelse( p == sel[5], -1, -5 ), y0+6, lab=goodnames[p], srt=15, cex=1.2, col=1 )
    #for( j in js[c(3,5)] ){
    for( j in js[c(4,6)] ){
        x1  <- which( phens[sort.list(y[j,],dec=T)] == p )
        y1  <- y[j,p]
        lines( c(x0-6,x1), c(y0,y1), col=cols[j] )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}
dev.off()
