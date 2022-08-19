rm( list=ls() )
load( file='data/allmaplesyrup_qcd.Rdata' )
rownames(y)[9] <- 'MTAG.All'

###### subset to phenotypes that are significantly predicted by any PRS
minps   <- apply( pvs, 2, min )
sub     <- which( minps < .05/ncol(pvs) )
y       <- y[,sub]
nphen   <- length(sub)
phens   <- phens[sub]

base    <- matrix( y['LifetimeMDD',], nmeth, nphen, byrow=T )
diffs	<- ( y - base ) / base * 100

pdf( 'figs/ExtendedDataFig4.pdf', width=12+2.5, height=4, useDingbats=FALSE )
layout( matrix(1:4,1,4), width=c(4,4,4,2.2) )
par( mar=c(4.8,4.9,.5,.5) )

plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65 )
js0 <- c(4,3,12)
for( i in js0 )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

sel         <- phens[ sort.list( y['GPpsy',], dec=T )[1:5] ]
goodnames   <- sel
names(goodnames)   <- sel
goodnames[4]    <- 'Neurot'
x0s     <- c( 5, 48, 20, 15, 30 )
y.off   <- c( 90, 50, -15, 55, 39 )
for( p in sel ){
    x0  <- x0s[which(sel==p)]
    y0  <- diffs['GPpsy',p]+y.off[which(sel==p)]+30
    text( x0 + ifelse( p == 'SelfRepDep', 4, 0 ), y0+7, lab=goodnames[p], srt=15, cex=1.5, col=1 )
    for( j in js0[-1] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

js_mtag	<- c( 4, 3, 7, 9, 8 )
plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65  )
for( i in js_mtag )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

sel         <- phens[ sort.list( diffs['MTAG.Envs',], dec=T )[1:5] ]
goodnames   <- sel
names(goodnames)    <- sel
goodnames[1:5]      <- c( 'College' ,'Alevels' ,'OtherQuals' ,'CurrSmoke' ,'MaternSmokeAtBirth' )
x0s     <- c( 5, NA, NA, 23, 44 )
y.off   <- c( -20, NA, NA, 55, 9 )
for( p in sel[c(1,4,5)] ){
    x0  <- x0s[which(sel==p)]
    y0  <- diffs['MTAG.Envs',p] +y.off[which(sel==p)]
    text( x0 + ifelse( p == sel[1], 6, 0 ), y0+ifelse( p == sel[5], 28, 11 ), lab=goodnames[p], srt=15, cex=1.5, col=1 )
    for( j in js_mtag[c(3:5)] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65  )
js		<- c( 4, 3, 2, 1, 11, 10 )
for( i in js )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

x0s     <- c( 4, NA, NA, 23, 44 )
for( p in sel[c(1,4,5)] ){
    x0  <- x0s[which(sel==p)]
    y0  <- diffs['MTAG.Envs',p] +y.off[which(sel==p)]
    text( x0 + ifelse( p == sel[1], 2, 0 ), y0+ifelse( p == sel[5], 22, 8 ), lab=goodnames[p], srt=15, cex=1.5, col=1 )
    for( j in js[c(4,6)] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

alljs  <- unique(c( js0, js_mtag, js ))
par( mar=c(0,0,0,0) )
plot.new()
legend( 'center', fill=cols[alljs], leg=rownames(y)[alljs], bty='n', cex=1.6 )

dev.off()
