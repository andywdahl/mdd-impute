rm( list=ls() )
load( file='data/allmaplesyrup_qcd.Rdata' ) #y, y1, y1.sd, phens, nphen, meths, nmeth, 
rownames(y)[9] <- 'MTAG.All'

###### subset to phenotypes that are significantly predicted by any PRS
minps   <- apply( pvs, 2, min )
sub     <- which( minps < .05/ncol(pvs) )
y       <- y[,sub] * 100
nphen   <- length(sub)
phens   <- phens[sub]

###### PRS sets
js0     <- c(4,3,12)            # LifetimeMDD, GPpsy, 23andMe
js_mtag	<- c( 4, 3, 7, 9, 8 )   # 6=AllDepEnvs excluded bc almost identical to 9=Multitraits
js_imp	<- c( 4, 3, 1, 2, 10, 11 )

pdf( 'figs/Fig6.pdf', width=12, height=4, useDingbats=FALSE )
par( mfrow=c(1,3) )
par( mar=c(4.2,4.2,.7,.7) )

## panel A: highlight top GPpsy phens
selA        <- phens[ sort.list( y['GPpsy',], dec=T )[1:5] ]
labsA       <- selA
names(labsA)<- selA
labsA[4]    <- 'Neurot'

## panel B, C: highlight most inflated phens by MTAG.Envs
base        <- matrix( y['LifetimeMDD',], nmeth, nphen, byrow=T )
diffs	    <- ( y - base ) / base * 100
selBC       <- phens[ sort.list( diffs['MTAG.Envs',], dec=T )[1:5] ]
labsBC      <- selBC
names(labsBC)<- selBC
labsBC[1:5]    <- c( 'College','Alevels','OtherQuals','CurrSmoke','MaternSmokeAtBirth' )

###### Panel A
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2 )
for( i in js0 )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )

for( p in selA ){
    x0  <- c( 10, 15, 20, 20, 30 )[which(selA==p)]+10
    y0  <- y['GPpsy',p]+c( 0, 0, 5, -3, 0 )[which(selA==p)]
    text( x0 + ifelse( p == 'SelfRepDep', 2, 0 ), y0+2, lab=labsA[p], srt=15, cex=1.2, col=1 )
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
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2  )
for( i in js_mtag )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )
legend( -6, 106, leg='B', bty='n', cex=1.9 )
legend( 'topright', fill=cols[js_mtag], leg=rownames(y)[js_mtag], bty='n', cex=1.3 )

for( p in selBC[c(1,4,5)] ){
    x0  <- c( 22, NA, NA, 39, 52 )[which(selBC==p)]
    y0  <- y['MTAG.Envs',p]+c( 0, 0, 5, -3, -15 )[which(selBC==p)]+30
    text( x0 + ifelse( p == selBC[5], -1, -5 ), y0+6, lab=labsBC[p], srt=15, cex=1.2, col=1 )
    for( j in js_mtag[c(1,2,3,5)] ){
        x1  <- which( phens[sort.list(y[j,],dec=T)] == p )
        y1  <- y[j,p]
        lines( c(x0-6,x1), c(y0,y1), col=cols[j] )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

###### Panel C
plot( c(1,nphen), c(0,100), type='n', ylab='PRS Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.2  )
for( i in js_imp )
    lines   ( 1:nphen, sort(y[i,],dec=T), col=cols[i], lwd=2.5 )
legend( -6, 106, leg='C', bty='n', cex=1.9 )
legend( 'topright', fill=cols[js_imp], leg=rownames(y)[js_imp], bty='n', cex=1.3 )

for( p in selBC[c(1,4,5)] ){
    x0  <- c( 22, NA, NA, 39, 52 )[which(selBC==p)]
    y0  <- y['MTAG.Envs',p]+c( 0, 0, 5, -3, -15 )[which(selBC==p)] + 30
    text( x0 + ifelse( p == selBC[5], -1, -5 ), y0+6, lab=labsBC[p], srt=15, cex=1.2, col=1 )
    for( j in js_imp[c(4,6)] ){
        x1  <- which( phens[sort.list(y[j,],dec=T)] == p )
        y1  <- y[j,p]
        lines( c(x0-6,x1), c(y0,y1), col=cols[j] )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}
dev.off()


###### Extended Data version plots the change relative to observed LifetimeMDD
pdf( 'figs/ExtendedDataFig4.pdf', width=12+2.5, height=4, useDingbats=FALSE )
layout( matrix(1:4,1,4), width=c(4,4,4,2.2) )
par( mar=c(4.8,4.9,.5,.5) )

plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65 )
for( i in js0 )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

for( p in selA ){
    x0  <- c( 5, 48, 20, 15, 30 )[which(selA==p)]
    y0  <- diffs['GPpsy',p]+c( 90, 50, -15, 55, 39 )[which(selA==p)]+30
    text( x0 + ifelse( p == 'SelfRepDep', 4, 0 ), y0+7, lab=labsA[p], srt=15, cex=1.5, col=1 )
    for( j in js0[-1] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65  )
for( i in js_mtag )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

for( p in selBC[c(1,4,5)] ){
    x0  <- c( 5, NA, NA, 23, 44 )[which(selBC==p)]
    y0  <- diffs['MTAG.Envs',p] +c( -20, NA, NA, 55, 9 )[which(selBC==p)]
    text( x0 + ifelse( p == selBC[1], 6, 0 ), y0+ifelse( p == selBC[5], 28, 11 ), lab=labsBC[p], srt=15, cex=1.5, col=1 )
    for( j in js_mtag[c(3:5)] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

plot( c(1,nphen), c(-80,340), type='n', ylab='Excess Pleiotropy (%)', xlab='Secondary Phenotypes', cex.lab=1.65  )
for( i in js_imp )
    lines   ( 1:nphen, sort(diffs[i,],dec=T), col=cols[i], lwd=4 )

for( p in selBC[c(1,4,5)] ){
    x0  <- c( 4, NA, NA, 23, 44 )[which(selBC==p)]
    y0  <- diffs['MTAG.Envs',p] +c( -20, NA, NA, 55, 9 )[which(selBC==p)]
    text( x0 + ifelse( p == selBC[1], 2, 0 ), y0+ifelse( p == selBC[5], 22, 8 ), lab=labsBC[p], srt=15, cex=1.5, col=1 )
    for( j in js_imp[c(4,6)] ){
        x1  <- which( phens[sort.list(diffs[j,],dec=T)] == p )
        y1  <- diffs[j,p]
        lines( c(x0,x1), c(y0-10,y1), col=cols[j], lwd=1.5 )
        points( x1, y1, col=cols[j], pch=16, cex=1.7 )
    }
}

alljs  <- unique(c( js0, js_mtag, js_imp ))
par( mar=c(0,0,0,0) )
plot.new()
legend( 'center', fill=cols[alljs], leg=rownames(y)[alljs], bty='n', cex=1.6 )

dev.off()
