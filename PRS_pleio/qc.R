rm( list=ls() )
load( 'data/allmaplesyrup_reshaped.Rdata' )

#### remove centre
y.centre<- y	[, grep( 'centre',phens )]
phens	<- phens[ -grep( 'centre',phens )]

js		<- c( 4, 3, 1, 2, 10, 11 )
cols	<- c( 'orange3', 'orange1', 'grey', 1, 'blue1', 6, 'blue4', 'cornflowerblue', 'cyan', 'purple3', 'purple1', 'pink1', 'pink4' )

#### remove PCs
phens	    <- setdiff( phens, c( 'array', paste0( 'PC', 1:20 ) ) )
y			<- y[,phens]
y.sd		<- y.sd[,phens]
nphen		<- length(phens)

### remove livew.grandchild -- has huge std error, 3x the next largest
tail( round( sort( colMeans( y.sd ) ), 2 ) )
phens   <- setdiff( phens, c( 'livew.grandchild.baseline', 'LifetimeMDD' ) )
y	    <- y    [,phens]
y.sd    <- y.sd [,phens]
pvs     <- pvs  [,phens]
nphen   <- length(phens)

save( y, y.sd, pvs, phens, nphen, meths, nmeth, js, cols, file='data/allmaplesyrup_qcd.Rdata' )
