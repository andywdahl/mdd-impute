rm( list=ls() ) 
if( file.exists( 'Rdata/setup_misc.Rdata' ) ) stop()
load('Rdata/setup_data.Rdata')

sink( 'Rout/setup_misc.Rout' )

np	<- nrow(dat)*ncol(dat)
nf	<- 10

sel	<- which( colMeans( is.na(dat) ) > .01 )

fs	<- c(.01,seq( .05, .95, len=nf-1 ))

colnames(dat)[sel]
colnames(dat)[-sel]
print( round( colMeans( is.na(dat) )[sel]  , 3 ))
print( round( colMeans( is.na(dat) )[-sel] , 3 ))

objs	<- function( X, Ximp, sel=1:ncol(X) ){

    ## subset to interesting columns in sel
	X		<- X	[,sel]
	Ximp	<- Ximp	[,sel]
	obs		<- which( X == Ximp )
	X[obs]  <- Ximp[obs]	<- NA

	cuts	<- apply( X, 2, function(x) mean( range(x,na.rm=T) ) )  ## potential thresholds for transforming quantitative imputations to {0,1}

	rbind( 
		colMeans( (X-Ximp)^2, na.rm=T ) / apply( X, 2, function(x) var(x,na.rm=T) ),    # MSE
		sapply( 1:ncol(X), function(p){
			a	<- NA
			try( a <- cor( X[,p], Ximp[,p], use='complete.obs' ) )
			a
		}), ## basically correlation, but with 'try' catches situations with all missing
		sapply( 1:ncol(X), function(p) mean(X[,p] == as.numeric( Ximp[,p] > cuts[p] ), na.rm=T ) )
	)
} 

save.image('Rdata/setup_misc.Rdata')
