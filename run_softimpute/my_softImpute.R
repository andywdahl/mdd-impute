library(softImpute)
my_softImpute <- function( Y, lambda.len=100, maxrank=TRUE, fixed.maxrank=ncol(Y)-1, nfolds=10, verbose=T, return_udv=FALSE ){

	## mus	<- colMeans( Y, na.rm=T )
	## Y		<- Y - matrix( mus, nrow(Y), ncol(Y), byrow=T ) 
	## sigs<- colMeans( Y^2, na.rm=T )
	## Y		<- Y / matrix( sqrt(sigs), nrow(Y), ncol(Y), byrow=T )

	obs	  <- which( ! is.na( Y ) )
	n_obs	  <- length(obs)
	folds	  <- split( obs, sample( rep( 1:nfolds, n_obs )[1:n_obs], n_obs ) )

	cv.loss <- matrix( NA, lambda.len, nfolds )

	lam0    <- lambda0(Y)
	lamseq  <- exp(seq(from=log(lam0+.2),to=log(.001),length=lambda.len))
	##### from Hastie blog: http://web.stanford.edu/~hastie/swData/softImpute/vignette.html
	#### lambda.len	<- 10
	#### lamseq			<- exp(seq(from=log(lam0+.2),to=log(1),length=lambda.len))

	for( f in 1:nfolds ){
		if( verbose ) print( f )

		Y.train	<- Y
		Y.train[ folds[[f]] ]	<- NA

		ranks     <- as.integer( lamseq )
		rank.max  <- ifelse( maxrank, 2, fixed.maxrank )
		warm      <- NULL
		for( i in seq(along=lamseq)){
			if( verbose ) cat( i, ' ' )
			out	    <- softImpute( x=Y.train, lambda=lamseq[i], rank=rank.max, warm=warm, maxit=1000 )

			Yhat	    <- complete( Y.train, out ) 
			cv.loss[i,f]  <- mean( (Yhat - Y)^2, na.rm=TRUE )

			warm      <- out
			if( maxrank ){
				ranks[i]  <- sum(round(out$d,4)>0)
				rank.max  <- min( ranks[i]+2, fixed.maxrank ) ### grows by at most 2, bounded by P/2
			}
		}
		if( verbose ) cat( '\n' )
	}
	loss		<- rowMeans( cv.loss )
	lambda	<- lamseq[ which.min(loss) ]

	out			<- softImpute( x=Y, rank.max = ncol(Y), lambda = lambda, maxit=1000 )
	Yhat		<- complete( Y, out ) 

	#Yhat		<- Yhat / matrix( sqrt(sigs), nrow(Y), ncol(Y), byrow=T )
	#Yhat		<- Yhat + matrix( mus       , nrow(Y), ncol(Y), byrow=T )
	if( return_udv ){
		u	<- out$u
		d	<- out$d
		v	<- out$v
	} else {
		u	<- NULL
		d	<- NULL
		v	<- NULL
	}

	return(list( 
		lambda  = lambda,
		lamseq	= lamseq,
		loss    = loss, 
		Y       = Yhat,
		u=u, d=d, v=v
	))

}
