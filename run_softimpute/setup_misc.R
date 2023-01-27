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
	X			<- X		[,sel]
	Ximp	<- Ximp	[,sel]
	obs		<- which( X == Ximp )
	X[obs]<- Ximp[obs]	<- NA

	cuts	<- apply( X, 2, function(x) mean( range(x,na.rm=T) ) )

	rbind( 
		colMeans( (X-Ximp)^2, na.rm=T ) / apply( X, 2, function(x) var(x,na.rm=T) ),
		sapply( 1:ncol(X), function(p){
			a	<- NA
			try( a <- cor( X[,p], Ximp[,p], use='complete.obs' ) )
			a
		}),
		sapply( 1:ncol(X), function(p) mean(X[,p] == as.numeric( Ximp[,p] > cuts[p] ), na.rm=T ) ) ### this metric is never used in our paper, but asks: How accurate is the best possible binarization?
	)
} 

phens					<- colnames(dat)
main_phens		<- c( 'LifetimeMDD', 'SelfRepDep', 'DepAll', 'Psypsy', 'GPpsy', 'neuroticismscore.baseline' )
neurot_phens	<- phens[ grep( 'neuroticism'	, phens ) ] 
age10_phens		<- phens[ grep( 'atage10'			, phens ) ] 
icd10_phens		<- phens[ grep( 'ICD10'				, phens ) ] 
edu_phens			<- phens[ grep( 'educ'				, phens ) ] 

save.image('Rdata/setup_misc.Rdata')
