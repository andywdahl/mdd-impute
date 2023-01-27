rm( list=ls() )

it				<- as.numeric( commandArgs(TRUE)[[1]] )
k					<- as.numeric( commandArgs(TRUE)[[2]] )
x					<- 1 

sinkfile	<- paste0( 'Rout/MTAG_variables_', it, '_', k, '_', x, '.Rout' )
savefile	<- paste0( 'Rdata/MTAG_variables_', it, '_', k, '_', x, '.Rdata' )
if( file.exists( savefile ) | file.exists( sinkfile ) ) stop()

sink( sinkfile )

load('Rdata/setup_misc.Rdata')
load('Rdata/setup_data.Rdata')
source( 'my_softImpute.R' )

set.seed(1234+it)

covars	<- c('age','sex', paste0( 'PC', 1:20 ) )
if( x == 1 ){
	phens	<- c(
		'LifetimeMDD',
		'GPpsy', 'Psypsy', 'DepAll', 'SelfRepDep', 'ICD10Dep',
		'townsend',
		'neuroticismscore.baseline', 'trauma',
		'stressbinary', 
		'mother.severedepression', 'father.severedepression', 'sibling.severedepression'
		, covars
	)
	#print( colnames(dat) )
	#print( phens[ ! phens %in% colnames(dat) ] )
	print( phens )
	print( dim(dat) )
	sel	<- names(sel)
	dat	<- dat[,phens]
	sel	<- which(colnames(dat) %in% phens)

	print( 'sel' ) 
	print( sel )
	print( dim(dat) )

} else {
	stop()
}

print( sel )
print( sel[ !sel %in% colnames(dat) ] )

N		<- nrow(dat)
np	<- length(c(dat)) 

allerrs		<- array( NA, c( 3, nf, length(sel) ), dimnames=list( c( 'mse', 'cor', 'auc' ), as.character(fs), colnames(dat)[sel] )  )
for( f.i in 1:2 ){

	print( f.i )
	f			<- fs[f.i]

	datit	<- dat
	if( k==1 ){	### copy-mask
		nobs	<- sum(!is.na(dat[,sel]))
		while( sum(!is.na(datit[,sel])) > ( nobs * (1-f) ) ){
			x1	<- sample(N,100)
			datit[x1,sel]	<- datit[x1,sel] + 0*datit[sample(N,100),sel]
		}
	} else if( k==2 ){	### completely at random
		datit[,sel]	<- datit[,sel] + sample( c(NA,0), prob=c(f,1-f), N*length(sel), rep=T )
	}

	allerrs[,f.i,]	<- objs( dat, my_softImpute( datit, lambda.len=100, nfolds=10 )$Y, sel )
	objs( dat, dat+rnorm(2), sel )
	save( allerrs, file=savefile )
	rm( datit ); gc()
}
save( allerrs, file=savefile )
sink()
