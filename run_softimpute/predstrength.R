rm( list=ls() )
load('Rdata/setup_misc.Rdata')

it				<- as.numeric( commandArgs(TRUE)[[1]] )
if( it > 2 ) stop()

savefile	<- paste0( 'Rdata/crossvalidated_2fold_', it, '.Rdata' )
sinkfile	<- paste0( 'Rout/crossvalidated_2fold_'	, it, '.Rout' ) 
if( !file.exists(savefile) & !file.exists(sinkfile) ){
	sink( sinkfile ) 
	load('Rdata/setup_misc.Rdata')
	source( 'my_softImpute.R' )
	set.seed(1234) 

	print( dim( dat ) )
	discard_IIDs	<- c()
	for( j in ( 5*(it-1) + 1:5 ) )
		discard_IIDs	<- c( discard_IIDs, read.table( paste0( 'test_sets/test', j, '.txt' ) )[,1] )

	print( dim( dat ) )

	dat		<- dat [ -match(discard_IIDs,IIDs),]
	IIDs	<- IIDs[ -match(discard_IIDs,IIDs) ]
	FIDs	<- FIDs[ -match(discard_IIDs,FIDs) ]

	print( dim( dat ) )
	print( dat[1:10,1:4] )

	out	<- my_softImpute( dat, lambda.len=100, nfolds=10, return_udv=TRUE )
	save( IIDs, FIDs, out, file=savefile ) 
	rm( out, dat ) 
	sink()
}
