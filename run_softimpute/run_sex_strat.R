rm( list=ls() )

it				<- as.numeric( commandArgs(TRUE)[[1]] )
k					<- as.numeric( commandArgs(TRUE)[[2]] )
sex				<- as.numeric( commandArgs(TRUE)[[3]] ) 

sinkfile	<- paste0( 'Rout/sexstrat_'	, it, '_', k, '_', sex, '.Rout' )
savefile	<- paste0( 'Rdata/sexstrat_'	, it, '_', k, '_', sex, '.Rdata' )
if( file.exists( savefile ) | file.exists( sinkfile ) ) stop() 
sink( sinkfile )

load('Rdata/setup_misc.Rdata')
load('Rdata/setup_data.Rdata')
source( 'my_softImpute.R' )

set.seed(1234+it)

if( sex == 0 ){
	sex.n	<- min(dat[,'sex'])
} else if( sex == 1 ){
	sex.n	<- max(dat[,'sex'])
}
nsex0	<- sum( dat[,'sex'] == min(dat[,'sex']) )
nsex1	<- sum( dat[,'sex'] == max(dat[,'sex']) )

print( dim (dat ) )
print( table(dat[,'sex'] ) ) 
dat		<- dat[ dat[,'sex'] == sex.n, ]
print( dim (dat ) )
dat	<- dat[ , -which( colnames(dat)=='sex' ) ]
dat	<- scale(dat)
print( dim (dat ) )
sel	 <- sel-1 ## bc sex is second column, first is age and not in sel

N		<- nrow(dat)
np	<- length(c(dat))

allerrs		<- array( NA, c( 3, nf, length(sel) ), dimnames=list( c( 'mse', 'cor', 'auc' ), as.character(fs), colnames(dat)[sel] )  )
for( f.i in 1:1 ){

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
	save( allerrs, nsex0, nsex1, file=savefile )
	rm( datit ); gc()
}
save( allerrs, nsex0, nsex1, file=savefile )
sink()
