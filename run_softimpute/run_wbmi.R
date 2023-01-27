rm( list=ls() )

it				<- as.numeric( commandArgs(TRUE)[[1]] )
k					<- as.numeric( commandArgs(TRUE)[[2]] )

sinkfile	<- paste0( 'Rout/wbmi_'	, it, '_', k, '.Rout' )
savefile	<- paste0( 'Rdata/wbmi_', it, '_', k, '.Rdata' )
if( file.exists( savefile ) | file.exists( sinkfile ) ) stop()

sink( sinkfile )

load('Rdata/setup_misc.Rdata')
load('Rdata/setup_data.Rdata')
source( 'my_softImpute.R' )

set.seed(1234+it)

N	<- nrow(dat)
np	<- length(c(dat))


#load bmi and match IDs
bmi	<- read.table( 'data/bmi.phen' ) 
ids	<- read.table( 'data/mddimpute.phenotypes.final.txt', header=TRUE, colClasses = c(rep("integer", 2), rep("NULL", ncol(dat))) )
#ids0	<- read.table( 'data/mddimpute.phenotypes.final.txt', header=TRUE )[,1:2]
#all.equal( ids, ids0 )

bmi	<- bmi[match( ids[,1], bmi[,1] ),]
stopifnot( all.equal( bmi[,1], ids[,1] ) )

bmi	<- scale( bmi[,3] )
dat	<- cbind( dat, bmi )

summary( lm( bmi ~ dat[,'age'] ) )
summary( lm( bmi ~ dat[,'sex'] ) )

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
	save( allerrs, file=savefile )
	rm( datit ); gc()
}
save( allerrs, file=savefile )
sink()
