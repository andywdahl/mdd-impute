rm( list=ls() )

it				<- as.numeric( commandArgs(TRUE)[[1]] )
sinkfile	<- paste0( 'Rout/saveimp_' , it, '.Rout' )
savefile	<- paste0( 'Rdata/saveimp_', it, '.Rdata') 
if( file.exists( savefile ) | file.exists( sinkfile ) ) stop() 
sink( sinkfile ) 

load('Rdata/setup_misc.Rdata')
load('Rdata/setup_data.Rdata')
source( 'my_softImpute.R' )

set.seed(1234+it)

N	<- nrow(dat)
np	<- length(c(dat))

f	<- fs[1]

datit	<- dat
nobs	<- sum(!is.na(dat[,sel]))
while( sum(!is.na(datit[,sel])) > ( nobs * (1-f) ) ){
	x1	<- sample(N,100)
	datit[x1,sel]	<- datit[x1,sel] + 0*datit[sample(N,100),sel]
} 
imp	<- my_softImpute( datit, lambda.len=100, nfolds=10 ) 
save( imp, file=savefile )

sink()
