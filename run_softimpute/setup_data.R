rm( list=ls() )
if( file.exists( 'Rdata/setup_data.Rdata' ) ) stop()
sink( 'Rout/setup_data.Rout' )

dat	<- read.table( 'data/mddimpute.phenotypes.final.txt', header=T )
dim( dat )
colnames( dat )

dat[1:4,1:4]
IIDs	<- dat[,'IID']
FIDs	<- dat[,'FID']
dat		<- dat[,-(1:2)]
dat[1:4,1:4]

dim( dat )
head( dat )
table( rowSums( is.na( dat ) ) )

N		<- nrow(dat)
P		<- ncol(dat)

print( sort( colMeans( is.na( dat ) ) ) )
dat	<- scale(dat)
badcols	<- which( colMeans( is.na( dat ) ) > .99 )
stopifnot( length(badcols) == 0 )
dim(dat)

save.image('Rdata/setup_data.Rdata')
sink()
