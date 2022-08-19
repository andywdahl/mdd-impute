rm( list=ls() )

x		<- read.table( 'data/allmaplesyrup.txt', head=T ) #columns: pheno,prs,meanratio,sdratio,meanr2,sdr2,meanmddr2,sdmddr2,lower,upper
x		<- x[sort.list(as.character(x[,1])),]
nmeth	<- 11
nphen   <- nrow(x)/nmeth
phens   <- as.character(x[1+nmeth*(1:nphen-1),'pheno'])
meths	<- as.character(x[1:nmeth,2])

## check name/pheno matching worked
stopifnot( all.equal( rep(phens,each=nmeth), as.character(x[,'pheno']) ) )
stopifnot( all.equal( rep(meths,nphen), as.character(x[,'prs']) ) )

#add individual runs of extra PRS
locread <- function(file,prs){
    x       <- read.table( file, head=T )
    pheno   <- x[,1]
    prs     <- rep( prs, nrow(x) )
    x       <- cbind( pheno, prs, x[,-1] )
    if( prs[1] == '23andMe' )
        x   <- x[-which( pheno == 'neuroticismscore' ),] # one extra run in file--innocuous
    x
}

x		<- rbind( x, locread( 'data/23andMe.txt', '23andMe' ), locread( 'data/PGC29.txt', 'PGC29' ) )
meths   <- c( meths, '23andMe', 'PGC29' )#, 'ipsych.txt' )
nmeth	<- nmeth+2

## check name/pheno matching worked
x	<- x[sort.list(as.character(x[,1])),]
stopifnot( all.equal( rep(meths,nphen), as.character(x[,'prs']) ) )
stopifnot( all.equal( rep(phens,each=nmeth), as.character(x[,'pheno']) ) )

#### step 1: reshapes the data into #methods by #phenotypes matrix
y       <- matrix( x[,'meanratio'], nmeth, nphen )
y.sd    <- matrix( x[,'sdratio'], nmeth, nphen )

rownames(y)     <- rownames(y.sd)  <- meths
colnames(y)     <- colnames(y.sd)  <- phens

y2      <- matrix( x[,'meanr2'], nmeth, nphen )
y2.sd   <- matrix( x[,'sdr2'], nmeth, nphen )
pvs     <- apply( y2/y2.sd, 1:2, function(z) pnorm(-z) ) ### test if the r2s are nonzero
rownames(pvs) <- meths
colnames(pvs) <- phens

save( y, y.sd, pvs, phens, nphen, meths, nmeth, file='data/allmaplesyrup_reshaped.Rdata' )
