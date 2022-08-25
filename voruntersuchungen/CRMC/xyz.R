#frist to sercond

if (!require( "statip", character.only = TRUE)) {
  install.packages("statip", dependencies = TRUE)     
}
if (!require( "provenance", character.only = TRUE)) {
    install.packages('provenance', dependencies = TRUE)
}

library(data.table)
library(transport)
library(provenance)

makejenshan <- function( v1, v2 ){
    v1[is.na(v1)] <- 0.0
  #print(v1)
v2[is.na(v2)] <- 0.0
  #sapply(x, class) 
  #print(v2)

    ll <- length(v1)
    sdiver = 0.0
    rdiver = 0.0
    m <- rep(0.0, ll)

    for( n in 1:ll ){
        m[n] = (v1[n]+v2[n])/2
    }
    
    for( n in 1:ll ){
        #message(sdiver, "  ",v1[n]*log10(v1[n]/m[n]), " v1n ",v1[n], " log ", log10(v1[n]/m[n]), " v d m ", v1[n]/m[n], " mn ",m[n], "\n" )
        
        if(v1[n] != 0){
            sdiver = sdiver+(v1[n]*log10(v1[n]/m[n]))/2
        }
        
        if(v2[n] != 0){ 
            rdiver = rdiver+(v2[n]*log10(v2[n]/m[n]))/2
        }
        
    }
    return( sdiver+rdiver )
}

hellinger <- function( x, y, lower = -Inf, upper = Inf, method = 1, ...){
    x[is.na(x)] <- 0.0
  print(x)
y[is.na(y)] <- 0.0
  #sapply(x, class) 
  print(y)
  fx <- densityfun(x, ...)
  fy <- densityfun(y, ...)
  if (method == 1) {
    g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2
    h2 <- stats::integrate(g, lower, upper, rel.tol=.Machine$double.eps^.05)$value/2
  } else if (method == 2) {
    g <- function(z) (fx(z)*fy(z))^0.5
    h2 <- 1 - stats::integrate(g, lower, upper, rel.tol=.Machine$double.eps^.05)$value
  } else {
    stop("incorrect 'method' argument", call. = FALSE)
  }
  sqrt(h2)
}

calchellingerdist <- function( da, clos ){
    rn <- as.matrix(da[1:nrow(da),1]) # select data of first column as names
    
    le <- nrow(rn) # get length of it to know how much itter and output size
    lo <- length(clos)
    
    print(le)
    print(length(rn))
    yy <- matrix(1:(le*le), nrow = le, dimnames = list(rn))

    #print(yy)
    for( n in 2:le ){
        for( m in 2:le ){
                print(n)
                print(m)
                #print(da[n,clos, with=FALSE])
                #print(as.numeric(as.matrix(da[n,clos, with=FALSE])))
                if( m != n ){
                   yy[n,m] = makejenshan( as.numeric(as.matrix(da[n,clos, with=FALSE])), as.numeric(as.matrix(da[m, clos, with=FALSE])))#hellinger( as.numeric(as.matrix(da[n,clos, with=FALSE])), as.numeric(as.matrix(da[m, clos, with=FALSE])), -Inf, Inf )#1-hellinger(x[n,], x[m,], -Inf, Inf)
                } else if( m < n ){ 
                   yy[n,m] = yy[n,m]
                } else {
                   yy[n,m] = 1.0 # unequal to self as a definition
                }
        }
    }
    return(yy)
}

wasser <- 
function( da, clos, exp = 1.5, scale = TRUE ){
    rn <- as.matrix(da[1:nrow(da),1]) # select data of first column as names
    
    le <- nrow(rn) # get length of it to know how much itter and output size
    lo <- length(clos)
    
    print(le)
    print(length(rn))
    y <- matrix(1:(le*le), nrow = le, dimnames = list(rn))

    for( n in 1:le ){
        for( m in 1:le ){
                if(m != n){
                   y[n,m] = wasserstein1d( as.numeric(as.matrix(da[n,clos, with=FALSE])), as.numeric(as.matrix(da[m, clos, with=FALSE])), p = 1, wa = NULL, wb = NULL)#1-wasserstein1d( x[n,], x[m,], p = 1, wa = NULL, wb = NULL)
                } else if( m < n ){ 
                   y[n,m] = y[n,m]
                } else {
                   y[n,m] = 1.0 # unequal to self as a definition
                }
        }
    }
    
    return(y)
    
    
}

memetadata <- fread("TKAMetadaten.csv", key="Lab.-No.")

xrfdata <- fread("TKAXRF.csv", key="Lab.-No.")

xrfdata <- fread("TKAXRF.csv", key="Lab.-No.")
mirdata <- merge( xrfdata, memetadata )
merdata<- mirdata[ware=="Terra Nigra" | ware=="Terra Sigillata", ]

#merdata[is.na(merdata)] <- 0.0 #remove NA fealds by setting them zero

print(nrow(merdata))
#columnsfordistcomp <- c("SiO2", "TiO2", "Al2O3", "Fe2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "V", "Cr", "Ni", "Cu", "Zn", "Rb", "Sr", "Y", "Zr", "Nb", "Ba", "La", "Ce", "Pb", "Th")
#columnsfordistcomp <- c("SiO2","TiO2", "Al2O3", "Fe2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5") #alle Oxide
#columnsfordistcomp <- c("TiO2", "Al2O3", "Fe2O3", "MnO", "MgO") #metalloxide
#columnsfordistcomp <- c("SiO2", "CaO", "Na2O", "K2O", "P2O5") #nm und halbm oxide
#columnsfordistcomp <- c("SiO2", "CaO", "Na2O", "K2O", "MgO", "Fe2O3", "Al2O3" ) #selection from paper
columnsfordistcomp <- c("SiO2", "CaO", "Na2O", "MgO", "Al2O3" )
#columnsfordistcomp <- c("V", "Cr", "Ni", "Cu", "Zn", "Rb", "Sr", "Y", "Zr", "Nb", "Ba", "La", "Ce", "Pb", "Th") #spurenelemente - weder canberra not euklide gut

#containning the abs messured values
dd <- as.matrix(merdata[ ,columnsfordistcomp, with=FALSE])

#dd <- CLR(bb) #change values to lagratios
print(nrow(dd))




distancematrix <- dist(dd, method = "manhattan") 
#distancematrix <- wasser( merdata, columnsfordistcomp )
print(nrow(distancematrix))
#print( is.na(distancematrix))
distancematrix[is.na(distancematrix)] <- 100000.0
#print( is.na(distancematrix))

labe <- as.matrix(merdata[ ,"ware", with=FALSE])[,1]
#print(labe)
hc <- hclust(distancematrix, method = "ward.D" )
hc$labels <- labe
png("plotdendogram.png",width=10000,height=800)

par(cex=1,font=5)
plot(hc, hang = -1, ylab = "Distance", main = "Dendrogramm der chem. Zusammensetzung, mit Warenarten gelabelt.")
dev.off()
