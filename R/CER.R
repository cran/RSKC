CER <-
function(ind, true.ind,nob=length(ind)){
	return(sum(abs(one(true.ind)-one(ind)))/choose(nob,2))}

one <-
function(index){
	on<-NULL 
	c<-combn(index,2)
	c<-t(c)
	on<-1*(c[,1]==c[,2])
	return(on)
	}

norm1 <-
function(y){sum(abs(y))}

norm2 <-
function(y){sqrt(sum(y^2))}


## function specific for the opt digits
## generate bitmap of given observation
showbitmap <-function(index)
  {
    ## data(bitmapMat) ## lazyloading
    ## data(bitmapLab) ## lazyloading
    Nbit=32
    for (iindex in 1 : length(index))
      {
        indivindex<- index[iindex]
        for (ibit in 1 : Nbit)
          {
            cat("\n")
            cat(bitmapMat[[indivindex]][ibit])
          }
        cat("\n obs=",indivindex," true digit=",bitmapLab[indivindex]," \n")
      }
  }


## function declaration

showDigit <- function(index,cex.main=1)
  {
    ## data(DutchUtility) ## lazyloading
    ## 4. DutchUtility-pix: 240 pixel averages in 2 x 3 windows; 
    ## 16 by 15
    ncols = 15 ## replace to 15 
    nrows = 16 ## replace to 16
    labels <- rep(0:9,each=200)
    plot(NA,xlim=c(0,ncols),ylim=c(0,nrows),axes=FALSE,xlab="",ylab="",cex.main=cex.main,
         main=paste("observation",index," True digit",labels[index],sep=""))
    abline(h=0:ncols,v=0:nrows,lty=2,col="gray70")
    axis(1,(1:ncols)-0.5,1:ncols,lty=0,cex=0.5)
    axis(2,(1:nrows)-0.5,nrows:1,padj=1,lty=0,cex=0.5)
    ##
    cols <- gray.colors(n=6,start=0.9,end=0.3)
    for (iy in 1 : nrows)
      {
        for (ix in 1 : ncols)
          {
            ## each row vector of a matrix DutchUtility contains the bitmap X; 
            ## DutchUtility[1,] = c(X[1,],X[2,],...) 
            Pickedcolor <- cols[DutchUtility[index,ix+ncols*(iy-1)]]
            polygon(x=c(ix-1,ix-1,ix,ix),
                    y=c(nrows-iy,nrows-iy+1,nrows-iy+1,nrows-iy),
                    col=Pickedcolor,
                    border=FALSE)
          }
      }
  }



Sensitivity <- function(L,trueL)
  { 
   unitrueL<- unique(trueL)
   K <- length(unitrueL)
   senst <- correspondClass <- rep(NA,K)
   tbl <- table(L,trueL)
   ## The true number of observations in each clusters
   trueTot<- colSums(tbl)
   prMat <- scale(tbl,scale=trueTot,center=FALSE)
   senst <- apply(prMat,2,max)*100
   correspondClass <-  apply(prMat,2,which.max)

   re <- rbind(senst,correspondClass)
   colnames(re) <- unitrueL
   rownames(re) <- c("Sensitivity (%)","Class label by alg.")
   return(re)
 }
