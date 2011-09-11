# Revised Silhouette

revisedsil<-function(
			d,W,C,out=NULL,CASEofINT=out,
			col1="black",CASEofINT2=NULL,
			col2="red",print.plot=TRUE){
        
        second.min<-function(x){
    	s<-sort(x,decreasing=FALSE)
    	return(s[2])}
    	c.size<-1.5
	    # out will be excluded for the calculation of cluster center
        trans.d<-sweep(d,2,sqrt(W),FUN="*")[,W!=0]
        K<-max(C);sumW<-sum(W)
        N<-nrow(d)
        if (is.null(out)) out<-N+1
            w.mu<-matrix(NA,K,sum(W!=0))
            C2<-C;C2[out]<--1;C2<-C2[1:N]
        for (k in 1 : K) 
             {
             w.mu[k,]<-colMeans(trans.d[C2==k,,drop=FALSE],na.rm=TRUE)
             }
        WdisC<-WDISC(trans.d,w.mu,K,N,W[W!=0],sumW)
              
        ai<-apply(WdisC,1,min)
        bi<-apply(WdisC,1,second.min)

        sil<-(bi-ai)/bi
        I<-1:N;obs.i<-sil.i<-list()
        for ( k in 1 : K){
         	sil.k<-sil[C==k]
         	sil.k.order<-sort(sil.k,decreasing=TRUE);
         	sil.rank<-rank(sil.k) # vector of rank 1=smallest 
         	I.k<-I[C==k]
         	I.k.order<-NULL
         	for ( i in 1 : length(I.k)) for ( j in 1 : length(I.k)) 
         	if(sil.rank[j]==i) I.k.order[length(I.k)-i+1]<-I.k[j]
                # sil ranked in decreasing order within k^th cluster
         	sil.i[[k]]<-sil.k.order 
         	# corresponding observation index
         	obs.i[[k]]<-I.k.order
         	}
    	
    	obs.index<-obs.i[[1]]
    	sil.index<-sil.i[[1]]
    	for ( k in 2 : K) {
    		obs.index<-c(obs.index,obs.i[[k]])
    		sil.index<-c(sil.index,sil.i[[k]])
    		}
    	
    	InterestOrNot<-rep("gray70",N)
    	for (o in 1 : length(CASEofINT))
    	InterestOrNot[which(obs.index==CASEofINT[o])]<-col1
    	
    	if(print.plot){
    		plot(sil.index,type="h",axes=FALSE,col=InterestOrNot,lwd=15,
                 xlab="",ylab="rev.silhouette",cex.main=1.5,cex.lab=1.5,
                 ylim=c(0,1))
            for ( i in 1 : length(CASEofINT2))
                text((1:N)[obs.index==CASEofINT2[i]],0.8,labels=paste(CASEofINT2)[i],col=col2)     
            axis(1,at=c(-10,100),cex=c.size)
            axis(1,at=(1:N)[InterestOrNot==col1],
                 labels=obs.index[InterestOrNot==col1]
                 ,cex=c.size)
            #axis(1,at=1:N,labels=obs.index,cex=0.5)     
            axis(2,seq(0,1,0.25),seq(0,1,0.25),cex=c.size)
            
           Nobsothers<-0
           for (k in 1 : K){
            	nk<-length(sil.i[[k]])
      	        sil.k.mean<-round(mean(sil.i[[k]]),2)
      	        text(Nobsothers+nk/2,0.5,sil.k.mean,cex=2.5)
      	        Nobsothers<-Nobsothers+nk
      	        }
      	}
    return(list(trans.mu=w.mu,WdisC=WdisC,sil.order=sil,sil.i=sil.i,obs.i=obs.i))
	}


	
