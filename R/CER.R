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

