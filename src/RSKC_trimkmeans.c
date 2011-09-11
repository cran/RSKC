

# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <R_ext/Linpack.h>

// RSKC_trimkmeans




void sample(const int k, int size, int rN[])
// This function act as R function sample(0:(size-1),k,replace=TRUE)
// k; the number of integers to be sampled between 0 and size - 1 WITH replacement;
// each integer has equal chances to be selected
// rN[k]; contains randomly selected k integers 
{
  int temp=0; 
  for (int i = 0; i < k; i++ )
    {
      rN[i] = (int) round(unif_rand()*100000) % size;
      
      //Rprintf("\v inside rN[i] %d",rN[i]);
      size--; // reduce the size by one

      if (i > 0)
	{
	  temp = rN[i]; 
	  for (int ii = 0 ; ii < i; ii ++)
	    {
	      if (rN[ii] <= rN[i])
		temp++;
	    }
	  rN[i] = temp;
	}
    }
  // e.g. select 4 integers from 0,1,...,9 (size = 10, n = 3)
  // for i=0, 8 was selected out of 0 to 9; rN[0] = 8; then size is reduced to 8
  // for i=1, we want to random sample from 0,1,..,7,9
  //          2 was selected out of 0 to 8; rN[1] = 2; then size is reduced to 7; 
  //          rN[ii] > rN[1] for all ii < 1, rN[1] does not change the value; rN[1] = 2
  // for i=2, we want to random sample from 0,1,3,..,7,9
  //          5 was selected out of 0 to 7; rN[2] = 5; then size is reduced to 5; 
  //          rN[ii] <= rN[2] is true when ii = 0 so rN[2] = rN[2] + 1; rN[2] = 6;
  // for i=3, we want to random sample from 0,1,3,4,5,7,9
  //          5 was selected out of 0 to 6; rN[3]=5; then size is reduced to 4;
  //          rN[ii] <= rN[3] is true when ii=1 and ii = 2 so rN[3] = rN[3] + 2; rN[3] = 7 
  // Final output should be: rN = {8,2,6,void}
} 


void realmin(const double iWdisC[],const int k, int *icls, double*themin) 
{
  // outputs 
  // icls: the position of the vector contiWdisCing the min
  // the min: the min value
  
  *icls=0;
  *themin=iWdisC[0];
  
  if (k == 1)
    {
      *themin=iWdisC[1];
      
    }else{
    
    for (int ik=1; ik < k;ik++ )
      {
	if (*themin > iWdisC[ik])
	  {
	    *themin=iWdisC[ik];
	    *icls=ik;
	  }
      }
  }
}

int identical(const int oldclass[],const int iclass[],const int n)
{
  int temp = 1; // two vectors are identical  
  for (int in=0;in < n;in++)
    {
      if (oldclass[in]!=iclass[in])
	{
	  temp = 0; // find one entry that oldclass[i]!= iclass
	  // the two vectors are NOT identical
	  break;
	}
    }
  return temp;
}


SEXP RSKC_trimkmeans(SEXP data_,
		     SEXP n_,
		     SEXP p_,
		     SEXP k_,
		     SEXP nin_,
		     //SEXP scaling,
		     SEXP runs_,
		     SEXP points_,
		     //SEXP countmode,
		     SEXP maxit_
		     )
{
  // the code of trimkmeans is modified so that
  // for alpha = 0, no case is trimmed. i.e. k-means clustering is performed     

    
  // === Importing values === // 
  const double *data=REAL(data_); // a vector of length n*p, containing the observations to be clustered 
  // data = c(Rdata[,1], Rdata[,2],...,Rdata[,ncol(Rdata)]); data[in+N*ip] = Rdata[in,ip]
  const int n = INTEGER(n_)[0];
  const int p = INTEGER(p_)[0]; 
  const int k = INTEGER(k_)[0];
  const int nin = INTEGER(nin_)[0];  // The number of trimmed observations
  const int runs = INTEGER(runs_)[0]; // The number of initial points
  const double *points = REAL(points_);
  const int maxit = INTEGER(maxit_)[0];
  const int nout = n - nin;
  double newcrit, crit=R_PosInf, themin=0.0;
  int oldclass0[n],iclass0[n], iclass[n], optclass0[n],optclass[n],I[n],rN[k],
    wend=0,itcounter=0,icls=0,nk,NfailC,AtLeastOneObs;
  int in,ik,iout,ip;
  double disttom[n],means[k][p], iWdisC[k],optmeans[k*p];
  // double optdisttom[n];
  
  GetRNGstate();

  for (int irun = 0 ; irun < runs ; irun ++)
    {
      //    ==== initial cluster center ====
      if (R_FINITE(points[0]))
	{
	  // use the initial cluster center selected by argument points
	  for (ik = 0 ; ik < k ; ik++) 
	    {
	      for (ip = 0 ; ip < p ; ip++)
		{
		  means[ik][ip] = points[p*ik+ip]; // The first p entries correspond to the first cluster center
		  //Rprintf("\v means[%d][%d]  %f",ik,ip,means[ik][ip]);
		}
	    }
	}else{
	
	// select randomly
	sample(k, n, rN); //void function  rN[k]; contains randomly selected n integers  
	for (ik = 0 ; ik < k ; ik++) 
	  {
	    for (ip = 0 ; ip < p ; ip++)
	      {
		means[ik][ip] = data[rN[ik]+n*ip];
		//Rprintf("\v means[%d][%d]  %f",ik,ip,means[ik][ip]);
	      }
	  }
	
      }

      //for (int i = 0 ; i < k ; i++) Rprintf(" selected %d",rN[i]) ;

      wend = 1; 
      itcounter = 0; 

      while(wend)
  	{
	  //Rprintf("\v\v\v itcounter %d \v\v\v",itcounter);
	  
  	  itcounter++;
  	  //Compute the distances between observations and the cluster centers 
  	  for (in=0; in < n;in++)
  	    {
	      //Rprintf("\v %d",in);
  	      for (ik=0 ; ik < k ; ik ++)
  		{
  		  iWdisC[ik] = 0.0;
  		  for (ip=0; ip<p; ip++ )
  		    iWdisC[ik] += pow( (data[in+n*ip] - means[ik][ip]), 2);

		  //Rprintf(" %f",iWdisC[ik]);
  		}
	     
  	      realmin(iWdisC, k, &icls, &themin); // void function
	      //Rprintf("cluster %d selected wdisc %f",icls,themin);
  	      disttom[in]= themin;
  	      iclass[in] = icls; //
  	      iclass0[in] = icls; // iclass0 == iclass at this point
  	      I[in] = in; // I = {0,1,2,3,...,n-1 }
  	    } // in


  	  // define outliers **it could be that some clusters do not have obs**
	  // NOTE: R_qsort_I is a void function which change the order of disttom and I
	  // disttom is sorted in increasing way.
	  // I originally contains the obs index (0-(n-1)). R_qsort_I changes the order of obs index
	  // according to the size of disttom: I[1] is the obs index corresponding to the smallest disttom
  	  //R_qsort_I(disttom, I, 1, n);

	  rsort_with_index(disttom, I, n);
	  //Rprintf("\v\v ");
	  //for (int i = 0 ; i < n; i++) Rprintf(" %d",I[i]);

  	  for (iout=0; iout < nout; iout++ ) // This loop does not run if nout = 0, nout + nin = n
  	    {
	      //Rprintf("out %d",I[nin+iout]);
  	      iclass0[I[nin+iout]] = -1; // iclass0[in] == -1 if the obs in is an outlier
	      
  	    }

  	  //stopping criteria if the class labels are the same then stop the iteration..
  	  if (itcounter >= maxit || identical(oldclass0,iclass0,n))
  	    {
	    
  	      wend = 0;
	    
  	    }else{

  	    NfailC = 0; // count the number of empty clusters 
  	    for (ik =0; ik < k; ik++ )
  	      {
  		// check if cluster ik is empty
  		// if a cluster is empty then cluster center is the first outlier obs
  		AtLeastOneObs = 0;
  		for (in = 0; in < n; in ++)
  		  {
  		    if (iclass0[in]==ik)
  		      {
  			AtLeastOneObs = 1; // find at least one observation in cluster k
  			break;
  		      }
  		  }

  		if (AtLeastOneObs==0)
  		  {
  		    // no observations found in this cluster ik
  		    for (ip = 0; ip < p ; ip++) means[ik][ip] = data[ I[n-1-NfailC]+ip*n ];
		    
  		    //NfailC++;

  		  }else{

  		  for (ip = 0; ip < p ; ip++)
  		    {
  		      nk = 0;
		      means[ik][ip] = 0.0;
  		      for (in = 0; in < n ; in++)
  			{
  			  if (iclass0[in]==ik) // iclass0 labels the trimmed obs index as -1
  			    {
  			      means[ik][ip] += data[in+ip*n];
  			      nk++;
  			    }
  			}
  		      means[ik][ip] /= nk;

		      //Rprintf("\v nk %d means[%d][%d]: %f",nk,ik,ip,means[ik][ip]);

  		    }
  		}
  		for (in = 0; in <n; in++) oldclass0[in] = iclass0[in];
  	      } //  ik
  	  } // if else(itcounter >= maxit || identical(oldclass0,iclass0))
  	} // while(wend)
      // check if this run with this starting point returns smaller within cluster sum of square
      newcrit = 0.0;
      for (in = 0 ; in < nin;in ++)
  	{
	  // note that disttom is ordered in increasing way. (the order does NOT correspond to the order of iclass0)
	  // top nout number of cases should not be included in the calculation of newcrit
	  newcrit += disttom[in];
  	}
      if (newcrit<=crit)
  	{
  	  // iclass does NOT contain -1. all the obs has cluster labels
  	  for (in = 0; in < n; in++ ) 
	    {
	      optclass[in] = iclass[in];
	      optclass0[in] = iclass0[in];
	      //optdisttom[in] = disttom[in];
	    }
  	  crit = newcrit;
  	  for (ik=0; ik < k ; ik++ )
  	    {
  	      for (ip=0;ip < p; ip++ ) optmeans[ik+k*ip] = means[ik][ip];
  	    }
  	}
    } //irun
  //

 
  SEXP res = PROTECT(allocVector(VECSXP, 5)); //result is stored in res
  SEXP labels = allocVector(INTSXP, n);
  SET_VECTOR_ELT(res, 0, labels);
  SEXP centers = allocVector(REALSXP, p*k);
  SET_VECTOR_ELT(res, 1, centers);
  SEXP oW = allocVector(INTSXP, nout);
  SET_VECTOR_ELT(res, 2, oW);
  SEXP WWSS = allocVector(REALSXP, 1);
  SET_VECTOR_ELT(res, 3, WWSS);
  SEXP classification = allocVector(INTSXP, n);
  SET_VECTOR_ELT(res, 4, classification);

  for (in = 0; in < n ; in++) 
    {
      INTEGER(labels)[in] = optclass[in];
      INTEGER(classification)[in] = optclass0[in];
    }
  
  for (ik = 0 ; ik < k ; ik++)
    {
      for (ip = 0; ip < p ; ip++)  REAL(centers)[ik+k*ip] = optmeans[ik+k*ip];
    }
  
  if (nout  > 0 )
    {
      iout=0;
      for (in = 0 ; in < n ; in++)
	{
	  if (optclass0[in] == -1)
	    {
	      INTEGER(oW)[iout] = in;
	      //Rprintf("outlied guy in %d",in);
	      iout++;
	      if (iout == nout) break;
	    }
	}
    } //if (nout  > 0 )
  REAL(WWSS)[0] = crit;
  
  PutRNGstate();

  UNPROTECT(1);
  return res;

} //trimkmeans function



SEXP RSKC_trimkmeans_missing(SEXP data_,
			     SEXP n_,
			     SEXP p_,
			     SEXP k_,
			     SEXP nin_,
			     SEXP runs_,
			     SEXP points_,
			     SEXP maxit_,
			     SEXP W_
			     )
{
  // the code of trimkmeans is modified so that
  // for alpha = 0, no case is trimmed. i.e. k-means clustering is performed     

    
  // === Importing values === // 
  const double *data=REAL(data_); // a vector of length n*p, containing the observations to be clustered 
  // data = c(Rdata[,1], Rdata[,2],...,Rdata[,ncol(Rdata)]); data[in+N*ip] = Rdata[in,ip]
  const int n = INTEGER(n_)[0];
  const int p = INTEGER(p_)[0]; 
  const int k = INTEGER(k_)[0];
  const int nin = INTEGER(nin_)[0];  // The number of trimmed observations
  const int runs = INTEGER(runs_)[0]; // The number of initial points
  const double *points = REAL(points_);
  const int maxit = INTEGER(maxit_)[0];
  const double *Ws = REAL(W_);
  const int nout = n - nin;
  double newcrit, crit=R_PosInf, themin=0.0;
  int oldclass0[n],iclass0[n], iclass[n], optclass0[n],optclass[n],I[n],rN[k],
    wend=0,itcounter=0,icls=0,nk,NfailC,AtLeastOneObs;
  int in,ik,iout,ip;
  double disttom[n],means[k][p], iWdisC[k],optmeans[k*p];
  double dist_ip,scaleW;
  // double optdisttom[n];
  
  GetRNGstate();

  // sum(Ws) does not change in this code, so we compute it at first
  double Wsum=0;
  for (ip=0;ip < p ; ip++) Wsum += Ws[ip];
  
  for (int irun = 0 ; irun < runs ; irun ++)
    {
      //    ==== initial cluster center ====
      if (R_FINITE(points[0]))
	{
	  // use the initial cluster center selected by argument points
	  for (ik = 0 ; ik < k ; ik++) 
	    {
	      for (ip = 0 ; ip < p ; ip++)
		{
		  means[ik][ip] = points[p*ik+ip]; // The first p entries correspond to the first cluster center
		  //Rprintf("\v means[%d][%d]  %f",ik,ip,means[ik][ip]);
		}
	    }
	}else{
	
	// select randomly
	sample(k, n, rN); //void function  rN[k]; contains randomly selected n integers  
	for (ik = 0 ; ik < k ; ik++) 
	  {
	    for (ip = 0 ; ip < p ; ip++)
	      {
		means[ik][ip] = data[rN[ik]+n*ip];
		//Rprintf("\v means[%d][%d]  %f",ik,ip,means[ik][ip]);
	      }
	  }
	
      }

      //for (int i = 0 ; i < k ; i++) Rprintf(" selected %d",rN[i]) ;

      wend = 1; 
      itcounter = 0; 


      while(wend)
  	{
	  //Rprintf("\v\v\v itcounter %d \v\v\v",itcounter);
	  
  	  itcounter++;
  	  //===Compute the distances between observations and the cluster centers 
	  // This part of codes is DIFFERENT from non-missing codes.
  	  for (in=0; in < n;in++)
  	    {
	      //Rprintf("\v %d",in);
  	      for (ik=0 ; ik < k ; ik ++)
  		{
  		  iWdisC[ik] = 0.0;
		  scaleW = 0.0;
  		  for (ip=0; ip<p; ip++ )
		    {
		      // dist_ip could be NA!! 
		      // note that: NA - double = NA
		      dist_ip = data[in+n*ip] - means[ik][ip];
		      if (!ISNA(dist_ip)){
			scaleW += Ws[ip];
			iWdisC[ik] += pow(dist_ip, 2);
		      }
		    }
		  if (scaleW == 0 ) 
		    {
		      Rprintf("WARNINGs: Too many missing values.. Try larger L1 value!!");
		    }
		  iWdisC[ik] *= Wsum/scaleW;
		  //Rprintf(" %f",iWdisC[ik]);
  		}
	     
  	      realmin(iWdisC, k, &icls, &themin); // void function
	      //Rprintf("cluster %d selected wdisc %f",icls,themin);
  	      disttom[in]= themin;
  	      iclass[in] = icls; //
  	      iclass0[in] = icls; // iclass0 == iclass at this point
  	      I[in] = in; // I = {0,1,2,3,...,n-1 }
  	    } // in


  	  // define outliers **it could be that some clusters do not have obs**
	  // NOTE: R_qsort_I is a void function which change the order of disttom and I
	  // disttom is sorted in increasing way.
	  // I originally contains the obs index (0-(n-1)). R_qsort_I changes the order of obs index
	  // according to the size of disttom: I[1] is the obs index corresponding to the smallest disttom
  	  //R_qsort_I(disttom, I, 1, n);

	  rsort_with_index(disttom, I, n);
	  //Rprintf("\v\v ");
	  //for (int i = 0 ; i < n; i++) Rprintf(" %d",I[i]);

  	  for (iout=0; iout < nout; iout++ ) // This loop does not run if nout = 0, nout + nin = n
  	    {
	      //Rprintf("out %d",I[nin+iout]);
  	      iclass0[I[nin+iout]] = -1; // iclass0[in] == -1 if the obs in is an outlier
	      
  	    }

  	  //stopping criteria if the class labels are the same then stop the iteration..
  	  if (itcounter >= maxit || identical(oldclass0,iclass0,n))
  	    {
	    
  	      wend = 0;
	    
  	    }else{

  	    NfailC = 0; // count the number of empty clusters 
  	    for (ik =0; ik < k; ik++ )
  	      {
  		// check if cluster ik is empty
  		// if a cluster is empty then cluster center is the first outlier obs
  		AtLeastOneObs = 0;
  		for (in = 0; in < n; in ++)
  		  {
  		    if (iclass0[in]==ik)
  		      {
  			AtLeastOneObs = 1; // find at least one observation in cluster k
  			break;
  		      }
  		  }

  		if (AtLeastOneObs==0)
  		  {
  		    // no observations found in this cluster ik
  		    for (ip = 0; ip < p ; ip++) means[ik][ip] = data[ I[n-1-NfailC]+ip*n ];
		    
  		    //NfailC++;

  		  }else{

  		  for (ip = 0; ip < p ; ip++)
  		    {
  		      nk = 0;
		      means[ik][ip] = 0.0;
  		      for (in = 0; in < n ; in++)
  			{
  			  if (iclass0[in]==ik) // iclass0 labels the trimmed obs index as -1
  			    {
			      if (!ISNA(data[in+ip*n]))
				{
				  means[ik][ip] += data[in+ip*n];
				  nk++;
				}
  			    }
  			}
		      if (nk > 0)
			means[ik][ip] /= nk;
		      else 
			means[ik][ip] = R_NaReal;
		      //Rprintf("\v nk %d means[%d][%d]: %f",nk,ik,ip,means[ik][ip]);

  		    }
  		}
  		for (in = 0; in <n; in++) oldclass0[in] = iclass0[in];
  	      } //  ik
  	  } // if else(itcounter >= maxit || identical(oldclass0,iclass0))
  	} // while(wend)
      // check if this run with this starting point returns smaller within cluster sum of square
      newcrit = 0.0;
      for (in = 0 ; in < nin;in ++)
  	{
	  // note that disttom is ordered in increasing way. (the order does NOT correspond to the order of iclass0)
	  // top nout number of cases should not be included in the calculation of newcrit
	  newcrit += disttom[in];
  	}
      if (newcrit<=crit)
  	{
  	  // iclass does NOT contain -1. all the obs has cluster labels
  	  for (in = 0; in < n; in++ ) 
	    {
	      optclass[in] = iclass[in];
	      optclass0[in] = iclass0[in];
	      //optdisttom[in] = disttom[in];
	    }
  	  crit = newcrit;
  	  for (ik=0; ik < k ; ik++ )
  	    {
  	      for (ip=0;ip < p; ip++ ) optmeans[ik+k*ip] = means[ik][ip];
  	    }
  	}
    } //irun
  //

 
  SEXP res = PROTECT(allocVector(VECSXP, 5)); //result is stored in res
  SEXP labels = allocVector(INTSXP, n);
  SET_VECTOR_ELT(res, 0, labels);
  SEXP centers = allocVector(REALSXP, p*k);
  SET_VECTOR_ELT(res, 1, centers);
  SEXP oW = allocVector(INTSXP, nout);
  SET_VECTOR_ELT(res, 2, oW);
  SEXP WWSS = allocVector(REALSXP, 1);
  SET_VECTOR_ELT(res, 3, WWSS);
  SEXP classification = allocVector(INTSXP, n);
  SET_VECTOR_ELT(res, 4, classification);

  for (in = 0; in < n ; in++) 
    {
      INTEGER(labels)[in] = optclass[in];
      INTEGER(classification)[in] = optclass0[in];
    }
  
  for (ik = 0 ; ik < k ; ik++)
    {
      for (ip = 0; ip < p ; ip++)  REAL(centers)[ik+k*ip] = optmeans[ik+k*ip];
    }
  
  if (nout  > 0 )
    {
      iout=0;
      for (in = 0 ; in < n ; in++)
	{
	  if (optclass0[in] == -1)
	    {
	      INTEGER(oW)[iout] = in;
	      //Rprintf("outlied guy in %d",in);
	      iout++;
	      if (iout == nout) break;
	    }
	}
    } //if (nout  > 0 )
  REAL(WWSS)[0] = crit;
  
  PutRNGstate();

  UNPROTECT(1);
  return res;

} //trimkmeans-missing function



SEXP tempC(SEXP vec1_, SEXP vec2_)
{
  // vector of length 4
   double *vec1 = REAL(vec1_); 
   double *vec2 = REAL(vec2_); 

  GetRNGstate();

  for (int i = 0 ;i < 4; i++)
    if (ISNA(vec1[i]-vec2[i]))
      {
	Rprintf("\v NA! i %d",i);
      }else{
      Rprintf("\v %f",vec1[i]-vec2[i]);
    }
  
  vec1[2] = R_NaReal;
  Rprintf("\v\v\v vec1[3] = R_NaReaL; \v\v\v");
    
  for (int i = 0 ;i < 4; i++)
    if (ISNA(vec1[i]-vec2[i]))
      {
	Rprintf("\v NA! i%d",i);
      }else{
      Rprintf("\v %f",vec1[i]-vec2[i]);
    }

  PutRNGstate();
  return (R_NilValue);
}
