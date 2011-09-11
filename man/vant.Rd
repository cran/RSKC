\name{v2002}
\alias{v2002}
\alias{ERv2002}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Contains the gene expression dataset of van 't et al. (2002)
}

\description{
This data set consisting on 4751 gene expressions for 78 primary breast tumors from
sporadic lymph-node-negative patients. 
The dataset was previously analyzed by van't Veer, Dal, van de Vijver, He, Hart, Mao, Peterse, van der Kooy, Marton, Witteveen, Schreiber,
Kerkhoven, Roberts, Linsley, Bernards, and Friend (2002) and it is publicly available at
http://bioinformatics.nki.nl/data.php. 
These authors applied a supervised classification
technique and found subset of 70 genes that helps to discriminate patients that develop
distant metastasis within 5 years from others.

The data matrix is found in \code{v2002}
The oestrogen receptor (ER) status of the 78 tumours are recorded in \code{ERv2002}
}
\usage{
data(v2002)
data(ERv2002)
}
%- maybe also 'usage' for other objects documented here.


\source{
The original dataset is freely available at
http://bioinformatics.nki.nl/data.php
}


\references{
Gene expression profiling predicts clinical outcome of breast cancer 
van 't Veer LJ, Dai H, van de Vijver MJ, He YD, Hart AA, Mao M, Peterse HL, van der Kooy K, 
Marton MJ, Witteveen AT, Schreiber GJ, Kerkhoven RM, Roberts C, Linsley PS, Bernards R, Friend SH 
Nature. 2002 Jan 31;415(6871):530-6. 
}
\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}



\examples{
data(v2002)
\dontrun{
data(v2002)
meds <- apply(v2002,1,median); mads<-apply(v2002,1,mad)
stand.d<-t(scale(t(v2002),center=meds,scale=mads))
set.seed(1)
(re <- RSKC(d=stand.d, ncl=2, alpha=10/78, L1 = 6, nstart=1000))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
