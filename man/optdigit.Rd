\name{optd}
\alias{optd}
\alias{bitmapLab}
\alias{bitmapMat}
\alias{showbitmap}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Optical Recognition of Handwritten Digits of Frank A, Asuncion A (2010).
}

\description{
The dataset describes n = 1797 digits from 0 to 9 (K = 10), handwritten by 13 subjects. 
Raw observations are 32x32 bitmaps, which are divided into nonoverlapping
blocks of 4x4 and the number of on pixels are counted in each block. 
This generates p = 64 (= 8x8) variable, recording the normalized counts of pixels in each block and each element
is an integer in the range 0 to 16.
The row names of the matrix optd contains the true labels (between 0 and 9), and the column names of it contains the position of the block in original bitmap.
}

\usage{
data(optd)
showbitmap(index)
}

\arguments{
\item{index}{
A vector containing integers between 1 and 1797.
Given the observation indices, the \code{showbitmap} returns their original 32 by 32 bitmaps on R console. 
}
}

\details{
The original dataset is freely available from USIMachine Learning Repository (Frank and Asuncion (2010)) website http://archive
.ics.uci.edu/ml/datasets.html.
}

\references{
Frank A, Asuncion A (2010). UCI Machine Learning Repository." \url{http://archive.ics.uci.edu/ml.}
}
\author{
Yumi Kondo <y.kondo@stat.ubc.ca>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
\dontrun{

data(optd)

truedigit <- rownames(optd)
(re <- RSKC(optd,ncl=10,alpha=0.1,L1=5.7,nstart=1000))
Sensitivity(re$labels,truedigit)
table(re$labels,truedigit)

## Check the bitmap of the trimmed observations 
showbitmap(re$oW)
## Check the features which receive zero weights
names(which(re$weights==0))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
