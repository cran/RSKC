\name{DutchUtility}
\alias{DutchUtility}
\alias{showDigit}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Multiple Features Data Set of Robert P.W. Duin.
}

\description{
This dataset consists of features of handwritten numerals (`0'--`9') (\emph{K}=10) extracted from a collection of Dutch utility maps.
%
Two hundred patterns per class (for a total of 2,000 (=\emph{N}) patterns)
have been digitized in binary images.
%
Raw observations are 32x45 bitmmaps, which are divided into
nooverlapping blocks of 2x3 and the number of pixels are counted in
each block. 
%
This generate \emph{p}=240 (16x15) variable, recodring the
normalized counts of pixels in each block and each element is an
integer in the range 0 to 6.
%
\code{rownames} of \code{DutchUtility} contains the true digits and \code{colnames} of it contains the position of the block matrix, from which the normalized counts of pixels are taken. 
}

\usage{
data(DutchUtility)
showDigit(index,cex.main=1)
}

\arguments{
\item{index}{
A scalar containing integers between 1 and 2000.
The function \code{ShowDigit} regenerates the sampled versions of the original images may be obtained (15x16 pixels).
(the source image (32x45) dataset is lost)
}
\item{cex.main}{
Specify the size of the title text with a numeric value of length 1.
}
}

\details{
The original dataset is freely available from USIMachine Learning Repository (Frank and Asuncion (2010)) website http://archive
.ics.uci.edu/ml/datasets.html.
}

\references{
Frank A, Asuncion A (2010). UCI Machine Learning Repository." \url{http://archive.ics.uci.edu/ml}.
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

data(DutchUtility)

truedigit <- rownames(DutchUtility)
(re <- RSKC(DutchUtility,ncl=10,alpha=0.1,L1=5.7,nstart=1000))
Sensitivity(re$labels,truedigit)
table(re$labels,truedigit)

## Check the bitmap of the trimmed observations 
showDigit(re$oW[1])
## Check the features which receive zero weights
names(which(re$weights==0))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
