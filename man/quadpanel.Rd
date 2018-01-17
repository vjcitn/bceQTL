\name{quadpanel}
\alias{quadpanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
display eQTL behavior in 4 expression panels, ER+, ER- and adjacent normal for each
}
\description{
display eQTL behavior in 4 expression panels, ER+, ER- and adjacent normal for each
}
\usage{
quadpanel(snp, gene, negES, posES, negnormES, posnormES, extraX = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{
string
}
  \item{gene}{
string
}
  \item{negES}{
ER negative ExpressionSet, tumor tissue
}
  \item{posES}{
ER positive ExpressionSet, tumor tissue
%%     ~~Describe \code{posES} here~~
}
  \item{negnormES}{
adjacent normal to ER- tissue ExpressionSet
}
  \item{posnormES}{
adjacent normal to ER+ tissue ExpressionSet
}
  \item{extraX}{
string for enhanced x label
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(nhscancpcns_ann)
data(nhsnormpcns_ann)
data(tcann)
library(EnsDb.Hsapiens.v79)
erneg = nhscancpcns_ann[, which(nhscancpcns_ann$erpos == 0 & !is.na(nhscancpcns_ann$erpos))]
erpos = nhscancpcns_ann[, which(nhscancpcns_ann$erpos == 1 & !is.na(nhscancpcns_ann$erpos))]
normerneg = nhsnormpcns_ann[, which(nhsnormpcns_ann$erpos == 0 & !is.na(nhsnormpcns_ann$erpos))]
normerpos = nhsnormpcns_ann[, which(nhsnormpcns_ann$erpos == 1 & !is.na(nhsnormpcns_ann$erpos))]
quadpanel("rs704010", "TC0701086", erneg, erpos, normerneg, normerpos, " intronic to ZMIZ1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
