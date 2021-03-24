\name{find.level3}
\alias{find.level3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Level 3 In-source Fragments
}
\description{
Creates a list of tables, where each table contains a parent feature and all its Level 3 in-source fragment features.
}
\usage{
find.level3(MS1directory, MS1.files, featureTable, type, peakCOR = 0.8, loss = 10, mz.tol = 0.02, rt.tol = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MS1directory}{
Directory containing .mzXML files in which MS1 scans were previously extracted from. This string is outputted by the generate.featuretable() function.
}
  \item{MS1.files}{
String vector containing names of .mzXML files in which MS1 scans were previously extracted from. This vector is outputted by the generate.featuretable() function.
}
  \item{featureTable}{
FeatureTable containg MS2 fragment and annotation information. This table is produced by the feature.annotation() function.
}
  \item{type}{
For single sample analysis: type = "single", for multiple samples analysis: type = "multi".
}
  \item{peakCOR}{
MS1 EIC Peak-peak Pearson correlation threshold to determine if two features are co-eluting, and thus a Level 1 in-source fragment. Default value is 0.8.
}
  \item{loss}{
Mass of in-source fragment has to be "loss" dalton lower than that of its parent. Default value for loss is 10 daltons.
}
  \item{mz.tol}{
Mass range used to extract feature EIC. Default value is 0.02 daltons.
}
  \item{rt.tol}{
Retention time range used to extract feature EIC. Default value is 15 seconds.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list object containing dataframes, where each dataframes contains a parent feature and all its level 3 in-source fragment features.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sam Shen, Jian Guo, Tao Huan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(ISFrag)
MS1directory <- "X:/Users/Sam_Shen/ISFtest20210127/HILIC(+)/HILIC(+)3/fullscan"
MS2directory <- "X:/Users/Sam_Shen/ISFtest20210127/HILIC(+)/HILIC(+)3/DDA"
type <- "single"
lib_directory <- "E:/SAM"
lib_name <- "convertedLibraryPos.msp"

featureTable <- generate.featuretable(MS1directory = MS1directory, type = type, ppm=10, peakwidth=c(10,120),
    mzdiff = 0.01, snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100, bw = 5, mzwid = 0.015,
    max = 100, CAMERA = F)
featureTable <- ms2.tofeaturetable(MS2directory = MS2directory, featureTable = featureTable)
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.7)
level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = type)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
