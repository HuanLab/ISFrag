\name{add.features}
\alias{add.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use or Add Existing Feature Table
}
\description{
Perform ISFrag analysis using user provided feature table in CSV format, or combine this feature table with XCMS feature table extracted in previous step.
}
\usage{
add.features(ft_directory, ft_name, featureTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ft_directory}{
Directory containing feature table in CSV format.
}
  \item{ft_name}{
Name of feature table in CSV format.
}
  \item{featureTable}{
Feature table from generate.featuretable() function. Optional argument.
}
}
\details{
The existing feature table must be in CSV format and contains only the following columns in order: "mz", "rt", "rtmin", "rtmax", "intensity". An additional column containing intensity values is needed for each additional sample in multi-sample analysis.
}
\value{
Returns a feature table in dataframe format, containing the following columns in order: "mz", "rt", "rtmin", "rtmax", "intensity". Multiple intensity columns for multi-sample analysis.
}
\references{
%% ~put references to the literature/web site here ~
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
MS1directory <- "X:/ISFtest20210127/HILIC(+)/HILIC(+)3/fullscan"
type <- "single"
ft_directory <- "X:/ISFtest20210127/HILIC(+)/HILIC(+)3/fullscan"
ft_name <- "HILIC(+)featuretable.csv"

# Combination of XCMS and additional user-provided feature table
featureTable <- generate.featuretable(MS1directory = MS1directory, type = type, ppm=10, peakwidth=c(10,120),
    mzdiff = 0.01, snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100, bw = 5, mzwid = 0.015,
    max = 100, CAMERA = F)
featureTable <- add.features(featureTable = featureTable, ft_directory = ft_directory, ft_name = ft_name)

# Only user-provided feature table
featureTable <- add.features(ft_directory = ft_directory, ft_name = ft_name)


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
