\name{generate.featuretable}
\alias{generate.featuretable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Feature Table
}
\description{
Generate feature table from raw .mzXML files using XCMS's CentWave peak detection algorithm. Perform feature alignment for multi-file analysis.
}
\usage{
generate.featuretable(MS1directory, type, CAMERA = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MS1directory}{
Directory of folder containing all .mzXML files containing MS1 scans. This folder must contain only relevant .mzXML files of either DDA, DIA, or Fullscan sample(s).
}
  \item{type}{
For single sample analysis: type = "single", for multiple samples analysis: type = "multi".
}
  \item{ppm}{
Numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition. Default = 10.
}
  \item{peakwidth}{
Numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds. Default = c(10,120).
}
  \item{mzdiff}{
Numeric(1) representing the minimum difference in m/z dimension required for peaks with overlapping retention times; can be negative to allow overlap. During peak post-processing, peaks defined to be overlapping are reduced to the one peak with the largest signal. Default = 0.01.
}
  \item{snthresh}{
Numeric(1) defining the signal to noise ratio cutoff. Default = 6.
}
  \item{integrate}{
Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact. Default = 1.
}
  \item{prefilter}{
Numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I. Default = c(3,100).
}
  \item{noise}{
Numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection). Default = 100.
}
  \item{bw}{
Sample alignment: bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram. Default = 5.
}
  \item{mzwid}{
Sample alignment: width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples. Default = 0.015.
}
  \item{max}{
Sample alignment: maximum number of groups to identify in a single m/z slice. Default = 100.
}
  \item{CAMERA}{
Whether to perform CAMERA analysis to identify adducts and isotopes. "T" or "F".
}
}
\details{
The arguments "bw", "mzwid", and "max" are used for sample alignemnt only when processing multiple samples. They will not be used for single sample analysis even if user assigns them values.
}
\value{
Returns a feature table in dataframe format, containing the following columns in order: "mz", "rt", "rtmin", "rtmax", "intensity1". An additional column containing intensity values will be added for each additional sample in multi-sample analysis.
Also returns a vector MS1.files containing name of all mzXML files and XCMS object MS1data.
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
type <- "single"

featureTable <- generate.featuretable(MS1directory = MS1directory, type = type, ppm=10, peakwidth=c(10,120),
    mzdiff = 0.01, snthresh = 6, integrate = 1, prefilter = c(3,100), noise = 100, bw = 5, mzwid = 0.015,
    max = 100, CAMERA = F)

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
