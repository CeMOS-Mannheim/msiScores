% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getScores.R
\name{getScores}
\alias{getScores}
\title{Statistical & multivariate scores for systematic workflow standardization}
\usage{
getScores(
  data.mat,
  sample.labels,
  n.pc = 3,
  n.samples = 3000,
  compute.r2 = TRUE,
  compute.fc = TRUE
)
}
\arguments{
\item{data.mat}{A matrix holding the intensity information of all spectra with rows having the samples (pixels or spectra) and columns
denoting variables (m/z features in case of MS).}

\item{sample.labels}{A character vector of length \code{nrow(data.matrix)} holding the labels of the samples (pixels or spectra).}

\item{n.pc}{Number of principal components to be included into the calculation of the scores. Defaults to 3. Note that the
visulaization will be always based on 3 PCs.}

\item{n.samples}{The size of the sampling pool for foldchange and correlation analysis. Defaults to 3000.}

\item{compute.r2}{A logical whether to compute coefficient of determination (r2) via \code{randR2} method. Default is TRUE.}

\item{compute.fc}{A logical whether to perform natural fold change analysis via \code{randR2} method. Default is TRUE.}
}
\value{
returns a list of four objects; the within-class scatter matrix \code{Ws}, the between-class scatter matrix \code{Bs} ,
a list of \code{R2} values per group and a list of \code{natFC} per group. Additionally it generates four plots; the variance per
principal component, a boxplot of \code{R2} values per group, a boxplot of \code{natFC} per group and an rgl 3D plot showing the
PCA space with the computed scores as a plot legend.
}
\description{
Computes statistical & multivariate scores to address the question of reproducibility of sample preparation workflows. It is compatible
with mass spectrometry imaging as well as fingerprinting data (or any other multivariate data). This function is the main function in the package
which calls internally the other functions available within this package.
}
\details{
This function tries to evaluate the relative similarity between multivariate samples and can be used to assess reproducibility
of data generation for method development and standardization questions. It generates a 3D plot of the supplied dataset \code{data.mat}
color-coded by the specific grouping provided by \code{sample.labels}. Along with the plot, it generates unbiased scores that
describe the relative scatter and overlap of the samples supplied as a plot legend. These scores include covariance per group, mean absolute
deviation (MAD) per group, euclidean distances between the centers of the groups, overall within-class scatter, between-class scatter and a measure of
group overlap (J-overlap). Additionally, the function computes the coefficient of determination based on randomized sub-sampling of the observations by
calling \code{randR2} generating a boxplot in the process. The function also calls \code{natFC} function to evaluate  the so called natural
foldchange on a paired randomly chosen combinations of observations sampled from \code{data.mat}. For each pair of samples the function computes
\code{max(|5th FC precentile|, |95th FC percentile|)}. The computed values are then pooled into a list named based on the supplied groups as indicated
in \code{sample.labels} and a boxplot is plotted.

Note that the within-class and between-class scatter as well as covariance scores are generated by taking the trace of the matrices rather than
the determinants as this is more stable.

For more details about the generated scores and for citation please refere to "Erich, Katrin, et al. Biochimica et Biophysica Acta (BBA)-
Proteins and Proteomics 1865.7 (2017): 907-915".
}
\examples{
data(exampleScores)
test <- msiScores::getScores(data.mat = mat, sample.labels = lab)


}
\references{
Erich, Katrin, et al. "Scores for standardization of on-tissue digestion of formalin-fixed paraffin-embedded tissue in
MALDI-MS imaging." Biochimica et Biophysica Acta (BBA)-Proteins and Proteomics 1865.7 (2017): 907-915.
}
\author{
Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
}
