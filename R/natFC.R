#' Pixel-wise foldchange analysis
#'
#' Performs foldchange analysis on randomly chosen, different combinations of observations of a given data matrix.
#'
#'
#' @param data.mat       A matrix holding the intensity information of all spectra with rows having the samples (pixels or spectra) and columns
#' denoting variables (m/z features in case of MS).
#'
#' @param sample.labels A character vector of length \code{nrow(data.matrix)} holding the labels of the samples (pixels or spectra).
#'
#' @param n.samples  The size of the sampling pool. Defaults to 3000.
#'
#' @param box.plot   A boolean used to plot the results in a boxplot. Defaults to TRUE.
#'
#' @param col.vec    A character vector specifying the colors of the groups with its length equal to number of unique groups in
#' \code{sample.labels}. If not provided, arbitrary colors will be assigned.
#'
#' @details The function evaluates the so called natural foldchange on a paired randomly chosen combinations of observations sampled from
#' \code{data.mat}. For each pair of samples the function computes \code{natFD = max(|5th FC precentile|, |95th FC percentile|)}. The computed
#' values are then pooled into a list named based on the supplied groups as indicated in \code{sample.labels}.
#'
#' For more details about the generated scores and for citation please refere to "Erich, Katrin, et al. Biochimica et Biophysica Acta (BBA)-
#' Proteins and Proteomics 1865.7 (2017): 907-915".
#'
#' @return returns a list of natFC values for each group indicated in the \code{sample.labels} object as well as overall natFC sampled from
#' all pixels of \code{data.mat}.
#'
#' @export
#'
#'
#'
#'
#' @author Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
#'
#'


#// compute the correlations per group

natFC                              = function(data.mat, sample.labels, n.samples = 3000, box.plot = TRUE,
                                              col.vec = NULL)
{

       #// initiate the return object
       uniqueLabels                = unique(sample.labels)
       numGroups                   = length(uniqueLabels)

       returnList                  = setNames(object = vector("list", length= numGroups + 1), nm = c(uniqueLabels, "overall"))

       for ( i in 1 : numGroups)
       {

              idx           = which(sample.labels == uniqueLabels[i])               # indices of the current group

              combs         = t(combn(idx, 2, simplify = TRUE))

              used.n.samples= ifelse(nrow(combs) < n.samples,
                                     nrow(combs),
                                     n.samples)

              combs         = combs[sample(seq(1, nrow(combs)), size = used.n.samples, replace = FALSE), ]




              for (ii in 1 : used.n.samples)
              {
                     fc                                 = numeric(ncol(data.mat))
                     sample1                            = data.mat[combs[ii, 1], ]
                     sample2                            = data.mat[combs[ii, 2], ]

                     biggerIdx                          = which(sample1 > sample2)
                     smallerIdx                         = which(sample1 < sample2)

                     fc[biggerIdx]                      = log2(sample1[biggerIdx] / sample2[biggerIdx])
                     fc[smallerIdx]                     = -log2(sample2[smallerIdx] / sample1[smallerIdx])


                     returnList[[i]]                    = c(returnList[[i]],
                                                            max(abs(quantile(fc, 0.05)), abs(quantile(fc, 0.95))))


              }




       }


       #// compute the correlations overall

       combs                              = t(combn(seq(1, nrow(data.mat)), 2, simplify = TRUE))

       combs                              = combs[sample(seq(1, nrow(combs)), size = n.samples, replace = FALSE), ]


       for (iii in 1 : n.samples)
       {
              fc                                 = numeric(ncol(data.mat))
              sample1                            = data.mat[combs[iii, 1], ]
              sample2                            = data.mat[combs[iii, 2], ]

              biggerIdx                          = which(sample1 > sample2)
              smallerIdx                         = which(sample1 < sample2)

              fc[biggerIdx]                      = log2(sample1[biggerIdx] / sample2[biggerIdx])
              fc[smallerIdx]                     = -log2(sample2[smallerIdx] / sample1[smallerIdx])


              returnList$overall                 = c(returnList$overall,
                                                     max(abs(quantile(fc, 0.05)), abs(quantile(fc, 0.95))))


       }




       if(box.plot)
       {
              #// plot boxplots
              windows()

              if(is.null(col.vec)) {col.vec = rainbow(numGroups)}

              par(mar = c(10, 4.1, 4.1, 2.1))

              boxplot(returnList, col = c(col.vec, "white"),
                      ylab = "Natural Foldchange (folds)", xlab = "", boxwex = 0.5, notch= T,
                      main = "Randomized natural Foldchange per Group ", las = 2)

              for (igroup in 1 : (numGroups + 1))
              {
                     points(y = returnList[[igroup]], x  = rnorm(length(returnList[[igroup]]), mean = igroup, sd = 0.1),
                            pch = 20, col = rgb(0,0,0,0.07))
              }
       }

       return(returnList)


}
