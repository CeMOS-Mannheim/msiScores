#' Pixel-wise correlation analysis
#'
#' Computes the coefficient of determination based on randomized sub-sampling of the observations.
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
#' @details For more details about the generated scores and for citation please refere to "Erich, Katrin, et al. Biochimica et Biophysica Acta (BBA)-
#' Proteins and Proteomics 1865.7 (2017): 907-915".
#'
#' @return returns a list of R2 values for each group indicated in the \code{sample.labels} object as well as overall R2 sampled from
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

randR2                             = function(data.mat, sample.labels, n.samples = 3000, box.plot = TRUE,
                                              col.vec = NULL)
{

       #// initiate the return object
       uniqueLabels                = unique(sample.labels)
       numGroups                   = length(uniqueLabels)

       returnList                  = setNames(object = vector("list", length= numGroups + 1), nm = c(uniqueLabels, "overall"))

       for ( i in 1 : numGroups)
       {

              idx                  = which(sample.labels == uniqueLabels[i])               # indices of the current group

              nSampleCutoff        = ifelse(length(idx) > 1000, 1000, length(idx))         # number of samples taken to speed up computation time

              idx                  = sample(x = idx, size = nSampleCutoff, replace = FALSE)      # sample the index vector

              corTable             = cor(x = t(data.mat[idx , ]))**2                          # find the correlations

              used.n.samples       = ifelse(length(corTable[lower.tri(corTable, diag = FALSE)]) < n.samples,
                                            length(corTable[lower.tri(corTable, diag = FALSE)]),
                                            n.samples)

              returnList[[i]]      = sample(corTable[lower.tri(corTable, diag = FALSE)], used.n.samples, replace = FALSE)

       }


       #// compute the correlations overall
       nSampleCutoff               = ifelse(length(sample.labels) > 1000, 1000, length(sample.labels))

       idx                         = sample(x = seq(1, length(sample.labels)), size = nSampleCutoff, replace = FALSE)      # sample the index vector

       corTable                    = cor(x = t(data.mat[idx , ]))**2                          # find the correlations

       used.n.samples       = ifelse(length(corTable[lower.tri(corTable, diag = FALSE)]) < n.samples,
                                     length(corTable[lower.tri(corTable, diag = FALSE)]),
                                     n.samples)




       returnList$overall          = sample(corTable[lower.tri(corTable, diag = FALSE)], used.n.samples, replace = FALSE)


       if(box.plot)
       {
              #// plot boxplots
              windows()

              if(is.null(col.vec)) {col.vec = rainbow(numGroups)}

              par(mar = c(10, 4.1, 4.1, 2.1))

              boxplot(returnList, col = c(col.vec, "white"),
                      ylab = "coeff. of Determination (R2)", xlab = "", boxwex = 0.5, notch= T,
                      main = "Randomized Coefficient of Determination per Group ", las = 2)

              for (igroup in 1 : (numGroups + 1))
              {
                     points(y = returnList[[igroup]], x  = rnorm(length(returnList[[igroup]]), mean = igroup, sd = 0.1),
                            pch = 20, col = rgb(0,0,0,0.07))
              }
       }



       return(returnList)


}
