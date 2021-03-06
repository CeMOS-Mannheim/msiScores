#' Statistical & multivariate scores for systematic workflow standardization
#'
#' Computes statistical & multivariate scores to address the question of reproducibility of sample preparation workflows. It is compatible
#' with mass spectrometry imaging as well as fingerprinting data (or any other multivariate data). This function is the main function in the package
#' which calls internally the other functions available within this package.
#'
#'
#' @param data.mat      A matrix holding the intensity information of all spectra with rows having the samples (pixels or spectra) and columns
#' denoting variables (m/z features in case of MS).
#'
#' @param sample.labels    A character vector of length \code{nrow(data.matrix)} holding the labels of the samples (pixels or spectra).
#'
#' @param n.pc              Number of principal components to be included into the calculation of the scores. Defaults to 3. Note that the
#' visulaization will be always based on 3 PCs.
#'
#' @param n.samples  The size of the sampling pool for foldchange and correlation analysis. Defaults to 3000.
#'
#' @param compute.r2 A logical whether to compute coefficient of determination (r2) via \code{randR2} method. Default is TRUE.
#'
#' @param compute.fc A logical whether to perform natural fold change analysis via \code{randR2} method. Default is TRUE.
#'
#' @details This function tries to evaluate the relative similarity between multivariate samples and can be used to assess reproducibility
#' of data generation for method development and standardization questions. It generates a 3D plot of the supplied dataset \code{data.mat}
#' color-coded by the specific grouping provided by \code{sample.labels}. Along with the plot, it generates unbiased scores that
#' describe the relative scatter and overlap of the samples supplied as a plot legend. These scores include covariance per group, mean absolute
#' deviation (MAD) per group, euclidean distances between the centers of the groups, overall within-class scatter, between-class scatter and a measure of
#' group overlap (J-overlap). Additionally, the function computes the coefficient of determination based on randomized sub-sampling of the observations by
#' calling \code{randR2} generating a boxplot in the process. The function also calls \code{natFC} function to evaluate  the so called natural
#' foldchange on a paired randomly chosen combinations of observations sampled from \code{data.mat}. For each pair of samples the function computes
#' \code{max(|5th FC precentile|, |95th FC percentile|)}. The computed values are then pooled into a list named based on the supplied groups as indicated
#' in \code{sample.labels} and a boxplot is plotted.
#'
#' Note that the within-class and between-class scatter as well as covariance scores are generated by taking the trace of the matrices rather than
#' the determinants as this is more stable.
#'
#' For more details about the generated scores and for citation please refere to "Erich, Katrin, et al. Biochimica et Biophysica Acta (BBA)-
#' Proteins and Proteomics 1865.7 (2017): 907-915".
#'
#' @return returns a list of four objects; the within-class scatter matrix \code{Ws}, the between-class scatter matrix \code{Bs} ,
#' a list of \code{R2} values per group and a list of \code{natFC} per group. Additionally it generates four plots; the variance per
#' principal component, a boxplot of \code{R2} values per group, a boxplot of \code{natFC} per group and an rgl 3D plot showing the
#' PCA space with the computed scores as a plot legend.
#'
#'
#'
#' @export
#'
#' @examples
#' data(exampleScores)
#' test <- msiScores::getScores(data.mat = mat, sample.labels = lab)
#'
#'
#' @author Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
#'
#' @references
#' Erich, Katrin, et al. "Scores for standardization of on-tissue digestion of formalin-fixed paraffin-embedded tissue in
#' MALDI-MS imaging." Biochimica et Biophysica Acta (BBA)-Proteins and Proteomics 1865.7 (2017): 907-915.
#'
#'
getScores     = function(data.mat, sample.labels, n.pc = 3, n.samples = 3000, compute.r2 = TRUE, compute.fc = TRUE)
{

       #// define the return list
       returnList                  = list()

       #// perfom pca

       cat("performing PCA .. ", "\n")
       pc                            = prcomp(x = data.mat, retx = TRUE, center = T, scale = F)


       #// plot the variances
       windows()
       pcChart(pc)




       coords                        = as.data.frame(pc$x[ , 1:n.pc]) # extract the points coordinates in pca



       #// RGL - plot 1st figure without ellipsoids



       #// find the center for each axis
       lim <- function(x){c(min(x), max(x)) * 1.1}
       xlim                          = lim(coords[,1])
       ylim                          = lim(coords[,2])
       zlim                          = lim(coords[,3])

       axes                          = rbind(c(median(xlim), ylim[1] *1, zlim[2]*1),
                                             c(xlim[2] * 1, median(ylim), zlim[1] * 1),
                                             c(xlim[2] * 1, ylim[1] * 1, median(zlim)))


       #// assign points to  groups
       pts                           = list()            # to hold the points coordinates in pca spcace
       ctrs                          = list()            # to hold the centers of groups
       uniqueLabels                  = unique(sample.labels)
       numGroups                     = length(uniqueLabels)

       for (igroup in 1 : numGroups)
       {
              pts[[igroup]]       = coords[which(sample.labels == uniqueLabels[igroup]), ]
              ctrs[[igroup]]      = colMeans(pts[[igroup]])
       }


       #// workout the colors
       #cdf                         = with(data.frame(labels = sample.labels),
       #                                   data.frame(labels = uniqueLabels,
       #                                              color = RColorBrewer::brewer.pal(n = length(uniqueLabels), name = "Set1")))

       tmpCol                      = ifelse(length(uniqueLabels) < 3, 3, length(uniqueLabels)) # used only for Ebrewer.pal function

       cdf                         = merge(data.frame(labels = sample.labels, stringsAsFactors = F),
                                           data.frame(labels = uniqueLabels,
                                           color = RColorBrewer::brewer.pal(n = tmpCol, name = "Set1")[1:numGroups],
                                           stringsAsFactors = F))

       cols                        = cdf$color
       uniqueCols                  = unique(cdf$color)











       #// find the set of all combinations of all groups - this is needed to draw the distances between groups
       if(numGroups > 1){
              distCombs            = combn(seq(1, numGroups), m = 2)

              #// append the above combination matirx with corresponding euclidean distances between respective centroids
              distCombs                     = data.frame(firstGroup = distCombs[1, ],
                                                         secondGroup = distCombs[2, ],
                                                         EuDistance = NA,
                                                         colour = sort(rainbow(ncol(distCombs) * 3), decreasing = TRUE)[1:ncol(distCombs)],

                                                         stringsAsFactors = FALSE)

              for (ictrs in 1 : nrow(distCombs))
              {
                     distCombs$EuDistance[ictrs] = sqrt(sum((ctrs[[distCombs$firstGroup[ictrs]]] - ctrs[[distCombs$secondGroup[ictrs]]])**2))
              }

       }







       withinGroupVariance           = vector("numeric", numGroups)

       #// To compute the within class scatter and between class scatter:


       Ws                            = DiscriMiner::withinSS(variables = coords, group = sample.labels)
       Bs                            = DiscriMiner::betweenSS(variables = coords, group = sample.labels)



       returnList$Ws               = Ws
       returnList$BS               = Bs



       #// MAD value

       MAD                         = c()

       for (iii in 1 : length(uniqueLabels))
       {
              idx           = which(sample.labels == uniqueLabels[iii])
              numPoints     = length(idx)
              coords_       = coords[idx, ]
              roiCenter     = colMeans(coords_[ , ])

              distances     = (rowSums(abs(coords_ - roiCenter))) #sqrt(rowSums((coords_ - roiCenter)**2))

              MAD[iii]      = mean(distances)
       }

       sdMAD                = sd(MAD, na.rm = TRUE)









       #// find the within-group variances
       for (icov in 1 : numGroups)
       {
              covMat                                  = cov(x = coords[which(sample.labels == uniqueLabels[icov]), ])
              withinGroupVariance[icov]               = sum(diag(covMat)) # based on the trace not deteminant

       }



       if(compute.r2)
       {
              #// compute the correlations per group
              cat("computing coefficient of determination on randomly chosen subset of spectra per group .. ", "\n")
              r2                          = randR2(data.mat = data.mat,
                                                              sample.labels = sample.labels,
                                                              n.samples = n.samples)




              returnList$R2List = r2

       }


       if(compute.fc)
       {
              #// compute foldchange
              cat("computing natural foldchange on randomly chosen subset of spectra per group .. ", "\n")
              fc                                               = natFC(data.mat = data.mat, sample.labels = sample.labels, n.samples = n.samples)


              returnList$FCList = fc
       }


       cat("plotting pca .. ", "\n")
       rgl::open3d()
       rgl::par3d(windowRect = c(0,23,1920,1040))
       rgl::rgl.viewpoint(theta = 15, phi = 3, fov = 50, zoom = 0.9)
       #par3d(cex = 2.5)

       rgl::par3d(cex = 2)
       rgl::plot3d(coords[ , 1:3], xlab = "", ylab = "", zlab = "", size = 20, pch = 19,
                   col = cols, box = F, alpha = 0.5, axes = F, type = "p")       # plot the pixel points - All

       rgl::axes3d(edges = c("x-+", "y+-", "z+-"), labels = TRUE, xlab = "PC1", ylab = "PC2", zlab = "PC3",tick = TRUE, nticks = 5,
                   box=F, expand = 1.03, lwd = 4, scale = 10, line = c(-60,0,0))
       rgl::grid3d(side = c( "y-+", "y++"))
       rgl::rgl.texts(axes[1, ], text = "PC1", adj = c(2, 2.8), color = "black", size = 3)
       rgl::rgl.texts(axes[2, ], text = "PC2", adj = c(-1.5, 0), color = "black", size = 3)
       rgl::rgl.texts(axes[3, ], text = "PC3", adj = c(-1.5, 2.5), color = "black", size = 3)
       #rgl.bg(color = "beige")

       if(numGroups > 1){
              palette(as.character(uniqueCols)) # adjust the palette cause for some reason legend3d depends on it
              #palette(as.character(cdf$color[which(cdf$labels %in% uniqueLabels)]))

       }


       if(numGroups > 1){
              rgl::legend3d("topleft" , horiz = FALSE,
                            title = "PCA Space",
                            c(uniqueLabels,
                              paste("covariance = ", round(withinGroupVariance, digits = 2)),
                              paste("MAD = ", round(MAD, digits = 2)),
                              paste("Eu distance = ", round(distCombs$EuDistance, digits = 2)),
                              paste("Within-class var = ", round(sum(diag(Ws)), digits = 2)),
                              paste("Between-class var = ", round(sum(diag(Bs)), digits = 2)),
                              paste("J-overlap = ", round(sum(diag(Ws))/sum(diag(Bs)), digits = 2))),
                            pch = c(rep(16, numGroups),
                                    rep(9, numGroups),
                                    rep(8, numGroups),
                                    rep(18, length(distCombs$EuDistance)),
                                    rep(15, 3)),
                            col = c(uniqueCols,
                                    uniqueCols,
                                    uniqueCols,
                                    distCombs$colour,
                                    rep("black", 3)),
                            cex=2.5, inset=c(0.02), bty = "n" ) #
       } else {

              rgl::legend3d("topleft" , horiz = FALSE,
                            title = "PCA Space",
                            c(uniqueLabels,
                              paste("covariance = ", round(withinGroupVariance, digits = 2)),
                              paste("MAD = ", round(MAD, digits = 2)),
                              paste("Within-class var = ", round(sum(diag(Ws)), digits = 2)),
                              paste("Between-class var is not available for 1 group"),
                              paste("J-overlap is not available for 1 group")
                              ),
                            pch = c(rep(16, numGroups),
                                    rep(9, numGroups),
                                    rep(8, numGroups),
                                    rep(15, 3)
                                    ),
                            col = c(uniqueCols,
                                    uniqueCols,
                                    uniqueCols,
                                    rep("black", 3)
                                    ),
                            cex=2.5, inset=c(0.02), bty = "n" )

       }




       #// draw the ellipsoids
       for (i in 1:numGroups) {

              ellips = rgl::ellipse3d(cov(pts[[i]][ , 1:3]), centre = ctrs[[i]][1:3], level = 0.95)


              rgl::shade3d(ellips, col = uniqueCols[i], alpha = 0.1, lit = 0.95)

              rgl::wire3d(ellips, col = "black", lit = FALSE, alpha = 0.2)


              rgl::rgl.spheres(ctrs[[i]][1:3], r = diff(xlim) / 40, col = uniqueCols[i] )

              #text3d(centers[[i]], texts = titles[i], cex = 2, adj = c(0,-1), col = "black" )
              #if (i == numGroups) {distin =  ctrs[[1]]} else {distin = ctrs[[i + 1]]}



       }

       #// draw the distance lines
       if(numGroups > 1) {
              for (i in 1 : nrow(distCombs))
              {
                     rgl::lines3d(rbind(ctrs[[distCombs$firstGroup[i]]][1:3], ctrs[[distCombs$secondGroup[i]]][1:3]),
                                  lwd = 5,  col = distCombs$colour[i], lty = "dashed", alpha = 0.8)
              }

       }



       rgl::play3d(rgl::spin3d(axis = c(0,1,0)), duration = 12 )

       #rgl::snapshot3d(filename = "pca.png")

       cat("\n")

       cat("Done. \n")

       return(returnList)


}
