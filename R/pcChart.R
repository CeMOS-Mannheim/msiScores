#' Variances chart of PCA
#'
#' Computes the variances (cumulative and individual) per principal component and plots them.
#'
#'
#' @param pca.object       The PCA list of class \code{prcomp}.
#'
#'
#' @return returns a plot illustrating the variance explained by the principal components.
#'
#' @export
#'
#'
#'
#'
#' @author Denis Abu Sammour, \email{d.abu-sammour@hs-mannheim.de}
#'
#'

pcChart = function(pca.object)
       {
       x.var <- pca.object$sdev ^ 2
       x.pvar <- x.var/sum(x.var)


       par(mfrow=c(2,2))
       plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b', xlim = c(0, 10))
       plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b', xlim = c(0, 10))
       screeplot(pca.object)
       screeplot(pca.object,type="l")
       par(mfrow=c(1,1))
}
