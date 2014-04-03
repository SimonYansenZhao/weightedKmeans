fgkm <- function(x, k, strGroup, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10) 
{
  if (missing(k))
    stop("the number of clusters 'k' must be provided")
    
  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type
  nc <-ncol(x) # integer
  
  # get numbers of feature group
  numGroups <- length( unlist( strsplit( strGroup, ':') ) ) 
  
  Z <- .C("fgkm",
          x = as.double(as.matrix(x)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          strGroup=as.character(strGroup),
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          cluster = integer(nr),
          centers = double(k * nc),
          featureWeight = double(k * nc),
          groupWeight = double(k * numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
		  withiness = double(k),
          PACKAGE="weightedKmeans"
          )
       
  centers <- matrix( Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- matrix(Z$featureWeight)
  dim(featureWeight) <- c(k, nc)
  colnames(featureWeight) <- vars
  
  groupWeight <- matrix(Z$groupWeight)
  dim(groupWeight) <- c(k, numGroups )
  colnames(groupWeight) <- 1:ncol(groupWeight)
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  rownames(featureWeight) <- 1:nrow(featureWeight)
  rownames(groupWeight) <- 1:nrow(groupWeight)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withiness), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
		 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  dim(result$featureWeight) <- c(k, nc)
  dim(result$groupWeight) <- c(k, numGroups)
  
  class(result) <- c("kmeans", "fgkm")
  return(result)
}
