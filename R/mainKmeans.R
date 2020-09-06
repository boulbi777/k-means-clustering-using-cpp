
###################################################################################
##' Cette fonction effectue la somme entre deux nombres
##'
##' @param x numeric. Matrice des donnees (1 ligne par observation).
##' @param K numeric. Nombre de clusters.
##' @param nbinit numeric. Nombre d'initalisations aléatoires de l'algorithme.
##' @param nbCPU numeric. Nombre de CPU (si nbCPU>1, on parallelise).
##' @param useC boolean. Si TRUE, on appelle le code C++.
##'
##' @examples
##' set.seed(123)
##' ech <- simuledata(100, 3)
##'
##' set.seed(123)
##' T1 <- Sys.time()
##' res <- myKmeans(ech$x, 4, 200, 1, FALSE)
##' T2 <- Sys.time()
##' set.seed(123)
##' res2 <- myKmeans(ech$x, 4, 200, 2, TRUE)
##' T3 <- Sys.time()
##' all.equal(res, res2)
##' difftime(T3, T2)
##' difftime(T2, T1)
##'
##' @return list
##' @export
##'
myKmeans <- function(x, K, nbinit, nbCPU, useC){
  checkInputs(x, K, nbinit, nbCPU, useC)
  all.centers <- giveStarting(x, K, nbinit)
  all.res <- NULL
  if (nbCPU == 1){
    if (useC)
      all.res <- lapply(1:nbinit, function(it) singleKmeansC(x, all.centers[[it]]))
    else
      all.res <- lapply(all.centers, singleKmeans, x=x)
  }else{
    grappe <- makeCluster(nbCPU)
    clusterExport(cl=grappe, varlist = c("x"), envir = environment())
    clusterEvalQ(grappe, {require(packageKmeans)})
    if (useC){
      all.res <- parLapply(grappe, all.centers, function(u) singleKmeansC(x, u))
    }else{
      all.res <- parLapply(grappe, all.centers, function(u) singleKmeans(x, u))
    }
    on.exit(stopCluster(grappe))
  }
  all.res <- all.res[[which.min(sapply(all.res, function(u) u$criterion))]]
  if (useC) all.res$zactu <- as.numeric(all.res$zactu)
  all.res
}

# K: nbr de classes
# nbinit: nbr d initialisations de l algo des Kmeans
# nbCPU: nbr de CPU (parallelisation si >1)
# useC: boolean (appel a du code C si TRUE)
# centers: une initialisation (matrice avec K lignes) en dont les centres de classes

# Fonction principale des Kmeans:
# Un Kmeans avec initialisation donnee
singleKmeans <- function(x, centers){
  K <- nrow(centers)
  zprec <- rep(0, nrow(x))
  zactu <- apply(sapply(1:K, function(k) rowSums(sweep(x, 2, centers[k,], "-")**2)), 1, which.min)
  while (any(zprec!=zactu)){
    centers <- t(sapply(1:K, function(k) colMeans(x[which(zactu==k), , drop = FALSE])))
    zprec <- zactu
    zactu <- apply(sapply(1:K, function(k) rowSums(sweep(x, 2, centers[k,], "-")**2)), 1, which.min)
  }
  list(zactu = zactu,
       centers = centers,
       criterion = sum(apply(sapply(1:K, function(k) rowSums(sweep(x, 2, centers[k,], "-")**2)), 1, min)))
}
