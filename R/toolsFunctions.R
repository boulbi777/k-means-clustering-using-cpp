###################################################################################
##' Cette fonction effectue la somme entre deux nombres
##'
##' @param n numeric. Taille de l'echantillon.
##' @param delta numeric. Parametre qui controle la separation des classes
##'
##' @examples
##' ech <- simuledata(50, 3)
##' @return list
##' @export
##'
simuledata <- function(n, delta){
  z <- sample(1:4, n, replace = TRUE)
  centers <- matrix(c(delta, delta, -delta, -delta, delta, -delta, delta, -delta), 4, 2)
  x <- matrix(rnorm(n * 2), n, 2)
  for (k in 1:4) x[which(z==k),] <- sweep(x[which(z==k),], 2, centers[k,], "+")
  list(x=x, z=z)
}

# n: taille echantillon
# delta: parametre pour separation des classes
# K: nbr de classes
# nbinit: nbr d initialisations de l algo des Kmeans
# nbCPU: nbr de CPU (parallelisation si >1)
# useC: boolean (appel a du code C si TRUE)

checkInputs <- function(x, K, nbinit, nbCPU, useC){
  if (!is.matrix(x)) stop("x must be a matrix")
  if (any(is.na(x))) stop("x cannot have missing values")
  if ((length(K)!=1) | (!is.numeric(K)) | (ceiling(K)!=K)) stop("K must be an integer of length 1")
  if ((length(nbinit)!=1) | (!is.numeric(nbinit))| (ceiling(nbinit)!=nbinit)) stop("nbinit must be an integer of length 1")
  if ((length(nbCPU)!=1) | (!is.numeric(nbCPU))| (ceiling(nbCPU)!=nbCPU)) stop("nbCPU must be an integer of length 1")
}

# Determine les points de departs
giveStarting <- function(x, K, nbinit)  replicate(nbinit, x[sample(1:nrow(x), K),], simplify = FALSE)
