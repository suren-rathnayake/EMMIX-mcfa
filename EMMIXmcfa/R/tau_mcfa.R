tau.mcfa <- function(Y, g, q, pivec, A, xi, omega, D, ...) {

p <- ncol(Y) 
if (is.null(p)) 
  p <- 1
n <- nrow(Y)
Fji <- array(NA, c(n, g))
inv_D <- diag(1 / diag(D))
for (i in 1 : g) {
  inv_S <- inv_D - inv_D %*% A %*%
            chol.inv(chol.inv(omega[,,i]) + t(A) %*% inv_D %*% A) %*%
            t(A) %*% inv_D
  logdetD <- log(det(as.matrix(omega[,,i]))) + sum(log(diag(D))) +
             log(det(chol.inv(omega[,,i]) + t(A) %*% inv_D %*% A))
  mahal_dist <- mahalanobis(Y, t(A %*% xi[, i, drop=FALSE]), inv_S, TRUE)
  Fji[,i] <- -0.5 * mahal_dist - (p / 2) * log(2 * pi) - 0.5 * logdetD
}

Fji <- sweep(Fji, 2, log(pivec), '+')
Fjmax <- apply(Fji, 1, max)
Fji <- sweep(Fji, 1, Fjmax, '-')
Fji <- exp(Fji)
tau <- sweep(Fji, 1, rowSums(Fji), '/')
return(tau)
}
