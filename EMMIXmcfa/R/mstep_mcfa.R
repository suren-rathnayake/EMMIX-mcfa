Mstep.mcfa <- function (Y, g, q, pivec, A, xi, omega, D, ...) {

p <- ncol(Y)
if (is.null(p))
  p <- 1

n <- nrow(Y)
tau <- tau.mcfa(Y, g, q, pivec, A, xi, omega, D)
inv_D <- diag(1 / diag(D))
A1 <- array(0, c(p,q))
A2 <- array(0, c(q,q))
Di <- array(0, c(p))

for (i in 1 : g) {

  gamma <- (inv_D - inv_D %*% A %*%
              chol.inv(chol.inv(omega[,, i]) +
                t(A) %*% inv_D %*% A) %*%
                t(A) %*% inv_D) %*% A %*% omega[,, i]
  ti <- sum(tau[, i])
  XI <- xi[, i, drop = F]
  tY <- sweep(Y, 1, tau[, i], '*')
  Y_AXI <- sweep(Y, 2, A %*% XI, '-')
  tY_AXI <- sweep(Y_AXI, 1, tau[, i], '*')
  xi[, i] <- XI + t(gamma) %*% as.matrix(colSums(tY_AXI)) / ti
  omega[,, i] <- (diag(q) - t(gamma) %*% A) %*% omega[,, i] +
                  t(gamma) %*% t(Y_AXI) %*% tY_AXI %*% gamma / ti -
                  (XI - xi[,i]) %*% t(XI - xi[, i])
  A1 <- A1 + colSums(tY) %*% t(XI) + t(Y) %*% tY_AXI %*% gamma
  A2 <- A2 + (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti
  Di <- Di + diag(t(Y) %*% tY)
  pivec[i] <- ti / n
}

A <- try(A1 %*% chol.inv(A2))
if (class(A) == "try-error") {

  model <- "tried to invert a ill-conditioned or a singular matrix"
  class(model) <- "error"
  return(model)
}

D <- diag(Di - rowSums((A %*% A2) * A )) / n
model <- list(g = g, q = q, pivec = pivec, A = A,
              xi = xi, omega = omega, D = D)
return(model)
}
