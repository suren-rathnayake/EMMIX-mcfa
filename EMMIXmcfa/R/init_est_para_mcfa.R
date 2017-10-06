init_est_para.mcfa <- function (Y, g, q, start, init_method = "eigenA", ...) {

p <- ncol(Y)
n <- nrow(Y)
xi <- array(NA, c(q, g))
omega <- array(NA, c(q, q, g))
pivec <- array(NA, c(1, g))

if (init_method == "randA") {

  A <- matrix(rnorm(p * q), nrow = p, ncol = q)
  C <- chol(t(A) %*% A)
  A <- A %*% solve(C)

  D <- numeric(p)
  for (i in 1 : g) {
    indices  <- which(start == i)
    n_i <- length(indices)
    D <- D + (n_i - 1) * apply(Y[indices, ], 2, var) / (n - g)
  }
  D <- diag(D)

  sqrt_D <- diag(sqrt(diag(D)))
  inv_sqrt_D <- diag(1 / diag(sqrt_D))

  for (i in 1 : g) {

    indices  <- which(start == i)
    pivec[i] <- length(indices) / n
    uiT <- Y[indices, ] %*% A
    xi[, i]  <- apply(uiT, 2, mean)
    Si <- cov(Y[indices, ])
    eig_list    <- try(eigen(inv_sqrt_D %*% Si %*% inv_sqrt_D), TRUE)
    H <- eig_list$vectors
    sort.lambda <- sort(eig_list$values, decreasing = TRUE, index.return=TRUE)
    lambda <- sort.lambda$x
    ix_lambda <- sort.lambda$ix
    if (q == p) {
      sigma2 <- 0
    } else {
      lamlast <- lambda[(q + 1) : p]
      lamlast <- lamlast[lamlast > 0]
      sigma2 <- mean(lamlast, na.rm = TRUE)
    }
    if (q == 1) {
      omega[,, i] <- t(A) %*% sqrt_D %*% H[, ix_lambda[1 : q]] %*%
                      diag((lambda[1 : q] - sigma2), q) %*%
                      t(H[, ix_lambda[1 : q]]) %*% sqrt_D %*% A
    } else {
      omega[,, i] <- t(A) %*% sqrt_D %*% H[, ix_lambda[1 : q]] %*%
                        diag((lambda[1 : q] - sigma2))%*%
                        t(H[, ix_lambda[1 : q]]) %*% sqrt_D %*% A
    }
  }
}

if (init_method == "eigenA") {

  svd_tY <- svd(t(Y)/ sqrt(n - 1))
  A <- svd_tY$u[, 1 : q, drop = FALSE]

  for (i in 1 : g) {
    
    indices  <- which(start == i)
    pivec[i] <- length(indices) / n
    uiT <- Y[indices, ] %*% A
    xi[, i]  <- apply(uiT, 2, mean)
    omega[,, i] <- cov(uiT)
  }
  sqrt_d <- svd_tY$d[(q + 1) : p]
  sqrt_d <- sqrt_d[!is.na(sqrt_d)]
  sqrt_d <- sqrt_d[sqrt_d > 0]
  D <- diag(mean(sqrt_d^2))
}

colnames (xi) <- paste0("xi_", 1 : g)
dimnames(omega)[[1]] <- paste0("u_", 1 : q)
dimnames(omega)[[2]] <- paste0("u_", 1 : q)
dimnames(omega)[[3]] <- paste0("omega_", 1 : g)
colnames(pivec) <- paste0("pi_", 1 : g)

model <- list(g = g, q = q, pivec = pivec, A = A, xi = xi,
                      omega = omega, D = D)
class(model) <- 'mcfa'
return(model)
}
