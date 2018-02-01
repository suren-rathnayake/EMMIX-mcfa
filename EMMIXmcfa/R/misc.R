chol.inv <- function(x, ...){
  C <- chol(x)
  inv_x <- chol2inv(C)
  return(inv_x)
}

perm <- function (n, r, v = 1 : n){
  if (r == 1)
    matrix(v, n, 1)
  else if (n == 1)
    matrix(v, 1, r)
  else {
    X <- NULL
    for (i in 1 : n)
      X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i])))
    X
  }
}

# the minimum number of misallocations
err <- function(cls, hat_cls) {

  if (length(cls) != length(hat_cls))
    stop('Length of the two arguments should be equal')

  tcls <- rep(0, length(cls))
  labs <- unique(cls)
  for (j in 1 : length(labs))
    tcls[cls == labs[j]] <- j
  cls <- tcls

  tcls <- rep(0, length(hat_cls))
  labs <- unique(hat_cls)
  for(j in 1 : length(labs))
    tcls[hat_cls == labs[j]] <- j
  hat_cls <- tcls

  labs <- unique(c(hat_cls, cls))
  g <- length(labs)
  x <- perm(g, g)
  min_err <- Inf
  for (j in 1 : nrow(x)) {

    new_cls <- rep(0, length(cls))
    for (k in 1 : g)
      new_cls[cls == labs[k]] <- x[j, k]

    e <- sum(new_cls != hat_cls)
    if (e < min_err)
      min_err <- e
  }
  return(min_err)
}

ari <- function (cls, hat_cls) {

  if (length(cls) != length(hat_cls))
    stop('Length of the two arguments should be equal')

  tab <- table(cls, hat_cls)
  if (sum(diag(tab)) == length(cls))
    return(1)

  A <- sum(choose(tab, 2))
  B <- sum(choose(rowSums(tab), 2))
  C <- sum(choose(colSums(tab), 2))
  D <- choose(sum(tab), 2)

  ARI <- (A - B * C / D) / ( (B + C) / 2 - B * C / D)
  return(ARI)
}

as.num.fac <- function(x) {as.numeric(as.factor(x))}

plot.emmix <- function(x, ...) {

  if (x$q == 1) {

    plot(x$Fmat, 1 : length(x$Fmat),  xlim = range(x$Fmat), axes = FALSE,
          xlab = expression(widehat(u)[1]), ylab = "", type = "p",
          pch = if (x$g <= 5) {20 + as.numeric(x$clust)} else {
          as.numeric(x$clust)}, col = as.numeric(x$clust),
          bg =  as.numeric(x$clust))
    axis(side = 1)
  }

  if (x$q == 2)

    plot(x$Fmat[, c(1, 2)], col = as.numeric(x$clust),
          ylim = range(x$Fmat[, 2]), pch = if (x$g <= 5) {
          20 + as.numeric(x$clust)} else{as.numeric(x$clust)},
          bg = as.numeric(x$clust), xlim = range(x$Fmat[, 1]),
          xlab = expression(widehat(u)[1]),
          ylab=expression(widehat(u)[2]))

  if (x$q > 2)

    pairs(x$Fmat, col = as.numeric(x$clust), bg = as.numeric(x$clust),
            pch = if(x$g <= 5) {
              20 + as.numeric(x$clust)} else{as.numeric(x$clust)})
}

predict.mcfa <- function(object, Y, ...) {
  
  tau <- do.call("tau.mcfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max) 
  clust
}

predict.mctfa <- function(object, Y, ...) {
  
  tau <- do.call("tau.mctfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max) 
  clust
}

summary.emmix <- function(object, ...) {
  
  cat("Call:\n")
  print(object$call)
  summ <- cbind(num_comp = object$g,
                num_fac = object$q,
                log_like = object$logL,
                BIC = object$BIC)
  cat("\n")
  print(summ)
}

print.mcfa <- function(x, ...) {
  
  g <- x$g
  q <- x$q
  
  cat("Call:\n")
  print(x$call)
  
  cat("\n", paste0("Coefficients: (for i in 1 to ", g, ")", "\n"))
      
  cat("\nMixing Proportions\n")
  cat("pi_i : ", round(x$pivec, 3), "\n")
  
  cat("\nLoading Matrix\n")
  colnames(x$A) <- paste0("q_",  1 : q)
  cat("A: \n")
  print(x$A)
  
  cat("\nFactor Means (in columns)\n")
  
  colnames(x$xi) <- paste0("xi_", 1 : g)
  rownames(x$xi) <- paste0("q_",  1 : q)
  print(x$xi)
  
  cat("\nFactor Covariance Matrices\n")
  colnames(x$omega) <- paste0("q_",  1 : q)
  rownames(x$omega) <- paste0("q_",  1 : q)
  
  for(j in 1 : g) {
    
    cat(paste0("omega_", j, ":"), "\n")
    print(x$omega[,, j])
    cat("\n")
  }
  
  cat("Diagonal of the Error Covariance Matrix\n")
  cat("diag D: \n", diag(x$D), "\n")
}

print.mctfa <- function(x, ...) {
  
  g <- x$g
  q <- x$q
  
  cat("Call:\n")
  print(x$call)
  
  cat("\n", paste0("Coefficients: (for i in 1 to ", g, ")", "\n"))
  
  cat("\nMixing Proportions\n")
  cat("pi_i : ", round(x$pivec, 3), "\n")
  
  cat("\nLoading Matrix\n")
  colnames(x$A) <- paste0("q_",  1 : q)
  cat("A: \n")
  print(x$A)
  
  cat("\nFactor Means (in columns)\n")
  
  colnames(x$xi) <- paste0("xi_", 1 : g)
  rownames(x$xi) <- paste0("q_",  1 : q)
  print(x$xi)
  
  cat("\nFactor Covariance Matrices\n")
  colnames(x$omega) <- paste0("q_",  1 : q)
  rownames(x$omega) <- paste0("q_",  1 : q)
  
  for(j in 1 : g) {
    
    cat(paste0("omega_", j, ":"), "\n")
    print(x$omega[,, j])
    cat("\n")
  }
  
  cat("Diagonal of the Error Covariance Matrix\n")
  cat("diag D: \n", diag(x$D), "\n")
  
  cat("nu: \n", x$v, "\n")
  
}

factor_scores <- function(model, Y, ...) UseMethod("factor_scores")

factor_scores.mcfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mcfa", c(list(Y = Y), model))
  scores
}

factor_scores.mctfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mctfa", c(list(Y = Y), model))
  scores
}
