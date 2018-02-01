rmcfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q 
  
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    
    stop("rmcfa require mvtnorm package. Please `install.packages(mvtnorm)`")
  }
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  for (i in 1 : g) {
    
    mu <- with(model, A %*% xi[, i])
    sigma <- with(model, A %*% omega[,, i] %*% t(A) + D)
    
    dat[which(n_mix == i), ] <- mvtnorm::rmvnorm(tabn[i], mean = mu, 
                                                 sigma = sigma, ...)
  }
  
  # n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  # tabn <- table(n_mix)
  # 
  # U <- matrix(NA, nrow = n, ncol = q)
  # for (i in 1 : g) {
  #  U[which(n_mix == i), ] <- rmvnorm(tabn[i], mean = model$xi[, i], sigma = model$omega[,, i], 
  #                                    method = "chol")
  # }
  # 
  # dat <- U %*% t(model$A) + rmvnorm(n, mean = rep(0, nrow(model$D)), sigma = model$D)
  dat
}


rmctfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q
  
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    
    stop("rmctfa require mvtnorm package. Please `install.packages(mvtnorm)`")
  }
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  for (i in 1 : g) {
    
    delta <- with(model, A %*% xi[, i])
    sigma <- with(model, A %*% omega[,, i] %*% t(A) + D)
    
    dat[which(n_mix == i), ] <- mvtnorm::rmvt(tabn[i], delta = delta, 
                                              sigma = sigma, 
                                              df = model$v[i], 
                                              type = "shifted", ...)
  }
  
 dat
}


# rmix <- function(n, model, ...) {
# 
# if (!requireNamespace("EMMIX", quietly = TRUE)) {
# 
#   stop("rmix require installing EMMIX package available from, \n",
#     "https://people.smp.uq.edu.au/GeoffMcLachlan/mix_soft/EMMIX_R/index.html")
# }
# 
# if (!(any(class(model) == "mcfa") || any(class(model) == "mctfa") ||
#       any(class(model) == "mfa") || any(class(model) == "mtfa"))) {
# 
#   stop("STOP:model must be of class mcfa, mctfa, mfa or mtfa")
# }
# 
# ncomp <- length(model$pivec)
# p <- dim(model$D)[1]
# sigma <- array(NA, c(p, p, ncomp))
# dat <- matrix(NA, nrow = n, ncol = p)
# 
# if ((any(class(model) == "mcfa")) || (any(class(model) == "mctfa"))) {
# 
#   if (any(class(model) == "mcfa")) {
#     dist <- "mvn"
#     dof <- NULL
#   }  
# 
#   if (any(class(model) == "mctfa")) {
#     dist <- "mvt"  
#     dof <- model$v
#   }
# 
#   mu <- with(model, A %*% xi)
#     
#   for (comp in 1 : ncomp) 
#     sigma[,, comp] <- with(model, A %*% omega[,, comp] %*% t(A) + D)
# 
# }
# 
# if ((any((class(model) == "mfa"))) || (any((class(model) == "mtfa")))) {
# 
#   if (any(class(model) == "mfa")) {
#     dist <- "mvn"
#     dof <- NULL
#   }  
# 
#   if (any(class(model) == "mtfa")) {
#     dist <- "mvt"  
#     dof <- model$v
#   }
# 
#   mu <- model$mu 
# 
#   if (model$sigma_type == "common") {
#      sigma[,, 1 : ncomp] <- with(model, B %*% t(B) + D)
# 
#   } else {
# 
#     if (model$D_type == "common") {
#       for (comp in 1 : ncomp)
#         sigma[,, comp] <- with(model, B[,, comp] %*%t(B[,, comp]) + D)  
#     } else {
#       for (comp in 1 : ncomp)
#         sigma[,, comp] <- with(model, B[,, comp] %*%t(B[,, comp]) + D[,, comp])   
#     }
#   }
# }
# 
# dat <- EMMIX::rdemmix2(n = n, p = p, g = ncomp, distr = dist, 
#                        pro = model$pivec, mu = mu, sigma = sigma, dof = dof)
# return(dat)
# }
