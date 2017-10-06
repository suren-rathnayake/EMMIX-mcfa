rmix <- function(n, model, ...) {

if (!requireNamespace("EMMIX", quietly = TRUE)) {

  stop("rmix require installing EMMIX package available from, \n",
    "https://people.smp.uq.edu.au/GeoffMcLachlan/mix_soft/EMMIX_R/index.html")
}

if (!(any(class(model) == "mcfa") || any(class(model) == "mctfa") ||
      any(class(model) == "mfa") || any(class(model) == "mtfa"))) {

  stop("STOP:model must be of class mcfa, mctfa, mfa or mtfa")
}

ncomp <- length(model$pivec)
p <- dim(model$D)[1]
sigma <- array(NA, c(p, p, ncomp))
dat <- matrix(NA, nrow = n, ncol = p)

if ((any(class(model) == "mcfa")) || (any(class(model) == "mctfa"))) {

  if (any(class(model) == "mcfa")) {
    dist <- "mvn"
    dof <- NULL
  }  

  if (any(class(model) == "mctfa")) {
    dist <- "mvt"  
    dof <- model$v
  }

  mu <- with(model, A %*% xi)
    
  for (comp in 1 : ncomp) 
    sigma[,, comp] <- with(model, A %*% omega[,, comp] %*% t(A) + D)

}

if ((any((class(model) == "mfa"))) || (any((class(model) == "mtfa")))) {

  if (any(class(model) == "mfa")) {
    dist <- "mvn"
    dof <- NULL
  }  

  if (any(class(model) == "mtfa")) {
    dist <- "mvt"  
    dof <- model$v
  }

  mu <- model$mu 

  if (model$sigma_type == "common") {
     sigma[,, 1 : ncomp] <- with(model, B %*% t(B) + D)

  } else {

    if (model$D_type == "common") {
      for (comp in 1 : ncomp)
        sigma[,, comp] <- with(model, B[,, comp] %*%t(B[,, comp]) + D)  
    } else {
      for (comp in 1 : ncomp)
        sigma[,, comp] <- with(model, B[,, comp] %*%t(B[,, comp]) + D[,, comp])   
    }
  }
}

dat <- EMMIX::rdemmix2(n = n, p = p, g = ncomp, distr = dist, 
                       pro = model$pivec, mu = mu, sigma = sigma, dof = dof)
return(dat)
}
