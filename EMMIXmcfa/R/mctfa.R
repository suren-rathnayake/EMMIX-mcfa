mctfa <- function(Y, g, q, ...) UseMethod("mctfa")

mctfa.default <- function (Y, g, q, itmax = 500, nkmeans = 20,
                          nrandom = 20, tol = 1.e-5, df_init = rep(30, g),
                          df_update = TRUE, init_clust = NULL,
                          init_para = NULL, init_method = 'randA',
                          conv_measure = 'diff', warn_messages = TRUE, ...) {

if (!is.matrix(Y))
  Y <- as.matrix(Y)
if (any(tolower(init_method) == c('e','eigen','eigena')))
  init_method <- 'eigenA'
if (any(tolower(init_method) == c('r', 'rand', 'randa')))
  init_method <- 'randA'
if (any(tolower(conv_measure) == c('d', 'diff')))
  conv_measure <- 'diff'
if (any(tolower(conv_measure) == c('r', 'ratio')))
  conv_measure <- 'ratio'

ERR <- is_valid_args.mctfa(Y, g, q, itmax, nkmeans, nrandom, tol,
                           df_init, df_update, init_clust, init_para,
                           init_method, conv_measure, warn_messages)
if (class(ERR) == "error") {
   stop(cat(ERR))
}
p <- ncol(Y)
n <- nrow(Y)

pb <- txtProgressBar(style = 3, char = '.')
prog <- 0
warn_msg <- NULL
maxLOGL <- -Inf
# Fit a model using given initial parameter estimates
if (!is.null(init_para)) {

  prog = 1/(1+nkmeans+nrandom)
  setTxtProgressBar(pb, prog)

  if (!check_para(p, q, g, init_para, "mctfa"))
    stop("incorrect specification of init_para", .call = FALSE)

  init_para <- init_para[c("g", "q", "pivec", "A", "xi", 
                            "omega", "D", "v")]
  init_para$df_update <- df_update

  estd_model <- est.mctfa(init_para = init_para, Y = Y, itmax = itmax,
                          tol = tol, conv_measure = conv_measure)

  if ((class(estd_model) == "mctfa")) {
    if (estd_model$logL > maxLOGL) {
      Hmodel <- estd_model
      maxLOGL <- Hmodel$logL
    }
    # cat(sprintf("\n g = %i, q = %i, init_para,  logL %8.4f \n",
    #                 g, q,  estd_model$logL))
  }
  if (class(estd_model) == "error") {
    when <- paste("init_para")
    what <- estd_model
    warn_msg <- cbind(when, what)
    colnames(warn_msg) <- c('when', 'what')
  }
}

# Use k-means, random, and given groupings of observations
# for parameter estimates.
if ((nkmeans!=0) || (nrandom!=0) || (!is.null(init_clust))) {

  # The set of all initial partitions
  initial_partitions <- start_clust(Y, g, init_clust, nkmeans, nrandom)
  maxinit <- ncol(initial_partitions)
  if (is.null(maxinit)) 
    maxinit <- 1
  # for progress bar
  if (prog == 0) {
    prog = 1/maxinit; tinit = maxinit
  } else {
    prog = prog + 1/maxinit; tinit <- 1 + maxinit
  }

  for (ii in 1 : maxinit) {

  # For some initialization methods, having a single sample cluster
  # can be a problem.
  if (min(table(initial_partitions[, ii])) == 1) {
    when <- paste("At start", ii)
    what <- "Initial partitioning has a single sample cluster."
    warn_msg <- rbind(warn_msg, cbind(when, what))
    next
  }

  # Initial estimates of model parameters (same as mcfa)
  init_model_para <- try(init_est_para.mcfa(Y, g, q,
                                            initial_partitions[, ii],
                                            init_method = init_method))

  if (class(init_model_para) == "try-error") {
    when <- paste("At start", ii)
    what <- "Failed to estimate initial parameters"
    warn_msg <- rbind(warn_msg, cbind(when, what))
    next
  }

  # Initial values for dof
  init_model_para$v <- df_init
  init_model_para$df_update <- df_update
  # EM steps
  estd_model <- est.mctfa(init_para = init_model_para, Y = Y,
                          itmax = itmax, tol = tol,
                          conv_measure = conv_measure)

  # keep the model with highest log-likelihood
  if (class(estd_model) == "mctfa") {
    if (estd_model$logL > maxLOGL) {
      Hmodel <- estd_model
      maxLOGL <- Hmodel$logL
    }
    # cat(sprintf("\n g = %i, q = %i, init %i logL %8.4f, maxlogL = %8.4f \n",
    #                 g, q, ii, estd_model$logL, maxLOGL))
  }

  if (class(estd_model) == "error") {
    when <- paste("At start", ii)
    what <- estd_model
    warn_msg <- rbind(warn_msg, cbind(when, what))
  }
  setTxtProgressBar(pb, prog)
  prog <- prog + 1 / tinit
  }
}
setTxtProgressBar(pb, 1)
close(pb)
if (!exists("Hmodel")) {
  cat("Failed to Estimate a Model. See Error Messages.")
  return(warn_msg)
}
# Make A^T A = I_q
CH <- chol(t(Hmodel$A) %*% Hmodel$A)
Hmodel$A <- Hmodel$A %*% solve(CH)
Hmodel$xi <- CH %*% Hmodel$xi
for (i in 1 : g) {
  Hmodel$omega[,, i] <- CH %*% Hmodel$omega[,, i] %*% t(CH)
}
d <- (g - 1) + p + q * (p + q) + g * q*  (q + 1) / 2 - q * q
Hmodel$BIC <- -2 * Hmodel$logL + d * log(n)
Hmodel$clust <- apply(Hmodel$tau, 1, which.max)
Hmodel <- append(Hmodel, do.call('factor_scores.mctfa',
                c(list(Y = Y), Hmodel)))
Hmodel$call <- match.call()
if (warn_messages == TRUE)
   Hmodel$warn_msg <- warn_msg
class(Hmodel) <- c("emmixmcfa", "mctfa")
return(Hmodel)
}
