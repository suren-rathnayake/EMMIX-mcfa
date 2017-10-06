is_valid_args.mcfa <- function (Y, g, q, itmax, nkmeans, nrandom,
																tol, initClust, init_para, init_method,
																conv_measure, warn_messages) {

#
if (any(is.na(Y))) {
	ERR <- "`Y' has missing value"
  class(ERR) <- "error"
	return(ERR)
}

if (!any(is.numeric(Y))) {
  ERR <- "`Y' has a non-numeric element"
  class(ERR) <- "error"
  return(ERR)
}

p <- ncol(Y)
if (is.null(p)) {
  ERR <- "the data must have more than one variable"
  class(ERR) <- "error"
  return(ERR)
}

if (p <= q) {
  ERR <- "The number of factors must be less than the number of variables"
  class(ERR) <- "error"
  return(ERR)
}

if (q < 0) {
  ERR <- "q must be a positive integer"
  class(ERR) <- "error"
  return(ERR)
}

if (g < 0) {
  ERR <- "g must be a positive integer"
  class(ERR) <- "error"
  return(ERR)
}

if (nkmeans < 0) {
  ERR <- "nkmeans must be a positive integer"
  class(ERR) <- "error"
  return(ERR)
}

if (nrandom < 0) {
  ERR <- "nrandom must be a positive integer"
  class(ERR) <- "error"
  return(ERR)
}

if (tol < 0) {
  ERR <- "tol must be a greater than zero"
  class(ERR) <- "error"
  return(ERR)
}

if ((init_method != "eigenA") && (init_method != "randA")) {
  ERR <- "init_method is either 'eigenA' and 'randA'"
  class(ERR) <- "error"
  return(ERR)
}

if ((conv_measure != "diff") && (conv_measure != "ratio")) {
  ERR <- "conv_measure needs to be either 'diff' or 'ratio'"
  class(ERR) <- "error"
	return(ERR)
}

if (class(warn_messages) != "logical") {
  ERR <- "warn_messages is either TRUE or FALSE'"
  class(ERR) <- "error"
  return(ERR)
}

ERR <- "TRUE"
return(ERR)
}

is_valid_args.mctfa <- function(Y, g, q, itmax, nkmeans, nrandom,
																tol, df_init, df_update, init_clust,
																init_para, init_method, conv_measure,
																warn_messages) {

#
ERR <- is_valid_args.mcfa(Y, g, q, itmax, nkmeans, nrandom, tol,
                          init_clust, init_para, init_method,
                          conv_measure, warn_messages)
if (class(ERR) == "error")
	 return(ERR)

ERR <- "TRUE"
return(ERR)
}
