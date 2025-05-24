#' Function for distributional regressions
#' @param formula
#' @param test.function The function that generates a vector of test function outputs from a given data
#' @param data The list of data sets
#' @examples
#'
#' @export
#' @importFrom formula.tools get.vars
dlm <- function(formula, test.function, data, whitening = FALSE){
  call = match.call()
  terms = formula.tools::get.vars(formula, data = names(data))

  for(term in terms){
    if(!(term %in% names(data))) stop(paste0("Data set named ", term, " is not found."))
  }

  phi_matrix = lapply(data[names(data) %in% terms], FUN=test.function)
  mat0 <- simplify2array(lapply(phi_matrix, colMeans))
  colnames(mat0) <- names(phi_matrix)
  phi_matrix = do.call(rbind, phi_matrix)
  sigma_phi = cov(phi_matrix)
  sds = sqrt(diag(sigma_phi))
  mat = apply(mat0, 2, function(x){return(x/sds)})
  mat = mat[,terms]

  # Whitening Transformation
  if (whitening){
    E = eigen(sigma_phi)
    D = E$values
    D[abs(E$values)>1e-5] = 1/sqrt(D[abs(E$values)>1e-5])
    D[abs(E$values)<=1e-5] = 0
    sqrtinv_sigma_phi = E$vectors %*% diag(D) %*% t(E$vectors)
    mat <- sqrtinv_sigma_phi %*% mat0
    mat = mat[,terms]
  }

  # Generate New Formula (Reparametrization)
  new_formula = paste(c(sapply(terms[-c(1,2)], function(x) paste0("I(", x, "-", terms[2],")")), 0), collapse = " + ")
  new_formula = as.formula(paste(c(terms[1], new_formula), collapse = " ~ "))

  # Apply LM Function
  lm.fit = lm(new_formula, offset = mat[,terms[2]], data = as.data.frame(mat))
  summ = summary(lm.fit)

  # Add / Change Statistics
  summ$call = call
  new_fstat = ((summ$coefficients[,1]- 1/length(terms[-1])) %*% solve(summ$cov.unscaled) %*% (summ$coefficients[,1] - 1/length(terms[-1])))/summ$sigma^2/(length(terms)-2)
  summ$fstatistic[1] = new_fstat
  new_rsq = 1 - sum((summ$residuals)^2)/ sum((mat[,terms[1]] - rowMeans(mat[,terms[-1]]))^2)
  new_adj_rsq = 1 - (1-new_rsq) * (summ$df[1] + summ$df[2]) / summ$df[2]
  summ$r.squared = new_rsq
  summ$adj.r.squared = new_adj_rsq
  est = 1 - sum(summ$coefficients[,1])
  sd =  sqrt(sum(summ$cov.unscaled)) * summ$sigma
  tstat = est/sd
  pval = 2*(1-pt(abs(tstat), df = summ$df[2]))
  summ$coefficients = rbind(c(est, sd, tstat, pval), summ$coefficients)
  dimnames(summ$coefficients)[[1]] = c(terms[-1])
  summ$terms <- NULL
  summ$aliased <- any(summ$aliased)
  summ$X = mat

  summ$lm.fit = lm.fit
  class(summ) <- c("dlm", "summary.lm")

  return(summ)
}

#' Plot method for dlm objects
#' @param x A dlm object (output from dlm function)
#' @param ... Additional arguments passed to plot.lm
#' @export
plot.dlm <- function(x, ...) {
  plot(x$lm.fit, which = c(1, 2), ...)
}

