library(formula.tools)
#' Function for distributional regressions
#' @param formula 
#' @param test.function The function that generates a vector of test function outputs from a given data point
#' @param data The list of data sets
#' @examples
#' 
#' @export
dlm <- function(formula, test.function, data, whitening = TRUE){
  
  call = match.call()
  terms = get.vars(formula, data = names(data))
  
  # Whitening Transformation
  phi_matrix = lapply(data[names(data) %in% terms], FUN=test.function)
  mat <- simplify2array(lapply(phi_matrix, colMeans))
  colnames(mat) <- names(phi_matrix)
  if (whitening){
    phi_matrix = do.call(rbind, phi_matrix)
    phi_matrix = apply(phi_matrix, 2, function(x){x/sd(x)})
    sigma_phi = cov(phi_matrix)
    E = eigen(sigma_phi)
    sqrtinv_sigma_phi = E$vectors %*% diag(1/sqrt(E$values)) %*% t(E$vectors)
    mat <- sqrtinv_sigma_phi %*% mat
  }
  mat <- mat[,terms]
  
  # Generate New Formula (Reparemetrization)
  new_formula = paste(c(sapply(terms[-c(1,2)], function(x) paste0("I(", x, "-", terms[2],")")), 0), collapse = " + ")
  new_formula = as.formula(paste(c(terms[1], new_formula), collapse = " ~ "))
  
  # Apply LM Function
  lm.fit = lm(new_formula, offset = mat[,terms[2]], data = as.data.frame(mat))
  summ = summary(lm.fit)
  summ$call = call
  #old_fstat = ((summ$coefficients[,1]) %*% solve(summ$cov.unscaled) %*% (summ$coefficients[,1]))/summ$sigma^2/(length(terms)-2)
  new_fstat = ((summ$coefficients[,1]- 1/length(terms[-1])) %*% solve(summ$cov.unscaled) %*% (summ$coefficients[,1] - 1/length(terms[-1])))/summ$sigma^2/(length(terms)-2)
  summ$fstatistic[1] = new_fstat
  new_rsq = 1 - sum((summ$residuals)^2)/ sum((mat[,terms[1]] - rowMeans(mat[,terms[-1]]))^2)
  new_adj_rsq = 1 - (1-new_rsq) * (summ$df[1] + summ$df[2]) / summ$df[2]
  summ$r.squared = new_rsq
  summ$adj.r.squared = new_adj_rsq
  
  # Include Offset
  est = 1 - sum(summ$coefficients[,1])
  sd =  sqrt(sum(summ$cov.unscaled)) * summ$sigma
  tstat = est/sd
  pval = 2*(1-pt(abs(tstat), df = summ$df[2]))
  summ$coefficients = rbind(c(est, sd, tstat, pval), summ$coefficients)
  dimnames(summ$coefficients)[[1]] = c(terms[-1])
  summ$X = mat
  
  # R-squared is based on using just the offset data
  return(summ)
}
