#' @title d score
#' @description  d score measures the similarity of two square matrix with same dimension. d_score equals 0 if the estimator is a column permutation of true value.
#' @param estimator A square matrix, usually an estimator of the \code{true_value} matrix.
#' @param true_value A square matrix, which the estimator is compared to.
#' @return A numeric value in [0,1].
#' @export
#' @examples d_score(diag(3), diag(3))

d_score <- function(estimator, true_value){
  d <- solve(estimator) %*% true_value
  q <- ncol(d)
  d_coeff <- 1/(2*q*(sqrt(q)-1))
  d_sum <- 0
  for(j in 1:q){
    max_ij <- max(abs(d[,j]))/sqrt(sum(d[,j]^2))
    max_ji <- max(abs(d[j,]))/sqrt(sum(d[j,]^2))
    temp <- 1/max_ij +1/max_ji-2
    d_sum <- d_sum + temp
  }
  d_score <- d_coeff*d_sum
  return(d_score)
}
