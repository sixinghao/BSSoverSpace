#' @title Generating Gaussian random fields with Matern covariance function
#' @description  Generate Gaussian random fields with Matern covariance function
#' @param coords coordinate of target randon field to be generated
#' @param dim dimension of target randon field to be generated
#' @param nu parameter of matern covariance function
#' @param kappa parameter of matern covariance function
#' @return A data matrix with number of rows equal to `coords`, and number of columns equal to `dim`.
#' @import rSPDE
#' @importFrom stats dist rnorm
#' @export
#'
gen_matern_gaussian_rf <- function(coords, dim, nu, kappa){
  n <- nrow(coords)
  field <- matrix(0, n, dim)
  dist_coords <- dist(coords, method="maximum")
  corr_list <- vector("list", dim)
  for (i in 1:dim) {
    ri <- unlist(lapply(dist_coords,  function(x) matern.covariance(x, nu = nu[i], kappa = kappa[i], sigma = 1)))
    r_matrix <- matrix(0, n, n)
    r_matrix[lower.tri(r_matrix)] <- ri
    r_matrix[upper.tri(r_matrix)] <- t(r_matrix)[upper.tri(r_matrix)]
    diag(r_matrix) <- 1
    r_svd <- svd(r_matrix)
    r_dsq <- diag(sqrt(r_svd$d))
    r_sqrt <- r_svd$u %*% r_dsq %*% t(r_svd$u)
    corr_list[[i]] <- r_sqrt
  }

  for (i in 1:dim) {
    field_rnrom <- rnorm(n)
    field[,i] <- corr_list[[i]] %*% field_rnrom
  }
  return(field)
}
