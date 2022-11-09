#' @title Blind Source Separation Over Space
#' @description  BSSS estimates the mixing matrix of blind source separation model for multivariate spatial data.
#' @param x A numeric matrix of dimension c(n, p), where the p columns correspond to the entries of the random field and the n rows are the observations.
#' @param coord A numeric matrix of dimension c(n,2) where each row represents the coordinates of a point in the spatial domain. Only needed if the argument kernel_list is NULL.
#' @param kernel_type A string indicating which kernel function to use. Either 'ring', 'ball' or 'gauss'.
#' @param kernel_parameter A numeric vector that gives the parameters for the kernel function. At least length of one for 'ball' and 'gauss' or two for 'ring' kernel.
#' @param kernel_list List of spatial kernel matrices with dimension c(n,n). Can be computed by the function \code{\link[SpatialBSS]{spatial_kernel_matrix}}.
#' @details BSSS estimates the mixing matrix by combining the information of all local covariance matrices together and conduct eigenanalysis.
#' @return BSSS returns a list, including the estimation of maxing matrix, the estimated latent field, and eigenvalues of matrix W for validating the estimation. Larger gaps among first few eigenvalues of matrix W strengthens the validity of estimation. See  Zhang, Hao and Yao (2022) <arXiv:2201.02023> for details.
#' @import SpatialBSS expm
#' @export
#' @examples
#'
#' \donttest{
#' sample_size <- 500
#' coords <- runif(sample_size * 2) * 50
#' dim(coords) <- c(sample_size, 2)
#' dim <- 5 # specify the dimensionality of random variable
#' nu <- runif(dim, 0, 6) # parameter for matern covariance function
#' kappa <- runif(dim, 0, 2) # parameter for matern covariance function
#' zs <- gen_matern_gaussian_rf(coords=coords, dim=dim, nu=nu, kappa=kappa)
#' mix_mat <- diag(dim) # create a diagonal matrix as the mixing matrix
#' xs <- t(mix_mat %*% t(zs))
#' example <- BSSS(xs, coords, 'ring', c(0,0.5,0.5,1,1,8))
#' d_score(example$mix_mat_est, mix_mat)
#' }



BSSS<- function(x, coord,  kernel_type, kernel_parameter, kernel_list = NULL){
  x <-x-mean(x)
  if (!missing(coord) &&
      !missing(kernel_parameter) && is.vector(kernel_parameter)){
    kernel_list <- SpatialBSS::spatial_kernel_matrix(coord, kernel_type = kernel_type,
                                                     kernel_parameters = kernel_parameter)
  }else{
    coord <- NULL
  }

  dim <- ncol(x)
  size <- nrow(x)
  B <- matrix(0, dim, dim)
  for (l in 1:size) {
    B <- B + tcrossprod(x[l,])
  }
  B <- B/size
  lambda_hat_b <- diag(eigen(B)$value)
  omega_hat_b <- eigen(B)$vectors
  ys <- t(solve(expm::sqrtm(lambda_hat_b)) %*% solve(omega_hat_b) %*% t(x))

  w_hat <- Reduce("+",lapply(SpatialBSS::local_covariance_matrix(x = ys, kernel_list = kernel_list, center=F),
                             function(lcov) tcrossprod(lcov)))/length(kernel_list)

  omega_est <- eigen(w_hat, symmetric=TRUE)
  omega_hat <- omega_est$vectors
  lambda_hat <- omega_est$values
  mix_mat_est <- omega_hat_b %*% expm::sqrtm(lambda_hat_b) %*% omega_hat
  latent_field_est <- tcrossprod(x, solve(mix_mat_est))
  return(list(mix_mat_est=mix_mat_est, latent_field_est=latent_field_est, w_eigenvalue=lambda_hat))
}
