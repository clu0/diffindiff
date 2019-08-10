#' @title generate hermite basis
#' @description  Generates hermite basis for a given order, which we use for
#' non-parametric regression
#'
#' @param X an n by p matrix containing the input features
#' @param order order of hermite polynomials that we will use
#' @return a n by q matrix, which is a basis expansion of the input X
#' @export
generate.basis = function(X, order=3) {
  H = lapply(1:ncol(X), function(j) {
    sapply(1:order, function(k) hermite(X[,j], k, prob = TRUE) / sqrt(factorial(k)))
  })
  polys = lapply(1:order, function(r) {
    partitions = combn(r + ncol(X) -1, ncol(X) - 1,
                       function(vec) c(vec, r + ncol(X)) - c(0, vec) - 1)
    elems = sapply(1:ncol(partitions), function(iter) {
      part = partitions[,iter]
      idx = which(part > 0)
      elem = H[[idx[1]]][,part[idx[1]]]
      if (length(idx) > 1) {
        for (id in idx[-1]) {
          elem = elem * H[[id]][,part[id]]
        }
      }
      elem
    })
    scale(elems) / sqrt(ncol(elems)) / r
  })
  Reduce(cbind, polys)
}
