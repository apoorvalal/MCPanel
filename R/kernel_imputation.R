# %% kernel KNN estimator (instead of matching on raw Y, match on kernel factors)
#' Kernel KNN imputation
#' @param Y0 N X T matrix with missing entries for treated periods (output from panelMat).
#' @param k scalar number of neighbours to average over
#' @param kern kernel for kernel PCA
#' @param nf number of factors in kernel PCA to match over
#' @return complete matrix with weighted Kernel-KNN imputation entries
#' @import kernlab
#' @importFrom FNN get.knnx
#' @export
KKNN = function(Y0,
      kpca_args = list(
        kernel = "rbfdot",
        kpar = list(sigma = 0.1),
        features = 5, k = 3)
  ) {
  if (is.null(rownames(Y0))) rownames(Y0) = 1:nrow(Y0)
  T0 = min(which(is.na(apply(Y0, 2, sum)))) - 1
  treated_units = which(is.na(Y0[, T0 + 1]))
  # compute kernel principal commponents
  kpca_args[['x']] = Y0[,1:T0]
  kpcs = do.call(kpca, kpca_args) |>
    rotated() |>
    set_rownames(rownames(Y0))
  # KPCs for two groups
  treated_Y0 = kpcs[treated_units, , drop = FALSE]
  untreated_Y0 = kpcs[-treated_units, , drop = FALSE]
  # do KNN on KPCA
  knn_res = FNN::get.knnx(untreated_Y0, treated_Y0, k = kpca_args$k)
  # for each treated unit, return weighted average of *raw* outcomes (not KPCAs)
  weighted_knn = sapply(1:length(treated_units),
  	\(i) {
      donors = rownames(untreated_Y0[knn_res$nn.index[i, ], ])
      distinv = 1 / knn_res$nn.dist[i, ]
      wt = distinv / sum(distinv)
      apply(Y0[donors, , drop = FALSE], 2, weighted.mean, wt)
  	}) |>
    t() |>
    set_rownames(names(treated_units))
  Y0[treated_units, ] = weighted_knn
  dimnames(Y0) = NULL
  return(Y0)
}
