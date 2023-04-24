#' Imputation Estimators for Panel Data
#'
#' @description Imputes missing elements in a N X T panel data matrix and returns a list of all imputations
#' @param Y0  N X T matrix with missing values (must be either simultaneous or staggered adoption)
#' @param methods character vectors with names of estimators. Must be a subset of ("dfm", "knn", "kknn", "nuclear", "optspace", "softimpute", "hardimpute", "mc", "sc", "env", "enh"). Defaults to all of the above
#' @param dfm_args arguments to pass to dfms call (must match dfms::DFM function documentation)
#' @param sdfm_args arguments to pass to sparseDFM call
#' @param knn_args arguments to pass to filling::KNNimpute call
#' @param kknn_args arguments to pass to kernlab::kpca call
#' @return list with imputations corresponding with methods
#' @import dfms sparseDFM filling
#' @export
imputationY = function(
    Y0, # matrix with missing entries
    methods = c("dfm", "sdfm", "kknn", "knn", "did", "mc", "sc", "env", "enh"),
    # model parameters
    # dfm args
    dfm_args = list(r = 3, p = 1),
    # sdfm
    sdfm_args = list(r = 3, q = 1),
    # knn
    knn_args = list(k = 5),
    # kernel knn
    kknn_args = list(
      kernel = "anovadot",
      kpar = list(sigma = 0.1),
      features = nrow(Y0)/2, k = 3
    )) {
  # initialise list
  result = list()
  ######################################################################
  ## dynamic factor models
  if ("dfm" %in% methods) {
    dfm_args[['X']] = t(Y0)
    dfm_status = try({
      est_model_DFM = do.call(DFM, dfm_args)
      est_model_DFM$anyNA = FALSE
      result[['DFM']] = t(fitted(est_model_DFM))
    })
    if (is.error(dfm_status)) result['DFM'] = NA
  }
  # sparse DFM
  if ("sdfm" %in% methods) {
    sdfm_args[['X']] = t(Y0)
    sdfm_status = try({
      est_mod_SDFM = do.call(sparseDFM, sdfm_args)
      result[['sparse DFM']] = t(est_mod_SDFM$data$fitted.unscaled)
    })
    if (is.error(sdfm_status)) result['sparse DFM'] = NA
  }
  ######################################################################
  # filling matrix completion
  if ("knn" %in% methods) result[['knn']] = filling::fill.KNNimpute(Y0, k = knn_args$k)$X
  # nuclear norm penalisation w/o FEs
  if ("nuclear" %in% methods) result[['nuclear norm']] = filling::fill.nuclear(Y0)$X
  if ("optspace" %in% methods) result[['optspace']] = filling::fill.OptSpace(Y0)$X
  if ("softimpute" %in% methods) result[['softimpute']] = filling::fill.SoftImpute(Y0)$X
  if ("hardimpute" %in% methods) result[['hardimpute']] = filling::fill.HardImpute(Y0)$X
  if ("svdimpute" %in% methods) result[['svdimpute']] = filling::fill.SVDimpute(Y0, k = knn_k)$X
  ## kernel knn - argument handling happens internally
  if ("kknn" %in% methods) result[['kernel KNN']] = KKNN(Y0, kknn_args)
  ######################################################################
  # MCNNM packages require complete matrices, so replace NAs with 0
  mask = 1 * (!is.na(Y0))
  Y_obs = replace(Y0, mask == 0, 0)
  # did
  if ("did" %in% methods) result[['did']] = MCPanel::DID(Y_obs, mask)
  ## MC
  if ("mc" %in% methods) {
    mc_status = try({
      result[['matrix completion']] = mcnnm_cv(Y_obs, mask) |> predict_mc()
    })
    if (is.error(mc_status)) result[['matrix completion']] = NA
  }
  ## synth
  if ("sc" %in% methods) {
    sc_status = try({
      result[['synthetic control']] = adh_mp_rows(Y_obs, mask)
    })
    if (is.error(sc_status)) result[['synthetic control']] = NA
  }
  ## EN - V
  if ("env" %in% methods) {
    env_status = try({
      result[['elastic net (V)']] = en_mp_rows(Y_obs, mask, num_alpha = 5)
    })
    if (is.error(env_status)) result[['elastic net (V)']] = NA
  }
  ## EN - H
  if ("enh" %in% methods) {
    enh_status = try({
      result[['elastic net (H)']] = t(en_mp_rows(t(Y_obs), t(mask), num_alpha = 5))
    })
    if (is.error(enh_status)) result[['elastic net (H)']] = NA
  }
  class(result) = "imputeY"
  return(result)
}

# %% ####################################################
#' Imputation estimators for ATT
#' @description meta-function to perform reshape and fit multiple imputation methods
#' @param df dataframe
#' @param unit_id unit identifier
#' @param time_id time identifier
#' @param treatment treatment identifier
#' @param outcome outcome identifier
#' @param imputation_methods list of methods to pass to imputationY
#' @param ... additional arguments to pass to imputationY function
#' @return event study coefficients from each method
#' @export
imputation_ATT = function(
    df, unit_id, time_id, treatment, outcome,
    imputation_methods = c("did", "dfm", "sdfm", "knn", "kknn", "mc", "sc", "env", "enh"),
    ...
  ) {
  # call panelMat
  setup = panelMat(df, unit_id, time_id, treatment, outcome)
  # keep track of treated unit indices
  setup[['treat_indices']] = which(rowSums(setup$W) > 0)
  # fit imputer
  setup[['imputations']] = imputationY(setup$Y0, methods = imputation_methods, ...)
  # time-average of treated outcomes for entire panel
  setup[['treat_time_avgs']] = with(setup, colMeans(Y[treat_indices, , drop = FALSE]))
  # difference this from imputation
  setup[['event_study_estimates']] = do.call(cbind, lapply(setup$imputations,
    \(x) setup$treat_time_avgs - colMeans(x[setup$treat_indices, , drop = FALSE])
  ))
  class(setup) = "imputationATT"
  return(invisible(setup))
}

# %% ####################################################
#' print method for imputationATT
#' @export
print.imputationATT = \(x, prec = 3){
  cat("ATT estimates \n")
  round(colMeans(with(x, event_study_estimates[-(1:T0),])), prec)
}

# %% ####################################################
#' plot method for imputationATT
#' @export
plot.imputationATT = \(x, est_names = colnames(x$event_study_estimates), prec = 2){
  # subset if necessary
  ests = x$event_study_estimates[, est_names]
  # compute ATTs (added to legend)
  atts = round(colMeans(with(x, ests[-(1:T0), ])), prec)
  colours = RColorBrewer::brewer.pal(ncol(ests), "Paired")
  matplot(ests,
    type = 'l', col = colours, lwd = 2, lty = 2,
    ylab = "Coef. Estimate", xlab = "Time", axes = F
  )
  matpoints(ests, pch = 16, col = colours, lwd = 2, lty = 2)
  abline(h = 0, lty = 3)
  abline(v = x$T0, lty = 4)
  legend("bottomleft", paste0(1:ncol(ests), " ", colnames(ests), " : ", atts),
    col = colours, lty = 1, ncol = 2, cex = 0.9,
    x.intersp = 0.8, text.width = 6, lwd = 5
  )
  axis(2)
  axis(side = 1, at = 1:nrow(ests), labels = as.numeric(rownames(ests)))
}

