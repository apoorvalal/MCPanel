#' Convert a long panel to wide matrix (cleaner data.table version)
#' @param dt data.table in long panel format
#' @param unit_id unit id name
#' @param time_id time id name
#' @param treat   treatement name
#' @param outcome outcome name
#' @return list with treatment matrix W, outcome matrix Y, N0 (number of control units), and T0 (number of untreated periods)
#' @import data.table
#' @export
panelMatrices = function(dt, unit_id, time_id, treat, outcome) {
  dt = as.data.table(dt)
  # function to extract first column, convert it to rownames for a matrix
  matfy = function(X) {
    idnames = as.character(X[[1]])
    X2 = as.matrix(X[, -1])
    rownames(X2) = idnames
    X2
  }
  # reshape formula
  fmla = as.formula(paste0(unit_id, "~", time_id))
  # treatment matrix
  kv = c(unit_id, time_id, treat)
  W = matfy(dcast(dt[, ..kv], fmla, value.var = treat))
  # outcome matrix
  kv = c(unit_id, time_id, outcome)
  Y = matfy(dcast(dt[, ..kv], fmla, value.var = outcome))
  # move treated units to bottom of W and Y matrix
  treatIDs = which(rowSums(W) > 1)
  W = rbind(W[-treatIDs, ], W[treatIDs, , drop = FALSE])
  Y = rbind(Y[-treatIDs, ], Y[treatIDs, , drop = FALSE])
  N0 = nrow(W) - length(treatIDs)
  T0 = min(which(colSums(W) > 0)) - 1
  list(W = W, Y = Y, N0 = N0, T0 = T0)
}
