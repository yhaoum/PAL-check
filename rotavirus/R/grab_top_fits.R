#' Grab rows for top `n` fits based on log likelihood
#'
#' @param x Object of class `EL_list`.
#' @param top_n Number of rows to grab.
#' @param se_penalty Penalizes log likelihood in ranking based on
#' `se*se_penalty`.
#'
#' @return Object of class `EL_list`.
#' @export
#'
#' @examples
#' \dontrun{
#' AK_mod = AK_model()
#' EL_out = eval_logLik(
#'   list(AK_mod, AK_mod, AK_mod),
#'   ncores = 1,
#'   np_pf = 1,
#'   nreps = 2
#' )
#' grab_top_fits(EL_out)$fits[, 1:4]
#' }
grab_top_fits = function(
    x,
    top_n = 1,
    se_penalty = 0
  ){
  score = x$fits$logLik - se_penalty*x$fits$se
  ranking = order(score, decreasing = TRUE)
  new_EL_list(
    fits = x$fits[ranking,, drop = FALSE][1:top_n,, drop = FALSE],
    ull = x$ull[ranking,, drop = FALSE][1:top_n,, drop = FALSE],
    se = x$se[ranking,, drop = FALSE][1:top_n,, drop = FALSE],
    cll = x$cll |>
      lapply(function(y) y[ranking,, drop = FALSE][1:top_n,, drop = FALSE]),
    cll_se = x$cll_se |>
      lapply(function(y) y[ranking,, drop = FALSE][1:top_n,, drop = FALSE]),
    np_pf = x$np_pf,
    nreps = x$nreps
  )
}
