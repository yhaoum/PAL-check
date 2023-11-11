#' Combine parameter estimates with best unit log likelihoods
#'
#' @param x Object of class `EL_list`.
#' @param top_n Number of top unit-specific parameters to use.
#' @param top_n_shared Number of top shared parameters to use.
#' @param se_penalty Penalizes log likelihood in ranking based on
#'   `se*se_penalty`.
#' @param omit A vector of names for unit-specific parameters that shouldn't be
#'   combined based on unit log likelihood. Instead, they are packaged with the
#'   shared parameters of the given fit and their efficacy is based on the
#'   overall log likelihood.
#' @param is_spat Indicates whether parameters in `x` are for a spatPOMP model.
#'   Function will fail to return desired results if wrong.
#'
#'
#' @return Object of class `EL_list`. Log likelihoods are not necessarily
#'   correct after combining.
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
#' EL_out$fits
#' combine_top_fits(EL_out, top_n = 2)
#' }
combine_top_fits = function(
    x,
    top_n = 1,
    top_n_shared = 1,
    se_penalty = 0,
    omit = NULL,
    is_spat = FALSE
){
  unit_names = colnames(x$ull)
  U = length(unit_names)
  N = max(top_n, top_n_shared)

  stopifnot(N <= nrow(x$fits))

  score_total = x$fits$logLik - x$fits$se*se_penalty
  ranking_total = order(score_total, decreasing = TRUE)[1:top_n_shared]
  best_fits = dplyr::select(
    x$fits[ranking_total,], -"logLik", -"se"
  )
  if(is_spat){
    expp = spatCoef_to_pparams(best_fits[1,], units = unit_names)
  } else {
    expp = coef_to_pparams(best_fits[1,])
  }
  shared_names = names(expp$shared)
  specific_names = setdiff(rownames(expp$specific), omit)

  fits_cols = vector(length = length(unit_names), mode = "list")
  ull_cols = vector(length = length(unit_names), mode = "list")
  se_cols = vector(length = length(unit_names), mode = "list")
  cll = vector(length = length(unit_names), mode = "list")
  names(cll) = unit_names
  cll_se = vector(length = length(unit_names), mode = "list")
  names(cll_se) = unit_names

  for(z in seq_along(unit_names)){
    un = unit_names[[z]]
    score = x$ull[[un]] - se_penalty*x$se[[un]]
    ranking = order(score, decreasing = TRUE)[1:top_n]
    pn = if(is_spat) paste0(specific_names, z) else paste0(specific_names, "[",un,"]")
    recycle_vec = sort(rep_len(1:top_n, N))

    fits_cols[[z]] = x$fits[ranking, pn][recycle_vec,]
    ull_cols[[z]] = x$ull[ranking, un][recycle_vec,]
    names(ull_cols)[[z]] = un
    se_cols[[z]] = x$se[ranking, un][recycle_vec,]
    names(se_cols)[[z]] = un
    cll[[unit_names[[z]]]] =
      x$cll[[unit_names[[z]]]][ranking,, drop = FALSE][recycle_vec,, drop = FALSE]
    cll_se[[unit_names[[z]]]] =
      x$cll_se[[unit_names[[z]]]][ranking,, drop = FALSE][recycle_vec,, drop = FALSE]
  }

  fits = dplyr::bind_cols(fits_cols)
  ull = dplyr::bind_cols(ull_cols)
  se = dplyr::bind_cols(se_cols)

  opn = if(is.null(omit)) NULL else paste0(rep(omit, each = U),"[",unit_names,"]")
  for(param_name in c(shared_names, opn)){
    fits[[param_name]] =
      best_fits[sort(rep_len(1:top_n_shared, N)), param_name, drop = TRUE]
  }
  fits = fits |>
    dplyr::mutate(
      logLik = rowSums(ull)[sort(rep_len(1:top_n, N))],
      se = sapply(1:nrow(se), function(z)
        sqrt(sum(se[z,]^2))
      )[sort(rep_len(1:top_n, N))]
    )

  new_EL_list(
    fits = fits[colnames(x$fits)],
    ull = ull[colnames(x$ull)],
    se = se[colnames(x$se)],
    cll = cll[names(x$cll)],
    cll_se = cll_se[names(x$cll_se)],
    np_pf = x$np_pf,
    nreps = x$nreps
  )
}
