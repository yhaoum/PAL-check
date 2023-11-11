#' Duplicate top fits from `EL_list`
#'
#' @param x Object of class `EL_list`.
#' @param out_length Number of parameter sets to output.
#' @param top_n Number of top fits to duplicate. If `combine = TRUE`, the number
#'   of top unit-specific parameters to combine.
#' @param top_n_shared Number of top shared parameters to use. Only used with
#'   `combine = TRUE`.
#' @param se_penalty Penalizes log likelihood in ranking based on
#'   `se*se_penalty`. Only used with `combine = TRUE`.
#' @param omit A vector of names for unit-specific parameters that shouldn't be
#'   combined based on unit log likelihood. Instead, they are packaged with the
#'   shared parameters of the given fit and their efficacy is based on the
#'   overall log likelihood. Only used with `combine = TRUE`.
#' @param combine Boolean specifying whether best specific fits should be
#'   combined.
#' @param units Character vector of unit names, which is necessary when
#'   duplicating the parameters of a `spatPomp` object. If the parameters belong
#'   to a `panelPomp` object, leave as `NULL`.
#'
#' @return List of parameters in the form of [panelPomp::pparams()].
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
#' duplicate_top_pparams(EL_out, out_length = 6, top_n = 2)
#' }
duplicate_top_pparams = function(
    x, #this is panelpomp object (a EL_list)
    out_length,
    top_n = 1,
    top_n_shared = 1,
    se_penalty = 0,
    omit = NULL,
    combine = FALSE,
    units = NULL
){
  is_spat = !is.null(units)
  if(ncol(x$ull) == 1 | combine == FALSE){
    grabbed_params = grab_top_fits(x, top_n = top_n)$fits |>
      dplyr::select(-"logLik", -"se")
  } else {
    grabbed_params = combine_top_fits(
      x,
      top_n = top_n,
      top_n_shared = top_n_shared,
      se_penalty = se_penalty,
      omit = omit,
      is_spat = is_spat
    )$fits |>
      dplyr::select(-"logLik", -"se")
  }
  top_params = dplyr::slice(
    grabbed_params,
    rep_len(1:nrow(grabbed_params), length.out = out_length)
  )
  lapply(1:nrow(top_params), function(z){
    if(!is_spat)
      coef_to_pparams(top_params[z,])
    else
      spatCoef_to_pparams(top_params[z,], units)
  })
}
