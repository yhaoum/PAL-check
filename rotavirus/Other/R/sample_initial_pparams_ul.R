#' Use upper and lower bounds to sample initial parameters from box
#'
#' @param shared_box_specs `tbl` with `param`, `lower`, and `upper` columns.
#' @param specific_box_specs `tbl` with `param`, `lower`, and `upper` columns.
#' format of pparams.
#' @param units Character vector of unit names.
#' @param n_draws Number of initial parameter sets to draw.
#'
#' @return A list of parameters sets in the `pparams()` format.
#' @export
#'
#' @examples
#' \dontrun{
#' AK_mod = AK_model()
#' shared_bounds = tibble::tribble(
#' ~param, ~lower,     ~upper,
#' "mu",        0.02,     0.02
#' )
#' specific_bounds = tibble::tribble(
#'   ~param,       ~lower,        ~upper,
#'   "R0",             10,            60,
#'   "rho",           0.1,           0.9,
#'   "sigmaSE",      0.04,           0.1,
#'   "amplitude",     0.1,           0.6,
#'   "S_0",          0.01,          0.07,
#'   "E_0",      0.000004,        0.0001,
#'   "I_0",      0.000003,         0.001,
#'   "R_0",           0.9,          0.99,
#'   "sigma",          25,           100,
#'   "iota",        0.004,             3,
#'   "psi",          0.05,             3,
#'   "alpha",       0.935,          1.05,
#'   "cohort",        0.1,           0.7,
#'   "gamma",          25,           320
#' )
#' sample_initial_pparams_ul(
#'   shared_box_specs = shared_bounds,
#'   specific_box_specs = specific_bounds,
#'   units = names(AK_mod),
#'   n_draws = 3
#' )
#' }
sample_initial_pparams_ul = function(
    specific_box_specs,
    units,
    n_draws
){
  helper_df = tidyr::expand_grid(
    x = specific_box_specs$param,
    y = units
  ) |>
    dplyr::rename(param = "x", unit = "y")

  expanded_specific = specific_box_specs |>
    dplyr::right_join(helper_df, by = "param") |>
    dplyr::mutate(`param[unit]` = paste0(.data$param,"[",.data$unit,"]"))

  to_named_vec = function(x, name_col, val_col){
    named_vec = x[[val_col]]
    names(named_vec) = x[[name_col]]
    named_vec
  }

  initial_parameters_tbl = dplyr::bind_cols(
    pomp::runif_design(
      lower = to_named_vec(expanded_specific, "param[unit]", "lower"),
      upper = to_named_vec(expanded_specific, "param[unit]", "upper"),
      nseq = n_draws
    )
  )
  lapply(1:nrow(initial_parameters_tbl), function(z)
    coef_to_pparams(initial_parameters_tbl[z,])
  )
}
