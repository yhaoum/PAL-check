#' Convert coef-style object to pparams-style object.
#'
#' Vector entries with name in the style of `param[unit]` are assumed to be
#' unit-specific whereas those with name in the style of `param` are assumed to
#' be shared.
#'
#' @param coef Vector in the style of `coef(panelPomp_obj)`. That is, a numeric
#'   vector with names styled as "`shared_parameter`" or
#'   "`specific_parameter[unit]`".
#'
#' @return A list of length 2 in the style of `pparams(panelPomp_obj)`. That is,
#'   a numeric vector with shared parameter names, and a matrix with specific
#'   parameters as row names and units as column names.
#' @export
#'
#' @examples
#' \dontrun{
#' coef_out = panelPomp::coef(AK_model())
#' coef_to_pparams(coef_out)
#' }
coef_to_pparams = function(coef){
  coef_tibble = tibble::tibble(
    full_param_name = names(coef),
    value = as.numeric(coef)
  ) |>
    tidyr::separate(
      .data$full_param_name,
      into = c("param", "unit"),
      sep = "\\[",
      fill = "right"
    ) |>
    dplyr::mutate(unit = gsub(pattern = "\\]", "", x = .data$unit))
  shared_tibble = dplyr::filter(coef_tibble, is.na(.data$unit))
  specific_tibble = dplyr::filter(coef_tibble, !is.na(.data$unit))
  shared_params = shared_tibble$value
  names(shared_params) = shared_tibble$param

  specific_params = specific_tibble |>
    tidyr::pivot_wider(names_from = "unit", values_from = "value") |>
    tibble::column_to_rownames(var = "param")

  list(shared = shared_params, specific = specific_params)
}
