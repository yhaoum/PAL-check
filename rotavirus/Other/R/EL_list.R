#' Make a list containing the results of `eval_logLik()`
#'
#' @name EL_list
#' @param fits Data frame of fit results. Columns should be named `logLik`,
#'   `se`, followed by the parameter names.
#' @param ull Data frame of unit log likelihoods. Column names are unit names.
#' @param se Data frame of unit standard errors. Column names are unit names.
#' @param cll List of matrices containing estimated conditional log likelihoods
#'   for each unit.
#' @param cll_se List of matrices containing standard errors for estimated
#'   conditional log likelihoods for each unit.
#' @param np_pf Number of particles used by `eval_logLik()`.
#' @param nreps Number of replications used by `eval_logLik()`.
#'
#' @return The arguments in list form with class `EL_list`.
#'
NULL

new_EL_list = function(
    fits,
    ull,
    se,
    cll,
    cll_se,
    np_pf,
    nreps
){
  stopifnot(is.data.frame(fits))
  stopifnot(is.data.frame(ull))
  stopifnot(is.data.frame(se))
  stopifnot(is.list(cll))
  stopifnot(is.list(cll_se))
  stopifnot(is.numeric(np_pf))
  stopifnot(is.numeric(nreps))
  out = list(
    fits = dplyr::as_tibble(fits),
    ull = dplyr::as_tibble(ull),
    se = dplyr::as_tibble(se),
    cll = cll,
    cll_se = cll_se,
    np_pf = np_pf,
    nreps = nreps
  )
  structure(out, class = "EL_list")
}

validate_EL_list = function(x){
  if(!all(c("fits", "ull", "se", "np_pf", "nreps") %in% names(x))){
    stop(
      "`x` must have names 'fits', 'ull', 'se', 'np_pf', and 'nreps'",
      call. = FALSE
    )
  }
  if(!all(is.data.frame(x$fits), is.data.frame(x$ull), is.data.frame(x$se))){
    stop(
      "'fits', 'ull', and 'se' must be of class data.frame"
    )
  }
  if(!all(is.numeric(x$np_pf), is.numeric(x$nreps))){
    stop(
      "'np_pf' and 'nreps' must be of class numeric"
    )
  }
  if(!all(c("logLik", "se") %in% colnames(x$fits))){
    stop(
      "'fits' must have 'logLik' and 'se' columns",
      call. = FALSE
    )
  }
  x
}

#' @rdname EL_list
#' @export
EL_list = function(
  fits,
  ull,
  se,
  cll,
  cll_se,
  np_pf,
  nreps
){
  x = new_EL_list(fits, ull, se, cll, cll_se, np_pf, nreps)
  validate_EL_list(x)
}
