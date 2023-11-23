new_EL_list = function(
    fits,
    ull,
    se,
    np_pf,
    nreps
){
  stopifnot(is.data.frame(fits))
  stopifnot(is.data.frame(ull))
  stopifnot(is.data.frame(se))
  stopifnot(is.numeric(np_pf))
  stopifnot(is.numeric(nreps))
  out = list(
    fits = dplyr::as_tibble(fits),
    ull = dplyr::as_tibble(ull),
    se = dplyr::as_tibble(se),
    np_pf = np_pf,
    nreps = nreps
  )
  structure(out, class = "EL_list")
}
