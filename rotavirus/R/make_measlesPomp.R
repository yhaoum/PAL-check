#' Make a panelPomp or spatPomp model using measles data
#'
#' @param model `panel_mechanics` or `spat_mechanics` object.
#' @param data List in the format of `twentycities`.
#' @param starting_pparams Parameters in the format of `pparams()` output. Set
#'   to NULL to assign NA values. Only for panelPomp models currently.
#' @param interp_method Method used to interpolate population and births.
#'   Possible options are `"shifted_splines"` and `"linear"`.
#' @param first_year Integer for the first full year of data desired.
#' @param last_year Integer for the last full year of data desired.
#' @param custom_obs_list List of observations where each element supplies
#'   observations for a different city. Useful when using simulated
#'   observations. Set to `NULL` to use real observations. Only for panelPomp
#'   models currently.
#' @param dt Size of the time step.
#'
#' @return A panelPomp or spatPomp object using the model and data supplied.
#' @export
#'
#' @examples
#' \dontrun{
#' make_measlesPomp(model_mechanics_001(), twentycities)
#' }
make_measlesPomp = function(
    model,
    data,
    starting_pparams = NULL,
    interp_method = c("shifted_splines", "linear"),
    first_year = 1950,
    last_year = 1963,
    custom_obs_list = NULL,
    dt = 1/365.25
){
  UseMethod("make_measlesPomp")
}
# TODO clean up arguments issue

#' @rdname make_measlesPomp
#' @export
make_measlesPomp.panel_mechanics = function(
    model,
    data,
    starting_pparams = NULL,
    interp_method = c("shifted_splines", "linear"),
    first_year = 1950,
    last_year = 1963,
    custom_obs_list = NULL,
    dt = 1/365.25
){
  rproc = model$rproc
  dmeas = model$dmeas
  rmeas = model$rmeas
  rinit = model$rinit
  pt = model$pt
  paramnames = model$paramnames
  states = model$states
  measles = data$measles
  demog = data$demog

  ## ----prep-data-------------------------------------------------
  units = unique(measles$unit)
  # Obs list
  dat_list = vector("list", length(units))
  # Population list
  demog_list = vector("list", length(units))
  for(i in seq_along(units)){
    dat_list[[i]] = measles |>
      dplyr::mutate(year = as.integer(format(date,"%Y"))) |>
      dplyr::filter(
        .data$unit == units[[i]] & .data$year >= first_year &
          .data$year < (last_year + 1)
      ) |>
      dplyr::mutate(
        time = julian(
          .data$date,
          origin = as.Date(paste0(first_year, "-01-01"))
        )/365.25 + first_year
      ) |>
      dplyr::filter(.data$time > first_year & .data$time < (last_year + 1)) |>
      dplyr::select("time", "cases")
    if(!is.null(custom_obs_list)) dat_list[[i]]$cases = custom_obs_list[[i]]

    demog_list[[i]] = demog |>
      dplyr::filter(.data$unit == units[[i]]) |>
      dplyr::select(-"unit")
  }
  ## ----prep-covariates-------------------------------------------------
  covar_list = vector("list", length(units))
  for(i in seq_along(units)){
    dmgi = demog_list[[i]]
    times = seq(from = min(dmgi$year), to = max(dmgi$year), by = 1/12)
    switch(interp_method[[1]],
      shifted_splines = {
        pop_interp = stats::predict(
          stats::smooth.spline(x = dmgi$year, y = dmgi$pop),
          x = times
        )$y
        births_interp = stats::predict(
          stats::smooth.spline(x = dmgi$year + 0.5, y = dmgi$births),
          x = times - 4
        )$y
      },
      linear = {
        pop_interp = stats::approx(
          x = dmgi$year,
          y = dmgi$pop,
          xout = times
        )$y
        births_interp = stats::approx(
          x = dmgi$year,
          y = dmgi$births,
          xout = times - 4
        )$y
      }
    )
    covar_list[[i]] = dmgi |>
    dplyr::reframe(
        time = times,
        pop = pop_interp,
        birthrate = births_interp
    )
    covar_list[[i]] = covar_list[[i]] |>
    dplyr::mutate(
      pop_1950 = dplyr::filter(
        covar_list[[i]], covar_list[[i]]$time == 1950
      )$pop
    )
  }
  for(i in seq_along(units)){
    log_pop_1950 = sapply(seq_along(units), function(x)
      log(covar_list[[x]][["pop_1950"]][[1]])
    )
    covar_list[[i]] = covar_list[[i]] |>
      dplyr::mutate(
        std_log_pop_1950 = (log(.data$pop_1950) - mean(log_pop_1950))/
          stats::sd(log_pop_1950),
        unit_num = i
      )
  }

  ## ----pomp-construction-----------------------------------------------
  lapply(seq_along(units), function(i){
    time = covar_list[[i]][[1]]
    dat_list[[i]] |>
      pomp::pomp(
        t0 = with(dat_list[[i]], 2*time[1] - time[2]),
        times = "time",
        rprocess = pomp::euler(rproc, delta.t = dt),
        rinit = rinit,
        dmeasure = dmeas,
        rmeasure = rmeas,
        covar = pomp::covariate_table(covar_list[[i]], times = "time"),
        accumvars = c("C","W"),
        partrans = pt,
        statenames = states,
        paramnames = paramnames
      )
  }) -> pomp_list
  names(pomp_list) = units

  ## ----panelPomp-construction-----------------------------------------------
  if(is.null(starting_pparams)){
    shared = as.numeric(rep(NA, length(model$shared_params)))
    specific = matrix(
      NA,
      nrow = length(model$specific_params),
      ncol = length(units)
    )
    class(specific) = "numeric"
    storage.mode(specific) = "numeric"
    rownames(specific) = model$specific_params
    colnames(specific) = units
    names(shared) = model$shared_params
  } else {
    shared = starting_pparams$shared
    specific = as.matrix(starting_pparams$specific)
    if(!setequal(names(shared), model$shared_params)){
      stop(
        "Starting shared parameters do not match parameters in model mechanics.",
        call. = FALSE
      )
    }
    if(!setequal(rownames(specific), model$specific_params)){
      stop(
        "Starting unit-specific parameters do not match parameters in model mechanics.",
        call. = FALSE
      )
    }
  }
  panelPomp::panelPomp(
    pomp_list,
    shared = shared,
    specific = specific
  )
}

#' @rdname make_measlesPomp
#' @export
make_measlesPomp.spat_mechanics = function(
    model,
    data,
    starting_pparams = NULL,
    interp_method = c("shifted_splines", "linear"),
    first_year = 1950,
    last_year = 1963,
    custom_obs_list = NULL,
    dt = 1/365.25
){
  if(!is.null(starting_pparams)){
    warning(
      "make_measlesPomp.spat_mechanics does not yet support non-default values for starting_pparams.",
      call. = FALSE
    )
  }
  if(!is.null(custom_obs_list)){
    warning(
      "make_measlesPomp.spat_mechanics does not yet support non-default values for custom_obs_list",
      call. = FALSE
    )
  }
  units = unique(data$measles$unit)
  U = length(units)
  cases_df = data$measles |>
    dplyr::mutate(year = as.integer(format(date,"%Y"))) |>
    dplyr::filter(.data$year >= first_year & .data$year < last_year + 1) |>
    dplyr::mutate(
      time = julian(.data$date, origin = as.Date(paste0(first_year,"-01-01")))/
        365.25 + first_year
    ) |>
    dplyr::filter(.data$time > first_year & .data$time < last_year + 1) |>
    dplyr::select("unit", "time", "cases") |>
    dplyr::arrange(.data$time, .data$unit)

  covar_list = vector("list", length(units))
  for(i in seq_along(units)){
    dmgi = dplyr::filter(data$demog, .data$unit == units[[i]])
    times = seq(from = min(dmgi$year), to = max(dmgi$year), by = 1/12)
    #times = dplyr::filter(cases_df, unit == units[[1]])$time
    switch(interp_method[[1]],
           shifted_splines = {
             pop_interp = stats::predict(
               stats::smooth.spline(x = dmgi$year, y = dmgi$pop),
               x = times
             )$y
             births_interp = stats::predict(
               stats::smooth.spline(x = dmgi$year + 0.5, y = dmgi$births),
               x = times - 4
             )$y
           },
           linear = {
             pop_interp = stats::approx(
               x = dmgi$year,
               y = dmgi$pop,
               xout = times
             )$y
             births_interp = stats::approx(
               x = dmgi$year,
               y = dmgi$births,
               xout = times - 4
             )$y
           }
    )
    covar_list[[i]] = dmgi |>
      dplyr::reframe(
        time = times,
        pop = pop_interp,
        birthrate = births_interp
      )
    covar_list[[i]] = covar_list[[i]] |>
      dplyr::mutate(
        pop_1950 = dplyr::filter(
          data$demog, .data$unit == units[[i]], .data$year == 1950
        )$pop,
        unit = units[[i]]
      )
  }
  for(i in seq_along(units)){
    log_pop_1950 = sapply(seq_along(units), function(x)
      log(covar_list[[x]][["pop_1950"]][[1]])
    )
    covar_list[[i]] = covar_list[[i]] |>
      dplyr::mutate(
        std_log_pop_1950 = (log(.data$pop_1950) - mean(log_pop_1950))/
          stats::sd(log_pop_1950)
      )
  }
  covar_df = dplyr::bind_rows(covar_list) |>
    dplyr::arrange(.data$time, .data$unit) |>
    dplyr::select("unit", "time", dplyr::everything())

  # Haversine formula for great circle distance between two points on a sphere
  # of radius r. Here, r defaults to a mean radius for the earth, in miles.
  distGreatCircle <- function(p1, p2, r = 3963.191) {
    Lon1 <- as.numeric(p1[,1])*pi/180
    Lat1 <- as.numeric(p1[,2])*pi/180
    Lon2 <- as.numeric(p2[,1])*pi/180
    Lat2 <- as.numeric(p2[,2])*pi/180
    a <- sin((Lat2-Lat1)/2)^2 + cos(Lat1)*cos(Lat2)*sin((Lon2-Lon1)/2)^2
    atan2(sqrt(a), sqrt(1 - a))*2*r
  }

  long_lat <- data$coord[,c("long","lat")]
  dmat <- matrix(0, U, U)
  for(u1 in 1:U) {
    for(u2 in 1:U) {
      dmat[u1,u2] = round(distGreatCircle(long_lat[u1,], long_lat[u2,]), 1)
    }
  }

  p <- data$demog |>
    dplyr::group_by(.data$unit) |>
    dplyr::summarize(mean_pop = mean(.data$pop)) |>
    tibble::deframe()
  v_by_g <- matrix(0, U, U)
  dist_mean <- sum(dmat)/(U*(U - 1))
  p_mean <- mean(p)
  for(u1 in 2:U){
    for(u2 in 1:(u1-1)){
      v_by_g[u1,u2] <- (dist_mean*p[u1]*p[u2]) / (dmat[u1,u2] * p_mean^2)
      v_by_g[u2,u1] <- v_by_g[u1,u2]
    }
  }
  to_C_array <- function(v) paste0("{", paste0(v, collapse = ","), "}")
  v_by_g_C_rows <- apply(v_by_g, 1, to_C_array)
  v_by_g_C_array <- to_C_array(v_by_g_C_rows)
  v_by_g_C <- pomp::Csnippet(
    paste0("const double v_by_g[",U,"][",U,"] = ",v_by_g_C_array,"; ")
  )

  set_unit_specific = pomp::Csnippet(
    paste0("const int ", model$specific_params,"_unit = 1;\n", collapse=" ")
  )
  set_shared = pomp::Csnippet(
    paste0("const int ", model$shared_params,"_unit = 0;\n", collapse=" ")
  )

  measles_globals = pomp::Csnippet(
    paste(v_by_g_C, set_unit_specific, set_shared, sep = "\n ")
  )

  # add a "1" for shared parameter names to make the pointers work
  # measles_paramnames = unlist(c(
  #   lapply(model$specific_params, function(x,U) paste0(x,1:U),U),
  #   lapply(model$shared_params, function(x) paste0(x,1))
  # ))

  first_unit_df = dplyr::filter(cases_df, .data$unit == cases_df$unit[[1]])
  model_out = spatPomp::spatPomp(
    cases_df,
    units = "unit",
    times = "time",
    t0 = with(first_unit_df, 2*time[1] - time[2]),
    unit_statenames = c('S','E','I','C'),
    covar = covar_df,
    rprocess = pomp::euler(model$rproc, delta.t = dt),
    unit_accumvars = 'C',
    paramnames =
      c(model$expanded_shared_params, model$expanded_specific_params),
    partrans = model$pt,
    globals = measles_globals,
    rinit = model$rinit,
    dmeasure = model$dmeas,
    rmeasure = model$rmeas,
    dunit_measure = model$dunit_measure
  )
  expanded_params =
    c(model$expanded_shared_params, model$expanded_specific_params)
  # TODO replace 0 with NA?
  dummy_params = rep(0, length(expanded_params))
  names(dummy_params) = expanded_params
  pomp::coef(model_out) = dummy_params
  model_out
}
