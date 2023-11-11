Sys.time()
devtools::session_info()
######### load-packages ###############################
# library(measlespkg)
library(foreach)
library(pomp)
source("model_mechanics_001_stocks.R")
source("sample_initial_pparams_ul.R")
source("make_measlesPomp.R")
source("run_round.R")
source("make_rw_sd.R")
source("combine_top_fits.R")
source("grab_top_fits.R")
source("plot_traces.R")
source("sim_plots.R")
source("panel_mechanics.R")
source("coef_to_pparams.R")
source("choose_units.R")
######## Source functions ############################
# invisible(sapply(list.files(path = "./R/functions", pattern = "*.R"),
#                  function(x) source(paste0("./R/functions/", x))))

######## Get arguments from command line #############
(out_dir = as.character(Sys.getenv("out_dir", unset = NA)))
## ############# OPTIONS #############################
# Set number of cores
ncores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE", unset = NA))
if(is.na(ncores)) ncores = parallel::detectCores()/4
print(ncores)

# Set fitting and filter parameters
RUN_LEVEL = 2
NP_MIF       = switch(RUN_LEVEL, 4, 5000)
NMIF         = switch(RUN_LEVEL, 4,  100)
NREPS_MIF    = switch(RUN_LEVEL, ncores,  ncores)
NP_EVAL      = switch(RUN_LEVEL, 4, 10000)
NREPS_EVAL   = switch(RUN_LEVEL, ncores,  ncores)
# TOP_N_FITS selects top fits from likelihood evaluation file specified in
# PREVIOUS_FIT_PATH. TOP_N_FITS must divide NREPS_MIF.
TOP_N_FITS   = switch(RUN_LEVEL, 2,  12)
DATA = ur_measles
MODEL = model_mechanics_001_stocks()
UNITS = c("Stocks")
#unique(twentycities$measles$unit)
BLOCK_MIF2 = TRUE
INTERP_METHOD = "shifted_splines"
# SIM_MODEL specifies whether simulated data from a given model should be used.
SIM_MODEL = NULL
COOLING_FRAC = 0.5
# EVAL_POINTS sets points to evaluate EVAL_PARAM at when performing profile
# search. Set both equal to NULL to do regular search.
EVAL_POINTS = NULL
EVAL_PARAM = NULL
### Seeds ###
MAIN_SEED = 169566674

(array_job_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)))

# Add to MAIN_SEED if running array job
if(!is.na(array_job_id)){
  MAIN_SEED = MAIN_SEED + array_job_id
  print(MAIN_SEED)
}
# Set PREVIOUS_FIT_PATH to NULL to choose starting parameters from a box
# instead of from a previous fit. Setting this equal to a path will nullify
# the portion of code which chooses starting parameters from a box.
PREVIOUS_FIT_PATH = NULL

# Use INITIAL_RW_SD to set random walk standard deviations for parameters.
DEFAULT_SD = 0.02
IVP_DEFAULT_SD = DEFAULT_SD*3
INITIAL_RW_SD = c(
  S_0 = IVP_DEFAULT_SD,
  E_0 = IVP_DEFAULT_SD,
  I_0 = IVP_DEFAULT_SD,
  R_0 = IVP_DEFAULT_SD,
  R0 = DEFAULT_SD,
  sigmaSE = DEFAULT_SD,
  amplitude = DEFAULT_SD,
  rho = DEFAULT_SD,
  gamma = DEFAULT_SD,
  psi = DEFAULT_SD,
  iota = DEFAULT_SD,
  sigma = DEFAULT_SD,
  cohort = DEFAULT_SD,
  alpha = DEFAULT_SD*10^(-2),
  mu = 0
)
if(!is.null(EVAL_PARAM))
  INITIAL_RW_SD[[EVAL_PARAM]] = 0

# Specify names of output files
RESULTS_FILE = "fit_results_out.rds"

### Diagnostic parameters ###
# For plots and final loglik calc, use combination of parameters which yields
# best sum of unit loglik?
USE_BEST_COMBO = TRUE
PLOT_TRACES = TRUE
PLOT_SIMS = FALSE
################## SETUP ###########################################
set.seed(MAIN_SEED)
# Create directory for output if it does not exist
write_path = switch(
  RUN_LEVEL,
  "./output/DUMMY/",
  out_dir
)
if(!dir.exists(write_path)) dir.create(write_path)
write_results_to = paste0(write_path, RESULTS_FILE)

# Use observations from simulation?
if(!is.null(SIM_MODEL)){
  sim = simulate(
    SIM_MODEL,
    seed = SIM_MODEL_SEED
  )
  sim_obs_list = lapply(seq_along(sim), function(x){
    pomp::obs(sim[[x]]) |>
      as.numeric()
  })
} else {
  sim_obs_list = NULL
}

###### Starting parameters #############################

# Specify box to sample initial specific parameters from
specific_radii = tibble::tribble(
  ~param, ~radius,
  "beta1",     2.5,
  "beta2",     0.1,
  "beta3",     0.1,
  "beta11",    0.025,
  "phi",       0.15,
  "od",        0.15,
  "sigma",     0.1,
)
if(!is.null(EVAL_PARAM))
  specific_radii[specific_radii[["param"]] == EVAL_PARAM, "radius"] = 0

specific_bounds = tibble::tribble(
  ~param,       ~lower,        ~upper,
  "beta1",        10,            15,
  "beta2",        0.2,           0.4,
  "beta3",        0.3,           0.5,
  "beta11",       0.11,          0.16,
  "phi",          0.01,          0.3,
  "od",           0.001,         0.3,
  "sigma",        0.001,         0.2,
)

if(!is.null(EVAL_PARAM)){
  eval_param_rows = specific_radii[["param"]] == EVAL_PARAM
  specific_bounds[eval_param_rows, "lower"] = EVAL_POINTS[[array_job_id]]
  specific_bounds[eval_param_rows, "upper"] = EVAL_POINTS[[array_job_id]]
}

# Sample initial parameters and place into lists
initial_pparams_list = sample_initial_pparams_ul(
  specific_box_specs = specific_bounds,
  units = UNITS,
  n_draws = NREPS_MIF
)

# Get starting parameters from previous fit
if(!is.null(PREVIOUS_FIT_PATH)){
  fit_results_in = readRDS(PREVIOUS_FIT_PATH)
  EL_in = fit_results_in$EL_out[[length(fit_results_in$EL_out)]]
  initial_pparams_list = duplicate_top_pparams(
    EL_in,
    out_length = NREPS_MIF,
    top_n = TOP_N_FITS,
    combine = TRUE
  )
}

################## Construct panelPomp object ##########################
measlesPomp_mod = make_measlesPomp(
  DATA |> choose_units(UNITS),
  starting_pparams = initial_pparams_list[[1]],
  model = MODEL,
  interp_method = INTERP_METHOD
)

if(!is.null(EVAL_POINTS)){
  coef_names = paste0(EVAL_PARAM, "[", UNITS, "]")
  coef(measles_ppomp_mod)[coef_names] = EVAL_POINTS[[array_job_id]]
}
###### MODEL FITTING #####################################
round_out = run_round(
  measlesPomp_mod,
  initial_pparams_list = initial_pparams_list,
  rw_sd = make_rw_sd(INITIAL_RW_SD),
  cooling_frac = COOLING_FRAC,
  nmif = NMIF,
  np_mif2 = NP_MIF,
  np_eval = NP_EVAL,
  nreps_eval = NREPS_EVAL,
  block = BLOCK_MIF2,
  ncores = ncores,
  write_results_to = write_results_to,
  print_times = TRUE
)


EL_final = round_out$EL_out[[length(round_out$EL_out)]]
print(as.data.frame(dplyr::arrange(EL_final$fits[,1:2], dplyr::desc(logLik))))

# Evaluate at parameters of best ULL combination
if(USE_BEST_COMBO){
  tictoc::tic()
  top_params = combine_top_fits(EL_final)$fits[-(1:2)]
  eval_model = measlesPomp_mod
  panelPomp::coef(eval_model) = top_params
  EL_out_best = pomp::bake(file = paste0(write_path, "best_eval.rds"),
    eval_logLik(
    model_obj_list = list(eval_model),
      ncores = ncores,
      np_pf = NP_EVAL,
      nreps = ncores*8,
      seed = NULL,
      divisor = 8164
    )
  )
  tictoc::toc()
}
EL_out_best$fits[,1:2]

################ diagnostics ###############################################
if(USE_BEST_COMBO == FALSE){
  top_params = grab_top_fits(EL_final)$fits %>%
    dplyr::select(-logLik, -se) %>%
    unlist()
}

if(PLOT_TRACES){
  plot_traces(
    round_out$mif2_out,
    plot_shared = "loglik",
    plot_specific = ".ALL",
    print_plots = FALSE,
    save_dir = paste0(write_path, "trace_plots/")
  )
}

if(PLOT_SIMS){
  measlesPomp_sim = measlesPomp_mod
  coef(measlesPomp_sim) = top_params
  AK_mod = AK_model()
  sim_plots(
    true_model = AK_mod,
    sim_model = measlesPomp_sim,
    print_plots = FALSE,
    save_dir = paste0(write_path, "sim_plots/")
  )
}



