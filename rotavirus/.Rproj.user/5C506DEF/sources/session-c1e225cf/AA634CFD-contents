library(doParallel)


cores <- 10
registerDoParallel(cores)

system.time(
  rnorm(10^8)
) -> time0

system.time(
  foreach(i=1:10^4) %dopar% rnorm(10^4)
  foreach(i=1:10) %dopar% rnorm(10^2)
) -> time4