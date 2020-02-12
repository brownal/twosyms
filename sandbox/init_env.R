#' Initialize abiotic environment.
#'
#' This function returns a list of vectors that define the external irradiance
#' (L), DIN concentration (N), and prey concentration (X) to be used as
#' subsequent inputs for a model simulation.
#'
#' The options for the functional form of each environmental factor are:
#' \itemize{
#'   \item 0=linear increase from min to max
#'   \item 1=linear decrease from max to min
#'   \item 2=sinusoid with period of one year and a range from min to max
#' }
#'
#' For light only, the following additional options are available:
#' \itemize{
#'   \item 3=sinusoid as in option 2 except that light is decreasing at time = 0
#'   \item 4=sinusoid with a higher peak the first year
#'   \item 5=sinusoid with a higher peak the first year (Either 1.5x the usual amplitude or L[4] if given)
#' }
#'
#' For prey only, functional form 2 will produce a sinusoid with
#' minimum = min + 10^-6 * min, maximum = max + 10^-6 * min
#' 
#' @param time A vector of time steps at which the model should be evaluated
#'   (units=days) (e.g., seq(0, 365, 0.1))
#' @param L A vector with length 3 defining external light (min, max, functional
#'   form (see details)). Units=mol photons m^-2 d^-1.
#' @param N A vector with length 3 defining external DIN concentration (min,
#'   max, functional form (see details)). Units=mol N L^-1.
#' @param X A vector with length 3 defining external prey concentration (min,
#'   max, functional form (see details)). Units=C-mol X L^-1.
#' @return A list of 3 numeric vectors (each with length=length(time))
#'   corresponding to light, DIN, and prey environmental forcing functions.
#' @examples
#' env1 <- init_env(time=seq(1,365,0.1), L=c(20,40,1), N=c(1e-7,1e-7,0), X=c(0,0,0))

init_env <- function(time, L, N, X) {

  # Light
  L <- if (L[3]==0) {
    seq(L[1], L[2], along=time)
  } else if (L[3]==1) {
    seq(L[2], L[1], along=time)
  } else if (L[3]==2) {
    0.5 * (L[2] - L[1]) * sin(0.0172*time) + L[1] + 0.5 * (L[2] - L[1])
  } else if (L[3]==3) {
    0.5 * (L[2] - L[1]) * sin(0.0172*(time-182)) + L[1] + 0.5 * (L[2] - L[1])
  } else if (L[3]==4) {
    f <- 0.5 * (L[2] - L[1]) * sin(0.0172*time)
    f <- scales::rescale(f, to=c(50, L[2]))
    ff <- scales::rescale(1 + (0.01 / (1 + exp(0.05*(time-200)))), to=c(0.6,1))
    f * ff
  } else if (L[3]==5) {
      # If there is an extra entry in L for the max light level the first year, use that
      # otherwise use 1.5*(L[2]-L[1]) + L[1] for the maximum light in the first year
      if (length(L) > 3) {
        year1max <- L[4]
      } else {
        year1max <- 1.5 * (L[2] - L[1]) + L[1]
      }
      
      # Check if input does not meet conditions required for form 5. If so, print a warning
      # and use functional form 3 instead.
    if ((L[1] > L[2]) || (year1max < L[2])) {
      if (L[1] > L[2]) {
        warning("Min light > max light. Using functional form 2 for light instead.")
      } else if (year1max < L[2]) {
        warning("Max light in year 1 < max light in other years. Using functional form 2 for light instead.")
      }
      0.5 * (L[2] - L[1]) * sin(0.0172*time) + L[1] + 0.5 * (L[2] - L[1])
    } else {
        basicCycle <- (sin(time*2*pi/365) + 1)/2
        ifelse(abs(time-365/4) < (365/2), (year1max - L[1])*basicCycle + L[1], (L[2] - L[1])*basicCycle + L[1])
    }
  }

  # DIN
  N <- if (N[3]==0) {
    seq(N[1], N[2], along=time)
  } else if (N[3]==1) {
    seq(N[2], N[1], along=time)
  } else {
    0.5 * (N[2] - N[1]) * sin(0.0172*time) + N[1] + 0.5 * (N[2] - N[1])
  }

  # Prey
  X <- if (X[3]==0) {
    seq(X[1], X[2], along=time)
  } else if (X[3]==1) {
    seq(X[2], X[1], along=time)
  } else {
    0.5 * (X[2] - X[1]) * sin(0.0172*time) + X[1]*10^-6 + 0.5 * (X[2] - X[1])
  }
  # Set environment specifications
  env <- list(L=L, N=N, X=X)
  return(env)
}