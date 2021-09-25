#' Run simulation of coral-\emph{Symbiodinium} model.
#'
#' This function runs a simulation of the coral-Symbiodinium bioenergetic model
#' using a specified time vector, environmental inputs, and parameters.
#'
#' @param time A vector of time steps at which the model should be evaluated
#'   (units=days) (e.g., seq(0, 365, 0.1))
#' @param env A list defining the constant environment in which to run model to steady state, with three named elements:
#' \describe{
#'     \item{L}{External light in mol photons m^-2 d^-1 (numeric).}
#'     \item{N}{External [DIN] in mol N L^-1 (numeric).}
#'     \item{X}{External [prey] in C-mol X L^-1 (numeric).}
#'   }
#' @param pars A list of (named) parameter values to use in the simulation.
#'   Parameter names must match those returned by \code{def_pars()}.
#' @return A list with the following named elements:
#'   \describe{
#'     \item{time}{The time vector used in the simulation}
#'     \item{env}{The environment object used in the simulation}
#'     \item{pars}{The parameter values used in the simulation}
#'     \item{H}{A tibble of all host biomass-specific model
#'     fluxes at each time step}
#'     \item{S}{A tibble of all symbiont biomass-specific
#'     model fluxes at each time step}
#'   }
#' @seealso \code{\link{init_env}}, \code{\link{def_pars}}
#' @examples
#' time1 <- seq(0, 365, 0.1)
#' pars1 <- def_pars()
#' env1 <- init_env(time=time1, L=c(20,40,1), N=c(1e-7,1e-7,0), X=c(0,0,0))
#' run1 <- run_coral(time=time1, env=env1, pars=pars1)

run_coral_comp_exclusion <- function(env, pars, max_time, dt, threshold) {
  vsynth <- Vectorize(synth)

  # Get number of symbionts in run
  nsym <- length(pars$initS)

  # Check that there are exactly 2 symbionts
  if(nsym != 2) {
    stop("Pars must have two symbionts.")
  }

  # Check that the threshold condition is not already met
  if(pars$initS[1]/pars$initS[2] < threshold){
    stop("S2/S1 ratio initialized below threshold.")
  }

  # Set up time to be the maximum time
  time=seq(0, max_time, dt)

  # Set up the environment over time
  envVtime <- coRal::init_env(time=time, L=c(env$L, env$L, 0), N=c(env$N, env$N, 0), X=c(env$X, env$X, 0))

  # Set initial values
  # ==================
  # Host fluxes
  for(x in c("jX", "jN", "rNH", "rhoN", "jeC", "jCO2", "jHG", "rCH", "dH.Hdt", "H")) {
    assign(x, rep(NA, length(time)))
  }
  jX <- (pars$jXm * envVtime$X / (envVtime$X + pars$KX))
  jN <- (pars$jNm * envVtime$N / (envVtime$N + pars$KN))
  rNH <- rep(pars$jHT0 * pars$nNH * pars$sigmaNH, length(time))
  rhoN[1] <- jN[1]
  jeC[1] <- 10
  jCO2[1] <- pars$kCO2 * jeC[1]
  jHG[1] <- 0.25
  #jHT=pars$jHT0
  rCH[1] <- pars$jHT0 * pars$sigmaCH
  dH.Hdt[1] <- pars$jHGm
  H[1] <- pars$initH
  # Symbiont fluxes
  for (x in c("rNS", "jL", "jCP", "jeL", "jNPQ", "jSG", "rhoC", "jST", "rCS", "cROS", "dS.Sdt", "S")) {
    assign(x, matrix(NA, ncol=nsym, nrow=length(time)))
  }
  rNS <- matrix(pars$jST0 * pars$nNS * pars$sigmaNS, ncol=nsym, nrow=length(time))
  jL[1,] <- envVtime$L[1] * pars$astar
  jCP[1,] <- pmax(0, vsynth(jL[1,] * pars$yCL, jCO2[1]*H[1]/pars$initS, pars$jCPm), na.rm=T)
  jeL[1,] <- pmax(jL[1,] - jCP[1,]/pars$yCL, 0)
  jNPQ[1,] <- pars$kNPQ
  #jCO2w=H$jCO2*H$H/S - jCP,
  jSG[1,] <- pars$jSGm/10
  rhoC[1,] <- jCP[1,]
  #jNw=0,
  jST[1,] <- pars$jST0
  rCS[1,] <- pars$jST0 * pars$sigmaCS
  cROS[1,] <- 1
  dS.Sdt[1,] <- pars$jSGm
  S[1,] <- pars$initS

  # Run simulation by updating
  # ==========================
  # Set the counter to know how many entries of the output lists/matrices we've filled
  counter <- 1

  for (t in 2:length(time)) {
    counter <- counter + 1

    # Symbiont fluxes
    S.t <- sum(S[t-1,])  # Get total symbiont abundance from prev time step
    # Photosynthesis
    # ==============
    # Light input flux
    jL[t,] <- (1.256307 + 1.385969 * exp(-6.479055 * (S.t/H[t-1]))) * envVtime$L[t] * pars$astar
    # CO2 input flux
    rCS[t,] <- pars$sigmaCS * (pars$jST0 + (1-pars$yC)*jSG[t-1,]/pars$yC)  # metabolic CO2 recycled from symbiont biomass turnover
    # Production flux (photosynthetic carbon fixation)
    jCP[t,] <- vsynth(jL[t,] * pars$yCL, (jCO2[t-1] + rCH[t-1])*H[t-1]/S.t + rCS[t,], pars$jCPm) / cROS[t-1,]
    # Rejection flux: CO2 (wasted to the envVtimemeironment)
    # jCO2w[t] <- max((H$jCO2[t-1] + H$rCH[t-1])*H$H[t-1]/S.t + rCS[t] - jCP[t], 0)
    # Rejection flux: excess light energy not quenched by carbon fixation
    jeL[t,] <- pmax(jL[t,] - jCP[t,]/pars$yCL, 0)
    # Amount of excess light energy quenched by NPQ
    jNPQ[t,] <- (pars$kNPQ^(-1)+jeL[t,]^(-1))^(-1/1)  # single substrate SU
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    cROS[t,] <- 1 + (pmax(jeL[t,] - jNPQ[t,], 0) / pars$kROS)^pars$k
    # Symbiont biomass
    # ================
    # Nitrogen input flux
    # rNS[t] <- pars$jST0[i] * pars$nNS[i] * pars$sigmaNS[i]  # Recylced N from symbiont biomass turnover.
    # Production flux (symbiont biomass formation)
    jSG[t,] <- vsynth(pars$yC*jCP[t,], (rhoN[t-1]*H[t-1]/S.t + rNS[t,])/pars$nNS, pars$jSGm)
    # Rejection flux: carbon (surplus carbon shared with the host)
    rhoC[t,] <- pmax(jCP[t,] - jSG[t,]/pars$yC, 0)
    # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
    # jNw[t] <- max(H$rhoN[t-1]*H$H[t-1]/S.t + rNS[t] - pars$nNS[i] * jSG[t], 0)
    # Symbiont biomass loss (turnover)
    jST[t,] <- pars$jST0 * (1 + pars$b * (cROS[t,] - 1))
    # State equations
    dS.Sdt[t,] <- jSG[t,] - jST[t,]  # Specific growth rates (Cmol/Cmol/d)
    S[t,] <- S[t-1,] + dS.Sdt[t,] * S[t-1,] * dt  # Biomass (Cmol)

    # Total amount of carbon shared by all symbionts
    rhoC.t <- sum(rhoC[t,]*S[t-1,])

    # Host fluxes
    # ============
    # Food input flux (prey=both carbon and nitrogen)
    # jX[t] <- (pars$jXm * env$X[t] / (env$X[t] + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    # jN[t] <- (pars$jNm * env$N[t] / (env$N[t] + pars$KN))  # N uptake from the environment
    # rNH[t] <- pars$jHT0 * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Production flux (host biomass formation)
    jHG[t] <- synth(pars$yC*(rhoC.t/H[t-1] + jX[t]), (jN[t] + pars$nNX*jX[t] + rNH[t]) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    rhoN[t] <- max(jN[t] + pars$nNX * jX[t] + rNH[t] - pars$nNH * jHG[t], 0)
    # Rejection flux: carbon -- given back to symbiont as CO2 input to photosynthesis
    jeC[t] <- max(jX[t] + rhoC.t/H[t-1] - jHG[t]/pars$yC, 0)
    # Host biomass loss
    # jHT[t] <- pars$jHT0
    rCH[t] <- pars$sigmaCH * (pars$jHT0 + (1-pars$yC)*jHG[t]/pars$yC)  # metabolic CO2 recycled from host biomass turnover
    jCO2[t] <- pars$kCO2 * jeC[t]  # carbon not used in host biomass is used to activate CCM's that deliver CO2 to photosynthesis
    # State equations
    dH.Hdt[t] <- jHG[t] - pars$jHT0  # Specific growth rates (Cmol/Cmol/d)
    H[t] <- H[t-1] + dH.Hdt[t] * H[t-1] * dt  # Biomass (Cmol)

    # Stop the simulation if the threshold has been reached
    if(S[t,2]/S[t,1] <= threshold) {
      break
    }
  }

  # Return results, trimming timepoints not simulated
  # ==============
  out <- data.frame(time=time[1:counter], env.L=envVtime$L[1:counter], env.N=envVtime$N[1:counter], env.X=envVtime$X[1:counter],
                    jX=jX[1:counter], jN=jN[1:counter], rNH=rNH[1:counter], rhoN=rhoN[1:counter], jeC=jeC[1:counter],
                    jCO2=jCO2[1:counter], jHG=jHG[1:counter], rCH=rCH[1:counter], dH.Hdt=dH.Hdt[1:counter], H=H[1:counter],
                    rNS=rNS[1:counter,], jL=jL[1:counter,], jCP=jCP[1:counter,], jeL=jeL[1:counter,], jNPQ=jNPQ[1:counter,],
                    jSG=jSG[1:counter,], rhoC=rhoC[1:counter,], jST=jST[1:counter,], rCS=rCS[1:counter,],
                    cROS=cROS[1:counter,], dS.Sdt=dS.Sdt[1:counter,], S=S[1:counter,])
  return(list(fluxes=out, threshold=threshold, reachedThreshold=(S[t,2]/S[t,1] <= threshold)))
}
