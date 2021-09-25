# A function that returns carbon/nitrogen limitation for host and symbiont biomass
# synthesizing units. Based on equation 23 of Cunning et al. 2017, J. Theor. Biol.
#
# Input:
#   - run, a coRal:run_coral or coRal::run_coral_oc simulation
#   - pars, the host and symbiont parameters used for run (coRal::def_pars format)
#
# Output:
#   A list whose entries are sequences of C/N limitation over time for the host ($H)
#   and symbiont ($S) or symbionts ($S.1, $S.2, ...)
#   C/N limitation returned is positive if biomass synthesis is nitrogen-limited.
#   Limitiation is negative if synthesis is carbon-limited.
#   Limitation is 0 if synthesis either co-limited or not limited by either C or N
#   (in this case, synthesis is limited by the max. synthesis rate, jHGm or jSGm).

cn_limitation <- function(run, pars) {
  # Initialize the list of limitation values to be returned
  limList <- list()

  # Get the number of symbionts
  nsym <- length(pars$initS)

  ### If there is only one symbiont
  if(nsym == 1) {

    ## Symbiont
    # inputs to symbiont growth SU
    symbC <- pars$yC * run$jCP # carbon input
    symbN <- (run$rhoN * run$H/run$S + run$rNS)/pars$nNS # nitrogen input

    # Symbiont C/N limitation (positive = N-limited, negative = C-limited)
    symbLim <- log(pmin(symbC, pars$jSGm) / pmin(symbN, pars$jSGm))
    limList$S <- symbLim # add symbiont C/N limitation to the output


    ## Host
    # inputs to the host growth SU
    hostC <- pars$yC * run$rhoC * run$S / run$H + run$jX + ifelse(is.null(run$jOC), 0, run$jOC) # carbon input
    hostN <- (run$jN + pars$nNX * run$jX + run$rNH) / pars$nNH # nitrogen input

    # Host C/N limitation
    hostLim <- log(pmin(hostC, pars$jHGm) / pmin(hostN, pars$jHGm))
    limList$H <- hostLim # add host C/N limitation to the output
  }

  ### If there are multiple symbionts
  if(nsym > 1) {

    ## Symbiont
    for(i in 1:nsym) {
      # Inputs symbiont i's growth SU
      symbC <- pars$yC * run[[paste0("jCP.", i)]] # carbon input
      symbN <- (run$rhoN * run$H/run[[paste0("S.", i)]] + run[[paste0("rNS.", i)]])/pars$nNS[i] # nitrogen input

      # Symbiont i's C/N limitation
      symbLim <- log(pmin(symbC, pars$jSGm[i]) / pmin(symbN, pars$jSGm[i]))
      limList[[paste0("S.", i)]] <- symbLim # add symb. i's C/N limitation to the output
    }

    ## Host
    # carbon from symbiont i to host
    csymbi <- sapply(1:nsym, function(i) {
        run[[paste0("rhoC.", i)]] * run[[paste0("S.", i)]]
      })

    # total carbon from symbionts to host
    csymbtot <- apply(csymbi, 1, sum)

    # carbon input to host growth SU
    hostC <- pars$yC * csymbtot / run$H + run$jX + ifelse(is.null(run$jOC), 0, run$jOC)

    # nitrogen input to host grwoth SU
    hostN <- (run$jN + pars$nNX * run$jX + run$rNH) / pars$nNH

    # Host C/N limitation
    hostLim <- log(pmin(hostC, pars$jHGm) / pmin(hostN, pars$jHGm))
    limList$H <- hostLim # add host C/N limitation to the output
  }

  return(limList)
}
