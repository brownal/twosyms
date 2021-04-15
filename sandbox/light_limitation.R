# A function that returns CO2/light limitation for host and symbiont biomass
# synthesizing units. Based on equation 23 of Cunning et al. 2017, J. Theor. Biol.
#
# Input:
#   - run, a coRal:run_coral or coRal::run_coral_oc simulation
#   - pars, the host and symbiont parameters used for run (coRal::def_pars format)
#
# Output:
#   A list whose entries are sequences of CO2/light limitation over time for the
#   symbiont ($S) or symbionts ($S.1, $S.2, ...)
#   CO2/light limitation returned is positive if photosynthesis is light-limited.
#   Limitiation is negative if photosynthesis is CO2-limited.
#   Limitation is 0 if synthesis either co-limited or not limited by either light
#   or CO2 (in this case, synthesis is limited by the max. photosynthesis rate, jCPm).

light_limitation <- function(run, pars) {
  # Initialize the list of limitation values to be returned
  limList <- list()

  # Get the number of symbionts
  nsym <- length(pars$initS)

  ### If there is only one symbiont
  if(nsym == 1) {

    # Inputs to photosynthesis SU
    light <- pars$yCL * run$jL # light input
    co2 <- (run$jCO2 + run$rCH) * (run$H / run$S) + run$rCS  # CO2 input

    # CO2/light limitation (positive = light-limited, negative = CO2-limited)
    clLim <- log(pmin(co2, pars$jCPm) / pmin(light, pars$jCPm))
    limList$S <- clLim # add CO2/light limitation to the output
  }

  ### If there are multiple symbionts
  if(nsym > 1) {

    for(i in 1:nsym) {
      # Inputs to symbiont i's photosynthesis SU
      light <- pars$yCL[[i]] * run[[paste0("jL.", i)]] # light input
      co2 <-  (run$jCO2 + run$rCH) * (run$H / run[[paste0("S.", i)]]) + run[[paste0("rCS.", i)]] # CO2 input

      # Symbiont i's CO2/light limitation (positive = light-limited, negative = CO2-limited)
      clLim <- log(pmin(co2, pars$jCPm[i]) / pmin(light, pars$jCPm[i]))
      limList[[paste0("S.", i)]] <- clLim # add symb. i's CO2/light limitation to the output
    }
  }

  return(limList)
}
