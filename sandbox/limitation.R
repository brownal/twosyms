
# Calculate the degree to which product formation is limited by each of its substrates
# Equation 23 in Cunning et al. 2017
# Input:
#     - jS1: Input flux of substrate 1
#     - jS2: Input flux of substrate 2
#     - jPm: Maximum rate of product formation
# Output:
#     - Negative if substrate 1 is "more limiting"
#     - Positive if substrate 2 is "more limiting"
#     - 0 if neither substrate is limiting (or they are equally limiting)
limitation <- function(jS1, jS2, jPm) {
      log(min(jS1, jPm)/min(jS2, jPm))
}


# Carbon input flux to host growth SU
# There is a 1-symbion and 2-symbiont version of this equation
# Input:
#     - yC: yield of biomass from organic carbon
#     - rhoC or rhoC.1, rhoC.2: fixed carbon shared by the symbiont
#     - S or S.1, S.2: symbiont biomass
#     - H: host biomass
#     - jX: Input of carbon from the coral feeding
# Output:
# -
carbonInput1S <- function(yC, rhoC, S, H, jX) {
      yC * rhoC * S/H + jX
}

carbonInput2S <- function(yC, rhoC.1, rhoC.2, S.1, S.2, H, jX) {
      yC * (rhoC.1 * S.1 + rhoC.2 * S.2)/H + jX
}


# Nitrogen input flux to host growth SU
# Input:
#     - jN: Input of environmental nitrogen
#     - nNX: Nitrogen per carbon in food
#     - jX: Input of carbon from coral feeding
#     - rNH: Recycled nitrogen from host turnover
#     - nNH: N:C molar ratio in host biomass
nitrogenInput <- function(jN, nNX, jX, rNH, nNH) {
      (jN + nNX * jX + rNH) * (1/nNH)
}




time <- seq(1,3000,0.1)
env <- init_env(time=time, X=c(5e-8,5e-8,0), N=c(5e-7,5e-7,0), L=c(15,15,0))

pars <- def_pars(nsym=2)

# Make symbiont 2 have lower jCPm and higher kROS (i.e., tolerant symbiont)
pars$jCPm[2] <- 1.0
pars$kROS[2] <- 250
pars$jSGm[2] <- 0.15
pars$initS[1:2] <- 0.5e-4

run <- run_coral(time=time, env=env,pars=pars)


limvstime <- mapply(function(H, S.1, S.2, jX, jN, rhoC.1, rhoC.2, rNH) {
            limitation(carbonInput(pars$yC, rhoC.1, rhoC.2, S.1, S.2, H, jX),
                  nitrogenInput(jN, pars$nNX, jX, rNH, pars$nNH), pars$Hgm)
            },
      run$H, run$S.1, run$S.2, run$jX, run$jN, run$rhoC.1, run$rhoC.2, run$rNH)

limvstime <- mapply(function(H, S.1, S.2, jX, jN, rhoC.1, rhoC.2, rNH, yC, nNX, nNH, Hgm) {
            limitation(carbonInput(yC, rhoC.1, rhoC.2, S.1, S.2, H, jX),
                  nitrogenInput(jN, nNX, jX, rNH, nNH), Hgm)
            },
      H=run$H, S.1=run$S.1, S.2=run$S.2, jX=run$jX, jN=run$jN, rhoC.1=run$rhoC.1, rhoC.2=run$rhoC.2, rNH=run$rNH,
      MoreArgs=list(yC=pars$yC, nNX=pars$nNX, nNH=pars$nNH, Hgm=pars$Hgm))


par(mfrow=c(2,1))
plot(time, run$S.1/run$H, type="l", col="blue", ylim=c(-3, 3), xlim=c(0, 800))
lines(time, run$S.2/run$H, col="red")
lines(time, 10*run$dH.Hdt, col="black")
lines(time, limvstime, col="green")

plot(time, run$S.1/run$H, type="l", col="blue", ylim=c(-3, 3), xlim=c(680, 720))
lines(time, run$S.2/run$H, col="red")
lines(time, 10*run$dH.Hdt, col="black")
lines(time, limvstime, col="green")
