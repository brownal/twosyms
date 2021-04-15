pars_both <- def_pars(nsym=2)
pars_both$jCPm[2] <- 1.0
pars_both$kROS[2] <- 250
pars_both$jSGm[2] <- 0.15
pars_both$initS[1:2] <- c(0.5e-4, 0.5e-4)

pars_none <- def_pars(nsym=1)
pars_none$initS <- 0

time <- seq(0, 365, 0.1)

grH <- rep(NA, nrow(rescue_sims))
grHnoST <- rep(NA, nrow(rescue_sims))
cnlimH <- rep(NA, nrow(rescue_sims))
cnlimS <- rep(NA, nrow(rescue_sims))
cnlimT <- rep(NA, nrow(rescue_sims))
stratio <- rep(NA, nrow(rescue_sims))
shratio <- rep(NA, nrow(rescue_sims))

for(i in 1:nrow(rescue_sims)) {

  Li <- rescue_sims[i, "L"]
  Ni <- rescue_sims[i, "N"]
  Xi <- rescue_sims[i, "X"]

  env <- init_env(time=time, L=c(Li, Li, 0), N=c(Ni, Ni, 0), X=c(Xi, Xi, 0))

  run <- run_coral_oc(time=time, env=env, pars=pars_both)
  #run0 <- run_coral_oc(time=time, env=env, pars=pars_none)
  
  grH[[i]] <- run$dH.Hdt[length(time)]
  grHnoST[[i]] <- run0$dH.Hdt[length(time)]

  cnlim <- cn_limitation(run, pars_both)
  cnlimH[[i]] <- cnlim$H[length(time)]
  cnlimS[[i]] <- cnlim$S.1[length(time)]
  cnlimT[[i]] <- cnlim$S.2[length(time)]

  stratio[[i]] <- (run$S.1/run$S.2)[length(time)]
  shratio[[i]] <- ((run$S.1 + run$S.2)/run$H)[length(time)]
}

summary(grHnoST)

summary(stratio)

hist(cnlimH, freq=FALSE)

sum(stratio > 0.99)/length(stratio)

plot(cnlimS[stratio > 0.99], grH[stratio > 0.99])

plot(cnlimS[stratio > 0.99], shratio[stratio > 0.99])
summary(stratio[stratio <= 0.99])

summary(shratio)

plot(cnlimS, cnlimT)

plot(cnlimS, shratio, xlab="C-N limitation of symbiont (positive = N-lim)",
	ylab="Total symbiont:host ratio", col="#648FFF80", pch=16)

points(cnlimT, shratio, col="#FFB00080", pch=15)
legend("topright", legend=c("Sensitive", "Tolerant"), pch=c(16, 15),
   col=c("#648FFF", "#FFB000"))

grHist <- hist(grHnoST, plot=FALSE, breaks=seq(-0.015, 0.005, 0.0025),
	right=FALSE)
grHist$density <- grHist$counts/sum(grHist$counts)


plot(grHist, freq=FALSE, ylab="Fraction of simulations",
  main="Aposymbiotic host growth rate after 365 days", col="gray", 
  xlab="Bars include only values < rightmost value")


sum(cnlimS >= 0)/140


plot(cnlimS[shratio < 0.2], shratio[shratio < 0.2])

## This probably due to transient dynamics; if I limit to simulations where
# the total symbiont:host ratio < 0.2, none of the sensitive symbionts are
# C-limited. (Tolerant symbionts are never N-limited in this set of sims)


## What if I run it for longer/shorter?

## The problem is that we thought the switch to N-limitation should happen
# immediately upon return to positive growth --> Are we just wrong about the
# mechanism? Or are we wrong about something more important?


# Also check recovery with just S




## What's the N-limitation when S:H becomes 0.2?
time <- seq(0, 365*2, 0.1)

for(i in 1:nrow(rescue_sims)) {

  Li <- 4 #rescue_sims[i, "L"]
  Ni <- rescue_sims[i, "N"]
  Xi <- rescue_sims[i, "X"]

  env <- init_env(time=time, L=c(Li, Li, 0), N=c(Ni, Ni, 0), X=c(Xi, Xi, 0))

  run <- run_coral_oc(time=time, env=env, pars=pars_both)

  cnlim <- cn_limitation(run, pars_both)

  sh <- ((run$S.1 + run$S.2)/run$H)


  cnlimH[[i]] <- cnlim$H[length(time)]
  cnlimS[[i]] <- cnlim$S.1[length(time)]
  cnlimT[[i]] <- cnlim$S.2[length(time)]

  stratio[[i]] <- (run$S.1/run$S.2)[length(time)]
  shratio[[i]] <- ((run$S.1 + run$S.2)/run$H)[length(time)]
}