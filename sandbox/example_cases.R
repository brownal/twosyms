time=seq(0, 365, by=0.1)
env <- init_env(L=c(30, 30, 0), N=c(3.4e-7, 3.4e-7, 0), X=c(3.4e-7, 3.4e-7, 0),
	time=time)

pars <- def_pars(nsym=1)

pars$kCO2 <- 40
pars$initS <- 1e-4

run <- run_coral(env=env, pars=pars, time=time)

plot(time[-1], run$dH.Hdt[-1], type="l", xlim=c(0.1, 365), ylim=c(-0.01, 0.2),
	xlab="Time (days since bleaching)", ylab="",
	main="Coral can recover with sensitive symbiont alone")
lines(time[-1], run$S[-1]/run$H[-1], col="blue")
legend("topleft",
	c("Host growth rate", "Sensitive symbiont : host ratio"), 
	col=c("black", "blue"), lty=1)

#####

env <- init_env(L=c(34, 34, 0), N=c(3.2e-7, 3.2e-7, 0), X=c(4e-8, 4e-8, 0),
	time=time)

pars <- def_pars(nsym=1)
pars$initS <- 1e-4

pars$kCO2 <- 40

pars2 <- def_pars(nsym=2)
pars2$kCO2 <- 40
pars2$initS[1:2] <- c(0.5e-4, 0.5e-4)
pars2$astar[2] <- 1
pars2$kNPQ[2] <- 210
pars2$jSGm[2] <- 0.25
pars2$jCPm[2] <- 2.4
pars2$kROS[2] <- 200
pars2$b[2] <- 3


#485 34.17253971 3.160585e-07 3.907876e-08 40.76301 1.0181615 210.7603
#0.2472269 2.413858 198.39149 3.034511
run <- run_coral(env=env, pars=pars, time=time)
run2 <- run_coral(env=env, pars=pars2, time=time)

par(mfrow=c(2, 2))
plot(time[-1], run$dH.Hdt[-1], type="l", xlim=c(0.1, 365), ylim=c(-0.1, 0.2),
	xlab="Time (days since bleaching)", ylab="",
	main="Coral cannot recover with sensitive symbiont alone")
lines(time[-1], run$S[-1]/run$H[-1], col="blue")
legend("topleft",
	c("Host growth rate with sensitive symb.", "Sensitive symbiont : host ratio"), 
	col=c("black", "blue"), lty=1)

plot(time[-1], run2$dH.Hdt[-1], type="l", xlim=c(0.1, 365), ylim=c(-0.1, 2),
	xlab="Time (days since bleaching)", ylab="",
	main="Host recovers with tolerant symbiont")
lines(time[-1], run2$S.1[-1]/run$H[-1], col="blue")
lines(time[-1], run2$S.2[-1]/run$H[-1], col="red")
legend("topleft",
	c("Host growth rate with sensitive symb.", "Sensitive symbiont : host ratio",
		"Tolerant symbiont : host ratio"), 
	col=c("black", "blue", "red"), lty=1)

pars2 <- def_pars(nsym=2)
pars2$initS[1:2] <- c(0.5e-4, 0.5e-4)
pars2$astar[2] <- 0.85
pars2$kNPQ[2] <- 122
pars2$jSGm[2] <- 0.25
pars2$jCPm[2] <- 2
pars2$kROS[2] <- 80
pars2$b[2] <- 1.32