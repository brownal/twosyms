library("coRal")
source("C://Users/Alexandra/Documents/coRal2/R/run_coral.R")

# shock takes a simulation of a healthy host and introduces a brief environmental shock, followed by recovery
# Input
# - healthy_run: a run_coral simulation
# - pars: the parameters used in healthy_run
# - shock_env: the light, nitrogen and food levels during the environmental shock; list of the form (L=light, N=nitrogen, X=food)
# - recovery_env: light, nitrogen, and food levels after the environmental shock; same form as shock_env
# - start_time: the time in healthy_run at which the environmental shock should occur
# - shock_duration: how long the environmental shock lasts for
# - recovery_duration: how long the host should spend in the recovery environment post-shock
# - dt: time step for run_coral
# Output
# - a run_coral simulation that beings at start_time + dt and simulates shock and recovery
shock <- function(healthy_run, pars, shock_env, recovery_env, start_time, shock_duration, recovery_duration, dt=0.1) {
	# Time to run simulation over
	time <- seq(start_time + dt, start_time + shock_duration + recovery_duration, dt)

	# Environmental conditions
	env <- list(
		L=c(rep(shock_env$L, length(seq(dt, shock_duration, dt))), rep(recovery_env$L, length(seq(dt, recovery_duration, dt)))),
		N=c(rep(shock_env$N, length(seq(dt, shock_duration, dt))), rep(recovery_env$N, length(seq(dt, recovery_duration, dt)))),
		X=c(rep(shock_env$X, length(seq(dt, shock_duration, dt))), rep(recovery_env$X, length(seq(dt, recovery_duration, dt)))))

# Get initial fluxes
	pre_fluxes <- lapply(healthy_run, function(x) {x[[start_time/dt]]})

	# Initial symbiont and host biomasses
	pre_pars <- pars
	pre_pars$initH <- healthy_run$H[[start_time/dt]]
	pre_pars$initS <- unlist(lapply(healthy_run[grep(paste0("^S", "\\.[[:digit:]]*"), names(healthy_run))],
		function(x) { x[[start_time/dt]] }))

	# Run simulation
	run <- run_coral(time=time, env=env, pars=pre_pars,
		fluxes=pre_fluxes)

	return(run)
}


# time_t.s finds the time that the tolerant:sensitive ratio first drops below a value of interest
time_t.s <- function(run, ratio) {
	run$time[diff(run$S.2/run$S.1 > ratio) != 0]
}

# time_t.h finds the time that the tolerant:host ratio first drops below a value of interest
time_t.h <- function(run, ratio) {
	run$time[diff(run$S.2/run$H > ratio) != 0]
}

# plot_symbs makes plots of S:H ratios for simulations with two symbionts
# Negative host growth shown with dashed line for S:H ratio
plot_symbs <- function(run, col1, col2, new_plot=FALSE, ...) {
if(new_plot){
plot(run$time, run$S.1/run$H, lwd=2, lty=2, col=col1, type="l",
	xlab="Time (days)", ylab="Symbiont:host ratio", ...)
} else {
lines(run$time, run$S.1/run$H, lwd=2, lty=2, col=col1)
}

lines(run$time[(run$dH.Hdt > 0)],
	run$S.1[(run$dH.Hdt > 0)]/run$H[(run$dH.Hdt > 0)], lwd=2, col=col1)

lines(run$time, run$S.2/run$H, lwd=2, lty=2, col=col2)
lines(run$time[(run$dH.Hdt > 0)],
	run$S.2[(run$dH.Hdt > 0)]/run$H[(run$dH.Hdt > 0)], lwd=2, col=col2)
}



################
# Know rescue cases: let's find a couple!
# 







################ A whole lot of runs and plots

pars <- def_pars(nsym=2)
pars$kROS[[2]] <- 250
pars$jCPm[[2]] <- 1
pars$jSGm[[2]] <- 0.15
pars$initS <- c(0.25, 0.75)



time <- seq(0, 300, 0.1)
env <- list(L=rep(20, length(time)), N=rep(2e-6, length(time)), X=rep(2e-7, length(time)))

time_preh0.1 <- seq(0, time_h0.1, 0.1)
env_preh0.1 <- list(L=rep(20, length(time_preh0.1)),
	N=rep(2e-6, length(time_preh0.1)), X=rep(2e-7, length(time_preh0.1)))

time_preh0.25 <- seq(0, time_h0.25, 0.1)
env_preh0.25 <- list(L=rep(20, length(time_preh0.25)),
	N=rep(2e-6, length(time_preh0.25)), X=rep(2e-7, length(time_preh0.25)))

run_preh0.1 <- run_coral(time=time_preh0.1, env=env_preh0.1, pars=pars)
run_preh0.25 <- run_coral(time=time_preh0.25, env=env_preh0.25, pars=pars)

healthy_run <- run_coral(time=time, env=env, pars=pars)
shock_100 <- shock(healthy_run, pars, 100, 20, 480)

time_1 <- time[diff(healthy_run$S.2/healthy_run$S.1 > 1) != 0]
shock_1 <- shock(healthy_run, pars, time_1, 14, 600-20-time_1)

time_0.1 <- time[diff(healthy_run$S.2/healthy_run$S.1 > 0.1) != 0]
shock_0.1 <- shock(healthy_run, pars, time_0.1, 20, 600-20-time_0.1)

time_0.01 <- time[diff(healthy_run$S.2/healthy_run$S.1 > 0.01) != 0]
shock_0.01 <- shock(healthy_run, pars, time_0.01, 20, 600-20-time_0.01)

time_h0.5 <- time[diff(healthy_run$S.2/healthy_run$H > 0.5) != 0]
shock_h0.5 <- shock(healthy_run, pars, time_h0.5, 20, 600-20-time_h0.5)

time_h0.25 <- time[diff(healthy_run$S.2/healthy_run$H > 0.25) != 0]
shock_h0.25 <- shock(healthy_run, pars, time_h0.25, 14, 300-14-time_h0.25)

time_h0.1 <- time[diff(healthy_run$S.2/healthy_run$H > 0.1) != 0]
shock_h0.1 <- shock(healthy_run, pars, time_h0.1, 14, time_h0.1+14+100)

time[diff(healthy_run$S.2/healthy_run$S.1 > 1) != 0]
shock_54.5 <- shock(healthy_run, pars, 54.5, 14, 600-74.5)

shock_10 <- shock(healthy_run, pars, 10, 7, 600-24)

plot_symbs(healthy_run, "black", "blue", TRUE)
plot_symbs(shock_100, "gray", "lightblue")
plot_symbs(shock_82.8, "pink", "green")
plot_symbs(shock_54.5, "orange", "darkgreen")
plot_symbs(shock_10, "yellow", "cyan")
plot_symbs(shock_0.1, "red", "lightblue")
plot_symbs(shock_0.01, "red", "lightblue")
plot_symbs(shock_1, "gray", "lightblue")
plot_symbs(shock_h0.5, "yellow", "cyan")
plot_symbs(shock_h0.25, "pink", "green")
plot_symbs(shock_h0.1, "pink", "green")


png("healthy_run.png")
plot_symbs(healthy_run, "black", "blue", TRUE)
legend("topright", c("sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("black", "blue"), lwd=2, lty=c(1, 1, 2, 2))
dev.off()


png("healthy_run_light_shock_1.png")
plot_symbs(healthy_run, "black", "blue", TRUE)
abline(v=time_h0.25, col="orange", lwd=2)
abline(v=time_h0.25+14, col="orange", lwd=2)
legend("topright", c("sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("black", "blue"), lwd=2, lty=c(1, 1, 2, 2))
dev.off()

png("healthy_run_light_shock_2.png")
plot_symbs(healthy_run, "black", "blue", TRUE)
abline(v=time_h0.1, col="orange", lwd=2)
abline(v=time_h0.1+14, col="orange", lwd=2)
legend("topright", c("sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("black", "blue"), lwd=2, lty=c(1, 1, 2, 2))
dev.off()



png("shock_h0-25.png")
plot_symbs(run_preh0.25, "black", "blue", TRUE)
abline(v=time_h0.25, col="orange", lwd=2)
abline(v=time_h0.25+14, col="orange", lwd=2)
plot_symbs(run_preh0.25, "black", "blue")
plot_symbs(shock_h0.25, "black", "blue")
legend("topright", c("sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("black", "blue"), lwd=2, lty=c(1, 1, 2, 2))
dev.off()

png("shock_h0-1.png")
plot_symbs(run_preh0.1, "black", "blue", TRUE)
abline(v=time_h0.1, col="orange", lwd=2)
abline(v=time_h0.1+14, col="orange", lwd=2)
plot_symbs(run_preh0.1, "black", "blue")
plot_symbs(shock_h0.1, "black", "blue")
legend("topright", c("sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("black", "blue"), lwd=2, lty=c(1, 1, 2, 2))
dev.off()


#### Plots for Roger 1/21/2021

pars <- def_pars(nsym=2)
pars$kROS[[2]] <- 250
pars$jCPm[[2]] <- 1
pars$jSGm[[2]] <- 0.15
pars$initS <- c(0.25, 0.75)


time <- seq(0, 300, 0.1)
env <- list(L=rep(20, length(time)), N=rep(2e-6, length(time)), X=rep(2e-7, length(time)))

healthy_run <- run_coral(time=time, env=env, pars=pars)

time_h0.1 <- time[diff(healthy_run$S.2/healthy_run$H > 0.1) != 0]
shock_h0.1 <- shock(healthy_run, pars, list(L=40, N=2e-6, X=2e-7),
	list(L=20, N=2e-6, X=2e-7), time_h0.1, 14, 1000-14-time_h0.1)

time_preh0.1 <- seq(0, time_h0.1, 0.1)
env_preh0.1 <- list(L=rep(20, length(time_preh0.1)),
	N=rep(2e-6, length(time_preh0.1)), X=rep(2e-7, length(time_preh0.1)))

run_preh0.1 <- run_coral(time=time_preh0.1, env=env_preh0.1, pars=pars)


png("shock_h0-1_extended.png")
plot_symbs(run_preh0.1, "black", "blue", TRUE, xlim=c(0, 1000), ylim=c(0, 1.2),
	main="Light shock and post-bleaching\nL=20, N=2e-6, X=2e-7, shock L=40",
	sub="Light shock occurs when tolerant:host = 0.1")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
plot_symbs(run_preh0.1, "black", "blue")
plot_symbs(shock_h0.1, "black", "blue")
legend("right", c("light shock",
	"sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("#FFFF0080", "black", "blue", "black", "blue"), 
	lwd=c(4, 2, 2, 2, 2), lty=c(1, 1, 1, 2, 2))
dev.off()


png("shock_h0-1_zoom_in.png")
plot_symbs(run_preh0.1, "black", "blue", TRUE, xlim=c(time_h0.1-14, time_h0.1+28), 
	ylim=c(0.05, 0.15),
	main="Zoomed-in on light shock\nL=20, N=2e-6, X=2e-7, shock L=40",
	sub="Light shock occurs when tolerant:host = 0.1")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
plot_symbs(run_preh0.1, "black", "blue")
plot_symbs(shock_h0.1, "black", "blue")
legend("bottomleft", c("light shock",
	"sens. symb, positive gr.", "tol. symb, positive gr.",
	"sens. symb, negative gr.",
	"tol. symb, negative gr."),
	col=c("#FFFF0080", "black", "blue", "black", "blue"), 
	lwd=c(4, 2, 2, 2, 2), lty=c(1, 1, 1, 2, 2))
dev.off()

# Check that survival is impossible with food alone (no light)
time_darkness <- seq(0, 1000, 0.1)
env_darkness <- list(L=rep(0, length(time_darkness)), 
	N=rep(2e-6, length(time_darkness)), X=rep(2e-7, length(time_darkness)))
darkness_run <- run_coral(time=time_darkness, env=env_darkness, pars=pars)
plot_symbs(darkness_run, "black", "blue", TRUE)
# Need light to survive!


# In the healthy host (no light shock) are the tolerant symbionts going to 0
# or to an equilibrium value > 0?

# How can I investigate this? 
# The derivative needs to be slowing down, obviously
# But it will be slowing down anyway

time <- seq(0, 600, 0.1)
env <- list(L=rep(20, length(time)), N=rep(2e-6, length(time)), X=rep(2e-7, length(time)))

healthy_run <- run_coral(time=time, env=env, pars=pars)

plot_symbs(healthy_run, "black", "blue", TRUE, ylim=c(0, 1))

plot(time, healthy_run$S.2/healthy_run$H, type="l", ylim=c(-0.1, 0.1))
sum(sign(diff(healthy_run$S.2/healthy_run$H)) == -1)/(length(time)-1)
sum(sign(diff(healthy_run$S.1[time > 100]/healthy_run$H[time > 100])) == 1)/(length(time[time > 100])-1)
sum(sign(diff((healthy_run$S.1[time > 100] + healthy_run$S.2[time > 100])/healthy_run$H[time > 100])) == -1)/(length(time[time > 100])-1)



# Don't know for sure how I would tell this, since tol:host shouldn't reach
# 0 even if it's asymptotically approaching it. However, the change in
# tol:host is negative for every time step of the simulation, which suggests
# it's at least never bouncing around something stable
# On the other hand, after 100 days, the change in sens:host is always positive
# And the change in (sens + tol):host is always positive after 100 days


##### For Ferdi
pars <- def_pars(nsym=2)
pars$kROS[[2]] <- 250
pars$jCPm[[2]] <- 1
pars$jSGm[[2]] <- 0.15
pars$initS <- c(0.25, 0.75)


time <- seq(0, 300, 0.1)
env <- list(L=rep(20, length(time)), N=rep(2e-6, length(time)), X=rep(2e-7, length(time)))

healthy_run <- run_coral(time=time, env=env, pars=pars)

time_h0.1 <- time[diff(healthy_run$S.2/healthy_run$H > 0.1) != 0]
shock_h0.1 <- shock(healthy_run, pars, list(L=40, N=2e-6, X=2e-7),
	list(L=20, N=2e-6, X=2e-7), time_h0.1, 14, 20)

time_preh0.1 <- seq(0, time_h0.1, 0.1)
env_preh0.1 <- list(L=rep(20, length(time_preh0.1)),
	N=rep(2e-6, length(time_preh0.1)), X=rep(2e-7, length(time_preh0.1)))

run_preh0.1 <- run_coral(time=time_preh0.1, env=env_preh0.1, pars=pars)

png("fluxes_pg1.png")
par(mfrow=c(3, 3))

plot(run_preh0.1$time, run_preh0.1$jX, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="jX")
lines(shock_h0.1$time, shock_h0.1$jX)
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)

plot(run_preh0.1$time, run_preh0.1$jN, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="jN", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jN)


plot(run_preh0.1$time, run_preh0.1$rNH, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="rNH", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rNH)


plot(run_preh0.1$time, run_preh0.1$rhoN, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="rhoN", type="l", ylim=c(0, 0.04))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rhoN)


plot(run_preh0.1$time, run_preh0.1$jeC, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="jeC", type="l", ylim=c(0, 0.15))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jeC)


plot(run_preh0.1$time, run_preh0.1$jCO2, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="jCO2", type="l", ylim=c(0, 1.2))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jCO2)


plot(run_preh0.1$time, run_preh0.1$jHG, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="jHG", type="l", ylim=c(0, 0.11))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jHG)


plot(run_preh0.1$time, run_preh0.1$rCH, xlim=c(time_h0.1-14, time_h0.1+28),  xlab="time",
ylab="", main="rCH", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rCH)

plot(run_preh0.1$time, run_preh0.1$rNS.1, xlim=c(time_h0.1-14, time_h0.1+28), xlab="time",
ylab="", main="rNS", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rNS.1)
lines(run_preh0.1$time, run_preh0.1$rNS.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$rNS.2, col="blue")


dev.off()

###

par(mfrow=c(3,3))

plot(run_preh0.1$time, run_preh0.1$jL.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jL", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jL.1)
lines(run_preh0.1$time, run_preh0.1$jL.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jL.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$jCP.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jCP", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jCP.1)
lines(run_preh0.1$time, run_preh0.1$jCP.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jCP.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$jeL.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jeL", type="l", ylim=c(15, 60))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jeL.1)
lines(run_preh0.1$time, run_preh0.1$jeL.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jeL.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$jNPQ.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jNPQ", type="l")
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jNPQ.1)
lines(run_preh0.1$time, run_preh0.1$jNPQ.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jNPQ.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$jSG.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jSG", type="l", ylim=c(0, 0.2))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jSG.1)
lines(run_preh0.1$time, run_preh0.1$jSG.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jSG.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$rhoC.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="rhoC", type="l", ylim=c(0, 1.5))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rhoC.1)
lines(run_preh0.1$time, run_preh0.1$rhoC.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$rhoC.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$jST.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="jST", type="l", ylim=c(0.03, 0.05))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$jST.1)
lines(run_preh0.1$time, run_preh0.1$jST.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$jST.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$rCS.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="rCS", type="l", ylim=c(0.03, 0.06))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$rCS.1)
lines(run_preh0.1$time, run_preh0.1$rCS.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$rCS.2, col="blue")


plot(run_preh0.1$time, run_preh0.1$cROS.1, xlim=c(time_h0.1-14, time_h0.1+28), typle="l", xlab="time",
ylab="", main="cROS", type="l", ylim=c(1, 1.12))
rect(xleft=time_h0.1, ybottom=-0.2, xright=time_h0.1+14, ytop=1.4,
	col="#FFFF0080", border=NA)
lines(shock_h0.1$time, shock_h0.1$cROS.1)
lines(run_preh0.1$time, run_preh0.1$cROS.2, col="blue")
lines(shock_h0.1$time, shock_h0.1$cROS.2, col="blue")











