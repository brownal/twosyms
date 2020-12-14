# This script makes some plots for presentations
#
# Ideally something like "symbiont_dominance_standard_sens_tol_initS_1.Rmd"
# would be run before running this script, because it requires a set of simulations
# and a DAPC analysis to already be run.
#
# Required:
# * healthy_sims: simulations of hosts with sensitive, tolerant, and both symbionts,
#	initialized with a "healthy" (non-bleached) symbiont:host ratio (e.g. 1).
#	healthy_sims should be a dataframe with the following columns:
#	- L and N: the light and nitrogen levels the simulation was run at
#	- X: the food level the simulation was run at. To match the current plots
#		the 4th unique value of X should be 1.2e-7
#	- kCO2: the host carbon-concentrating ability. To match the current plots
#		kCO2 should at least take a value of 10 (can take other values too)
# 	- survival_sen: whether a host with a sensitive symbiont survives for the
#		given values of L, N, X, and kCO2
#	- survival_tol, survival_both: same as survival_sen but for hosts with
#		only tolerant symbionts and both sensitive and tolerant symbionts
#	- tol_per_sen: how many tolerant symbionts per sensitive symbiont at the end
#		of each simulation run
#
# * mydapc_thinned: a dapc object produced by running a DAPC on which symbiont dominates
#	a host when both are present, based on environmental conditions (light, nitrogen,
#	and food) and host traits (kCO2 = carbon-concentrating ability). The possible
#	host states are sensitive-dominated, tolerant-dominated, and no survival (host
#	cannot survive in these environmental conditions)

#########################################################################################

# Example plot of host survival with sensitive symbionts

my_sims <- healthy_sims[healthy_sims$kCO2 == 10 &
	healthy_sims$X == unique(healthy_sims$X)[[4]],]
my_mat <- matrix(FALSE, nrow=length(unique(healthy_sims$N)),
	ncol=length(unique(healthy_sims$L)), dimnames = list(unique(healthy_sims$N),
	unique(healthy_sims$L)))

for (i in 1:length(unique(my_sims$N))) {
      for (j in 1:length(unique(my_sims$L))) {
             my_mat[[i, j]] <- my_sims[(my_sims$L == unique(my_sims$L)[[j]]) &
			(my_sims$N == unique(my_sims$N)[[i]]), "survival_sen"]
      }
}

png("Survive_sens_alone.png")
plot(my_mat[nrow(my_mat):1,2:ncol(my_mat)], xlab="Light", ylab="Nitrogen",
	col=c("gray", "#648FFF"), main="kCO2 = 10, X=1.2e-7")
dev.off()


# Example plot of host survival with tolerant symbionts

my_mat2 <- matrix(FALSE, nrow=length(unique(healthy_sims$N)),
	ncol=length(unique(healthy_sims$L)), dimnames = list(unique(healthy_sims$N),
	unique(healthy_sims$L)))

for (i in 1:length(unique(my_sims$N))) {
      for (j in 1:length(unique(my_sims$L))) {
             my_mat2[[i, j]] <- my_sims[(my_sims$L == unique(my_sims$L)[[j]]) &
			(my_sims$N == unique(my_sims$N)[[i]]), "survival_tol"]
      }
}

png("Survive_tol_alone.png")
plot(my_mat2[nrow(my_mat2):1,2:ncol(my_mat2)], xlab="Light", ylab="Nitrogen",
	col=c("gray", "#FFB000"), main="kCO2 = 10, X=1.2e-7")
dev.off()


# Example plot of host dominance by sensitive/tolerant symbionts when both are present

my_mat3 <- matrix("Neither", nrow=length(unique(healthy_sims$N)),
	ncol=length(unique(healthy_sims$L)), dimnames = list(unique(healthy_sims$N),
	unique(healthy_sims$L)))

for (i in 1:length(unique(my_sims$N))) {
      for (j in 1:length(unique(my_sims$L))) {
		my_mat3[[i, j]] <- "No survival"
		liveQ <-
             my_sims[(my_sims$L == unique(my_sims$L)[[j]]) &
			(my_sims$N == unique(my_sims$N)[[i]]), "survival_both"]
		tolPsen <- my_sims[(my_sims$L == unique(my_sims$L)[[j]]) &
			(my_sims$N == unique(my_sims$N)[[i]]), "tol_per_sen"]
		if(liveQ) {
			my_mat3[[i, j]] <- ifelse(tolPsen > 1, "Tolerant",
				"Sensitive")
		}
      }
}

png("dominance.png")
plot(my_mat3[nrow(my_mat3):1,2:ncol(my_mat3)], xlab="Light", ylab="Nitrogen",
	col=c("gray", "#648FFF", "#FFB000"), main="kCO2 = 10, X=1.2e-7")
dev.off()


# Example plot of host survival with both symbionts present

my_mat4 <- matrix(FALSE, nrow=length(unique(healthy_sims$N)),
	ncol=length(unique(healthy_sims$L)), dimnames = list(unique(healthy_sims$N),
	unique(healthy_sims$L)))

for (i in 1:length(unique(my_sims$N))) {
      for (j in 1:length(unique(my_sims$L))) {
             my_mat4[[i, j]] <- my_sims[(my_sims$L == unique(my_sims$L)[[j]]) &
			(my_sims$N == unique(my_sims$N)[[i]]), "survival_both"]
      }
}

png("Survive_both.png")
plot(my_mat4[nrow(my_mat4):1,2:ncol(my_mat4)], xlab="Light", ylab="Nitrogen",
	col=c("gray", "pink"), main="kCO2 = 10, X=1.2e-7")
dev.off()


# Example simulation time series

simple_ts <- coRal::run_coral(time=seq(0, 365, 0.1), env=coRal::init_env(time=seq(0, 365, 0.1), L=c(20, 20, 0), N=c(2e-6, 2e-6, 0), X=c(2e-7, 2e-7, 0)),
      pars=pars_both)

png("simple_ts.png")
plot(simple_ts$time, simple_ts$S.1/simple_ts$H, type="l", col="#648FFF",
      xlab="Time (days)", ylab="Symbiont:host ratio", ylim=c(0, 0.2), xlim=c(1, 365))
lines(simple_ts$time, simple_ts$S.2/simple_ts$H, col="#FFB000", lwd=2)
legend("right", c("Sensitive symb.", "Tolerant symb."), col=c("#648FFF", "#FFB000"), lwd=c(1, 2))
dev.off()


# Plot which symbiont dominates each host projected onto the first 2 axes found by DAPC

png("dapc.png")
scatter(mydapc_thinned, cell = 0, cstar = 0, mstree = FALSE, lwd = 2, lty = 2, axesell=TRUE, scree.da=FALSE,
      col=c("#000000", "#648FFF", "#FFB000"), pch=c(0, 2, 16))
arrows(x0=0, y0=0, x1=mydapc_thinned$var.load[1,1], y1=mydapc_thinned$var.load[1,2])
text(x=mydapc_thinned$var.load[1,1], y=mydapc_thinned$var.load[1,2], labels="Light", pos=1)
arrows(x0=0, y0=0, x1=mydapc_thinned$var.load[2,1], y1=mydapc_thinned$var.load[2,2])
text(x=mydapc_thinned$var.load[2,1], y=mydapc_thinned$var.load[2,2], labels="Nitrogen", pos=1)
arrows(x0=0, y0=0, x1=mydapc_thinned$var.load[3,1], y1=mydapc_thinned$var.load[3,2])
text(x=mydapc_thinned$var.load[3,1], y=mydapc_thinned$var.load[3,2], labels="Food", pos=2)
arrows(x0=0, y0=0, x1=mydapc_thinned$var.load[4,1], y1=mydapc_thinned$var.load[4,2])
text(x=mydapc_thinned$var.load[4,1], y=mydapc_thinned$var.load[4,2], labels="KCO2", pos=3)
summary(mydapc_thinned)
dev.off()
