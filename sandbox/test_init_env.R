library("coRal")
source("init_env.R")

# Function for checking that two outputs of init_env are the same
checkingfun <- function(envcoRal, envmyfun, summary=FALSE) {
    testResults <- c(
        (all(envcoRal$X == envmyfun$X) && any(envcoRal$X == envmyfun$X)),
        (all(envcoRal$N == envmyfun$N) && any(envcoRal$N == envmyfun$N)),
        (all(envcoRal$L == envmyfun$L) && any(envcoRal$L == envmyfun$L)))
    if (summary) {
        cat("X values match?", testResults[1], "\n")
        cat("N values match?", testResults[2], "\n")
        cat("L values match?", testResults[3], "\n")
    } else {
        testResults
    }
}


# test that functional form options other than 5 (the new one), produce the same results for the new init_env and the old one
print("Test that option 1 produces the same values of (X, N, L) as the old function")
env1coRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,1))
env1new <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,1))
checkingfun(env1coRal, env1new, TRUE)

print("Test that option 2 produces the same values of (X, N, L) as the old function")
env2coRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,2))
env2new <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,2))
checkingfun(env2coRal, env2new, TRUE)

print("Test that option 3 produces the same values of (X, N, L) as the old function")
env3coRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,3))
env3new <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,3))
checkingfun(env3coRal, env3new, TRUE)

print("Test that option 4 produces the same values of (X, N, L) as the old function")
env4coRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,4))
env4new <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,4))
checkingfun(env4coRal, env4new, TRUE)


# Test that functional form option 5 (the new one), produces the same results as option 2 when L[1] > L[2] or L[4] < L[2]
print("Test that option 5 produces the same results as option 2 when L[1] > L[2]")
env2AcoRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(20,15,2))
env5Anew <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(20,15,5))
checkingfun(env2AcoRal, env5Anew, TRUE)

print("Test that option 5 produces the same results as option 2 when L[4] < L[2]")
env2BcoRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,2))
env5Bnew <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15,20,5,16))
checkingfun(env2BcoRal, env5Bnew, TRUE)

print("Test that option 5 produces the same results as option 2 L[1] > L[2] and L[4] < L[2]")
env2CcoRal <- coRal::init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(20,15,2))
env5Cnew <- init_env(time=seq(1,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(20,15,5,14))
checkingfun(env2CcoRal, env5Cnew, TRUE)


# Test that without a fourth entry for L, the first year has the same min as the others and an amplitude
# that is 1.5x the other years
print("Option 5, no year 1 max specified: first year should have peak with 1.5x amplitude of other years")
plot(seq(-1000,1000,0.1), init_env(time=seq(-1000,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15, 25, 5))$L)

# Test that with a fourth entry for L, the first year has the same min as the others and a max given by L[4]
print("Option 5, year 1 max specified: first year should peak at year 1 max (35)")
plot(seq(-1000,1000,0.1), init_env(time=seq(-1000,1000,0.1), X=c(2e-7,2e-7,0), N=c(1e-7,1e-7,0), L=c(15, 25, 5, 35))$L)