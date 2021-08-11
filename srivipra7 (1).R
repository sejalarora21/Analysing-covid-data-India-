
library(deSolve)


################### Phase 1 ##########################################
Infctd <- read.csv("https://raw.githubusercontent.com/sejalarora21/SRI-VIPRA/main/covid9.csv",header = F,sep = ",")
Infected <- Infctd[,6]
Day <- 1:(length(Infected))
N <- 1000000000 # population of the India


old <- par(mfrow = c(1, 2))

plot(Day, Infected, type ="b")
plot(Day, Infected, log = "y")
abline(lm(log10(Infected) ~ Day))
title("Confirmed infections COVID-19 in India (Oct-Nov)", outer = TRUE, line = -2)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init <- c(S = N-6312584, I = Infected[1], R = 5273201)
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}


Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
print(Opt_par)

r0<- as.numeric(Opt_par[1]/Opt_par[2])
r0

t <- 1:120 #time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
col <- 1:3 # colour
print(fit) 
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col)
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col, log = "y")

legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.05)
title("Predicted Cases 2019-nCoV UK (worst case)", outer = TRUE, line = -2)


####################### Phase 2: May-July #########################################################
####################### Without Vaccination: Worst Case scenario #####################################


Infctd2 <- read.csv("https://raw.githubusercontent.com/sejalarora21/SRI-VIPRA/main/Covid13%20-%20may-july.csv",header = F,sep = ",")
Infected2 <- Infctd2[,6]
Day2 <- 1:(length(Infected2))
N2 <- 1000000000 # population of the India

old <- par(mfrow = c(1, 2))


SIR2 <- function(time, state, parameters2) {
  par <- as.list(c(state, parameters2))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init2 <- c(S = N2-19164969, I = Infected2[1], R = 15684406)
RSS2 <- function(parameters2) {
  names(parameters2) <- c("beta", "gamma")
  out <- ode(y = init2, times = Day2, func = SIR2, parms = parameters2)
  fit2 <- out[ , 3]
  sum((Infected2 - fit2)^2)
}


Opt2 <- optim(c(0.5, 0.5), RSS2, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt2$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

Opt_par2 <- setNames(Opt2$par, c("beta", "gamma"))
print(Opt_par2)

r0_new2<- as.numeric(Opt_par2[1]/Opt_par2[2])
r0_new2

t2 <- 1:120 #time in days
fit2 <- data.frame(ode(y = init2, times = t2, func = SIR2, parms = Opt_par2))
col <- 1:3 # colour
print(fit2) 
matplot(fit2$time, fit2[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col)
matplot(fit2$time, fit2[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col, log = "y")


####################### Phase 2: may june #########################################################
####################### With Vaccination: 50% Effective #####################################

Infctd3 <- read.csv("https://raw.githubusercontent.com/sejalarora21/SRI-VIPRA/main/Covid13%20-%20may-july.csv",header = F,sep = ",")
Infected3 <- Infctd3[,6]
Day3 <- 1:(length(Infected3))
N3 <- 1000000000 # population of the India

old <- par(mfrow = c(1, 2))


SIR3 <- function(time, state, parameters3) {
  par <- as.list(c(state, parameters3))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init3 <- c(S = N3-32813647, I = Infected3[1], R = 29333084)
RSS3 <- function(parameters3) {
  names(parameters3) <- c("beta", "gamma")
  out <- ode(y = init3, times = Day3, func = SIR3, parms = parameters3)
  fit3 <- out[ , 3]
  sum((Infected3 - fit3)^2)
}

parameters3<- c(beta = 0.25106961, gamma = 0.1357772)

r0_new3<- 0.25106961/0.1357772
r0_new3

t3 <- 1:120 #time in days
fit3 <- data.frame(ode(y = init3, times = t3, func = SIR3, parms = parameters3))
col <- 1:3 # colour
print(fit3) 
matplot(fit3$time, fit3[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col)
legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.005)
matplot(fit3$time, fit3[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col, log = "y")
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.005)
title("Predicted Cases of Covid-19, India in August'21 (50% effectiveness of vaccine)", outer = TRUE, line = -2)
####################### Phase 2: may june #########################################################
####################### With Vaccination: 70% Effective #####################################

Infctd4 <- read.csv("https://raw.githubusercontent.com/sejalarora21/SRI-VIPRA/main/Covid13%20-%20may-july.csv",header = F,sep = ",")
Infected4 <- Infctd4[,6]
Day4 <- 1:(length(Infected4))
N4<- 1000000000 # population of the India

old <- par(mfrow = c(1, 2))


SIR4 <- function(time, state, parameters4) {
  par <- as.list(c(state, parameters4))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init4 <- c(S = N4-38273119, I = Infected4[1], R = 34792556)
RSS4<- function(parameters4) {
  names(parameters4) <- c("beta", "gamma")
  out <- ode(y = init4, times = Day4, func = SIR4, parms = parameters4)
  fit4 <- out[ , 3]
  sum((Infected4 - fit4)^2)
}


parameters4<- c(beta = 0.23785542, gamma = 0.1357772)

r0_new4<- 0.23785542/0.1357772
r0_new4

t4 <- 1:120 #time in days
fit4 <- data.frame(ode(y = init4, times = t4, func = SIR4, parms = parameters4))
col <- 1:3 # colour
print(fit4) 
matplot(fit4$time, fit4[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col)
legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.005)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col,
       xpd=TRUE, inset=c(0, -.38), cex=.8)
matplot(fit4$time, fit4[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col, log = "y")
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.005)
title("Predicted Cases of Covid-19, India in August'21 (70% effectiveness of vaccine)", outer = TRUE, line = -2)

####################### Phase 2: may june july #########################################################
####################### With Vaccination: 95% Effective #####################################

Infctd5 <- read.csv("https://raw.githubusercontent.com/sejalarora21/SRI-VIPRA/main/Covid13%20-%20may-july.csv",header = F,sep = ",")
Infected5 <- Infctd5[,6]
Day5 <- 1:(length(Infected5))
N5<- 1000000000 # population of the India

old <- par(mfrow = c(1, 2))


SIR5 <- function(time, state, parameters5) {
  par <- as.list(c(state, parameters5))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init5 <- c(S = N5-45097458, I = Infected5[1], R = 41616895)
RSS5<- function(parameters5) {
  names(parameters5) <- c("beta", "gamma")
  out <- ode(y = init5, times = Day5, func = SIR5, parms = parameters5)
  fit5 <- out[ , 3]
  sum((Infected5 - fit5)^2)
}


parameters5<- c(beta = 0.21142704, gamma = 0.1357772)

r0_new5<-0.21142704/0.1357772
r0_new5

t5<- 1:120 #time in days
fit5 <- data.frame(ode(y = init5, times = t5, func = SIR5, parms = parameters5))
col <- 1:3 # colour
print(fit5) 
matplot(fit5$time, fit5[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col,
       xpd=TRUE, inset=c(0, -.38), cex=.8)
matplot(fit5$time, fit5[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 1.5, lty = 1, col = col, log = "y")
title("Predicted Cases of Covid-19, India in August'21 (95% effectiveness of vaccine)", outer = TRUE, line = -2)
