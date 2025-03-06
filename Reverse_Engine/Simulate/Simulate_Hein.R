library(tidyverse)
library(readxl)
library(survival)

#####
## First: check whether we can reconstruct the original (incorrect!) one minus Kaplan-Meiers
#####

cvdeath <- read_excel("CVDeath IPD.xlsx")
sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = cvdeath)
plot(sf, fun="event", col = c("red", "blue"))

mi <- read_excel("FirstMI IPD.xls")
sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = mi)
plot(sf, fun="event", col = c("red", "blue"))

revas <- read_excel("FirstRevas IPD.xls")
sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = revas)
plot(sf, fun="event", col = c("red", "blue"))

## Also: check the numbers of events
table(cvdeath$event_ipd, cvdeath$arm)
table(mi$event_ipd, mi$arm)
table(revas$event_ipd, revas$arm)

#####
## Second: estimate constant hazard rates for each of the events
#####

# CV death 
tmp <- aggregate(cbind(t_ipd, event_ipd) ~ arm, data = cvdeath, sum)
tmp$rate <- tmp$event_ipd / tmp$t_ipd
tmp
tmp1 <- tmp
tseq <- seq(0, 40, by = 0.05)

sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = cvdeath)
plot(sf, fun="event", col = c("red", "blue"))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# MI
tmp <- aggregate(cbind(t_ipd, event_ipd) ~ arm, data = mi, sum)
tmp$rate <- tmp$event_ipd / tmp$t_ipd
tmp
tmp2 <- tmp

sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = mi)
plot(sf, fun="event", col = c("red", "blue"))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# Revascularization
tmp <- aggregate(cbind(t_ipd, event_ipd) ~ arm, data = revas, sum)
tmp$rate <- tmp$event_ipd / tmp$t_ipd
tmp
tmp3 <- tmp

sf <- survfit(Surv(t_ipd, event_ipd) ~ arm, data = revas)
plot(sf, fun="event", col = c("red", "blue"))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# We also need a censoring distribution???



#####
## Third: generate data with multiple events under assumptions
#####

# Assume that for each of the events considered, the rate of first event is
# the estimated rate of event multiplied by an (as yet unknown) constant, called a

# Also assume that for each event type, after an initial first event
# other than death, the event is double the first event rate

a <- 0.8

# Now we are going to generate multiple events for each individual, according
# to these assumptions, reestimate cumulative events and check with the original
# data

#
# Placebo
#

n <- 13780
res <- matrix(NA, n, 11)
res[, 1] <- 1:n
t1 <- rexp(n, rate = a * tmp1$rate[1])
t2 <- rexp(n, rate = a * tmp2$rate[1])
t3 <- rexp(n, rate = a * tmp3$rate[1])
tnext <- pmin(t1, t2, t3)
dnext <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t1] <- 1
res[, 2] <- tnext
res[, 3] <- dnext
whnext <- which(dnext > 1)
nnext <- length(whnext)

# After first event, generate second event, except if dnext = 1
t1 <- rexp(nnext, rate = 2 * a * tmp1$rate[1])
t2 <- rexp(nnext, rate = 2 * a * tmp2$rate[1])
t3 <- rexp(nnext, rate = 2 * a * tmp3$rate[1])
tnext <- pmin(t1, t2, t3)
dnext <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t1] <- 1
res[whnext, 4] <- res[whnext, 2] + tnext
res[whnext, 5] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After second event, generate third event, except if dnext = 1
t1 <- rexp(nnext, rate = 2 * a * tmp1$rate[1])
t2 <- rexp(nnext, rate = 2 * a * tmp2$rate[1])
t3 <- rexp(nnext, rate = 2 * a * tmp3$rate[1])
tnext <- pmin(t1, t2, t3)
dnext <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t1] <- 1
res[whnext, 6] <- res[whnext, 4] + tnext
res[whnext, 7] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After third event, generate fourth event, except if dnext = 1
t1 <- rexp(nnext, rate = 2 * a * tmp1$rate[1])
t2 <- rexp(nnext, rate = 2 * a * tmp2$rate[1])
t3 <- rexp(nnext, rate = 2 * a * tmp3$rate[1])
tnext <- pmin(t1, t2, t3)
dnext <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t1] <- 1
res[whnext, 8] <- res[whnext, 6] + tnext
res[whnext, 9] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After fourth event, generate fifth (final) event, except if dnext = 1
t1 <- rexp(nnext, rate = 2 * a * tmp1$rate[1])
t2 <- rexp(nnext, rate = 2 * a * tmp2$rate[1])
t3 <- rexp(nnext, rate = 2 * a * tmp3$rate[1])
tnext <- pmin(t1, t2, t3)
dnext <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t1] <- 1
res[whnext, 10] <- res[whnext, 8] + tnext
res[whnext, 11] <- dnext

#
# Evolo
#
n1 <- 13784
t1 <- rexp(n, rate = a * tmp1$rate[2])
t2 <- rexp(n, rate = a * tmp2$rate[2])
t3 <- rexp(n, rate = a * tmp3$rate[2])
