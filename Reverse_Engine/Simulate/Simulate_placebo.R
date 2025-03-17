library(tidyverse)
library(readxl)
library(survival)
library(writexl)

#####
## First: check whether we can reconstruct the original (incorrect!) one minus Kaplan-Meiers
#####
cvdeath <- read_excel("Reverse_Engine/Simulate/cvd_ipd.xlsx")
mi <- read_excel("Reverse_Engine/Simulate/mi_ipd.xlsx")
revas <- read_excel("Reverse_Engine/Simulate/revasc_ipd.xlsx")
stroke <- read_excel("Reverse_Engine/Simulate/stroke_ipd.xlsx")
hosp <- read_excel("Reverse_Engine/Simulate/hosp_ipd.xlsx")

## Also: check the numbers of events
cvd_event <- table(cvdeath$status, cvdeath$treat)[2,]
mi_event <- table(mi$status, mi$treat)[2,]
revas_event <- table(revas$status, revas$treat)[2,]
stroke_event <- table(stroke$status, stroke$treat)[2,]
hosp_event <- table(hosp$status, hosp$treat)[2,]
t(rbind(mi_event, revas_event,cvd_event,stroke_event,hosp_event))

#####
## Second: estimate constant hazard rates for each of the events
#####

# CV death 
tmp <- aggregate(cbind(time, status) ~ treat, data = cvdeath, sum)
tmp$rate <- tmp$status / tmp$time
tmp
tmp1 <- tmp
tseq <- seq(0, 42, by = 0.05)

# Plot cumulative incidence function
sf_cvd <- survfit(Surv(time, status) ~ treat, data = cvdeath)
plot(sf_cvd, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of CvDeath", xaxt = "n")
axis(1, at = c(0,6,12,18,24,30,36,42))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# MI 
tmp <- aggregate(cbind(time, status) ~ treat, data = mi, sum)
tmp$rate <- tmp$status / tmp$time
tmp
tmp2 <- tmp
tseq <- seq(0, 42, by = 0.05)

# Plot cumulative incidence function
sf_mi <- survfit(Surv(time, status) ~ treat, data = mi)
plot(sf_mi, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of MI", xaxt = "n")
axis(1, at = c(0,6,12,18,24,30,36,42))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# Revascularization 
tmp <- aggregate(cbind(time, status) ~ treat, data = revas, sum)
tmp$rate <- tmp$status / tmp$time
tmp
tmp3 <- tmp
tseq <- seq(0, 42, by = 0.05)

# Plot cumulative incidence function
sf_revas <- survfit(Surv(time, status) ~ treat, data = revas)
plot(sf_revas, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of Revasc.", xaxt = "n")
axis(1, at = c(0,6,12,18,24,30,36,42))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# Stroke 
tmp <- aggregate(cbind(time, status) ~ treat, data = stroke, sum)
tmp$rate <- tmp$status / tmp$time
tmp
tmp4 <- tmp
tseq <- seq(0, 42, by = 0.05)

# Plot cumulative incidence function
sf_stroke <- survfit(Surv(time, status) ~ treat, data = stroke)
plot(sf_stroke, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of Stroke", xaxt = "n")
axis(1, at = c(0,6,12,18,24,30,36,42))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

# 5. Hospitalization
tmp <- aggregate(cbind(time, status) ~ treat, data = hosp, sum)
tmp$rate <- tmp$status / tmp$time
tmp
tmp5 <- tmp
tseq <- seq(0, 42, by = 0.05)

# Plot cumulative incidence function
sf_hosp <- survfit(Surv(time, status) ~ treat, data = hosp)
plot(sf_hosp, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of Hospitalization", xaxt = "n")
axis(1, at = c(0,6,12,18,24,30,36,42))
lines(tseq, 1 - exp(-tmp$rate[1] * tseq), type="l", col = "red", lty=3)
lines(tseq, 1 - exp(-tmp$rate[2] * tseq), type="l", col = "blue", lty=3)

####
## We also need a censoring distribution
####
sfcens <- survfit(Surv(time, status==0) ~ treat, data=cvdeath)
plot(sfcens)
# This suggests uniform censoring between 20 and 36 months
censseq <- punif(tseq, 17, 36, lower.tail = FALSE)
lines(tseq, censseq, type="l", col="blue", lty=3)

#####
## Third: generate data with multiple events under assumptions
#####

# Assume that for each of the events considered, the rate of first event is
# the estimated rate of event multiplied by an (as yet unknown) constant, called a
# by cause
a_cvd <- c(0.598, 12,  4, 5, 15)
a_mi <- c(0.636,  8,  4, 6, 8.5)
a_revas <- c(0.412, 13,  10, 20, 31)
a_stroke <- c(0.812,  3,  2,  4, 1.5)
a_hosp <- c(0.835,  9, 7, 12, 19.5)
# by time
#a1 <- c(0.585, 0.624, 0.4062, 0.8, 0.82)
#a2 <- c(12, 8, 13, 4, 9)
#a3 <- c(8, 4,  9,  1, 10)
#a4 <- c(14, 12, 20,  8, 30)
#a5 <- c(12.4, 15.0, 13.6, 10.0, 32.0)
# Also assume that for each event type, after an initial first event
# other than death, the event is double the first event rate
# maybe also need to assume another constant vector r for relapse chance 

# r <- c(1, 10, 10, 20, 20) 
# Now we are going to generate multiple events for each individual, according
# to these assumptions, re-estimate cumulative events and check with the original
# data

#
# Placebo
#
set.seed(2025)
n <- 13780
res <- matrix(NA, n, 11)
colnames(res) <- c("id", "time1",'type1',"time2",'type2',"time3",'type3',
                    "time4",'type4','time5',"type5")
res[, 1] <- 1:n
t1 <- rexp(n, rate = a_cvd[1] * tmp1$rate[1])
t2 <- rexp(n, rate = a_mi[1] * tmp2$rate[1])
t3 <- rexp(n, rate = a_revas[1] * tmp3$rate[1])
t4 <- rexp(n, rate = a_stroke[1] * tmp4$rate[1])
t5 <- rexp(n, rate = a_hosp[1] * tmp5$rate[1])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res[, 2] <- tnext
res[, 3] <- dnext
whnext <- which(dnext > 1)
nnext <- length(whnext)

# After first event, generate second event, except if dnext = 1
t1 <- rexp(nnext, rate = a_cvd[2] * tmp1$rate[1])
t2 <- rexp(nnext, rate = a_mi[2] * tmp2$rate[1])
t3 <- rexp(nnext, rate = a_revas[2] * tmp3$rate[1])
t4 <- rexp(nnext, rate = a_stroke[2] * tmp4$rate[1])
t5 <- rexp(nnext, rate = a_hosp[2] * tmp5$rate[1])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res[whnext, 4] <- res[whnext, 2] + tnext
res[whnext, 5] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After second event, generate third event, except if dnext = 1
t1 <- rexp(nnext, rate = a_cvd[3] * tmp1$rate[1])
t2 <- rexp(nnext, rate = a_mi[3] * tmp2$rate[1])
t3 <- rexp(nnext, rate = a_revas[3] * tmp3$rate[1])
t4 <- rexp(nnext, rate = a_stroke[3] * tmp4$rate[1])
t5 <- rexp(nnext, rate = a_hosp[3] * tmp5$rate[1])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res[whnext, 6] <- res[whnext, 4] + tnext
res[whnext, 7] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After third event, generate fourth event, except if dnext = 1
t1 <- rexp(nnext, rate = a_cvd[4] * tmp1$rate[1])
t2 <- rexp(nnext, rate = a_mi[4] * tmp2$rate[1])
t3 <- rexp(nnext, rate = a_revas[4] * tmp3$rate[1])
t4 <- rexp(nnext, rate = a_stroke[4] * tmp4$rate[1])
t5 <- rexp(nnext, rate = a_hosp[4] * tmp5$rate[1])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res[whnext, 8] <- res[whnext, 6] + tnext
res[whnext, 9] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After fourth event, generate fifth (final) event, except if dnext = 1
t1 <- rexp(nnext, rate = a_cvd[5] * tmp1$rate[1])
t2 <- rexp(nnext, rate = a_mi[5] * tmp2$rate[1])
t3 <- rexp(nnext, rate = a_revas[5] * tmp3$rate[1])
t4 <- rexp(nnext, rate = a_stroke[5] * tmp4$rate[1])
t5 <- rexp(nnext, rate = a_hosp[5] * tmp5$rate[1])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res[whnext, 10] <- res[whnext, 8] + tnext
res[whnext, 11] <- dnext

#####
## Fourth: Censor data
#####

#
# placebo
#

set.seed(2025)
ctime <- runif(n, 17, 36)
res <- cbind(res,ctime)
res <- data.frame(res)

# censor first event
res$type1 <- ifelse(res$time1 < ctime, res$type1, 0)
res$time1 <- ifelse(res$time1 > res$ctime, res$ctime, res$time1)
res[which(res$type1 == 0), 4:(ncol(res)-1)] <- NA
# censor second event
res$type2 <- ifelse(res$time2 < ctime, res$type2, 0)
res$time2 <- ifelse(res$time2 > res$ctime, res$ctime, res$time2)
res[which(res$type2 == 0), 6:(ncol(res)-1)] <- NA
# censor third event
res$type3 <- ifelse(res$time3 < ctime, res$type3, 0)
res$time3 <- ifelse(res$time3 > res$ctime, res$ctime, res$time3)
res[which(res$type3 == 0), 8:(ncol(res)-1)] <- NA
# censor fourth event
res$type4 <- ifelse(res$time4 < ctime, res$type4, 0)
res$time4 <- ifelse(res$time4 > res$ctime, res$ctime, res$time4)
res[which(res$type4 == 0), 10:(ncol(res)-1)] <- NA
# censor fifth event
res$type5 <- ifelse(res$time5 < ctime, res$type5, 0)
res$time5 <- ifelse(res$time5 > res$ctime, res$ctime, res$time5)

# number of events each time
n_1st <- sum(res$type1 > 0, na.rm = TRUE)
n_1st

n_2nd <- sum(res$type2 > 0, na.rm = TRUE)
n_2nd

n_3rd <- sum(res$type3 > 0, na.rm = TRUE)
n_3rd

n_4th <- sum(res$type4 > 0, na.rm = TRUE)
n_4th

n_5th <- sum(res$type5 > 0, na.rm = TRUE)
n_5th

n_total <- n_1st + n_2nd + n_3rd + n_4th + n_5th
n_total

# number of causes specific event frequencies
count_1st <- table(res$type1)
count_1st

count_2nd <- table(res$type2)
count_2nd

count_3rd <- table(res$type3)
count_3rd

count_4th <- table(res$type4)
count_4th

count_5th <- table(res$type5)
count_5th

write_xlsx(res,"censor_pla.xlsx")
#####
## Fifth: Extract Key counting numbers
#####

# Numbers from the clinical report
ncum_repo_pla <- c(240,639,965,262,239)
n1_repo_pla <- c(142,423,394,226,160)
n_event_repo_pla <- c(806,543,214)

# Numbers from simulation
n_multi <- c(1563-sum(count_1st[2:6]),757-sum(count_2nd[2:6]),214-sum(count_3rd[2:6]),178-sum(count_4th[2:6])-sum(count_5th[2:6]))
n_total_pla <- n_1st + n_2nd + n_3rd + n_4th + n_5th
n_1event_pla <- (n_1st + n_multi[1]) - (n_2nd + n_multi[2])
n_2event_pla <- (n_2nd + n_multi[2]) - (n_3rd + n_multi[3])
n_3pevent_pla <- n_3rd + n_multi[3]

n_events_sim = c(n_1event_pla,n_2event_pla,n_3pevent_pla)
count_sim <- rbind(count_1st,count_2nd, count_3rd,count_4th,count_5th)
count_cum_sim <- colSums(count_sim)

# output a comparison table
compare_pla <- data.frame(cumu_events = c(ncum_repo_pla,2345),
                          cumu_sim = c(count_cum_sim[2:6],sum(count_cum_sim[2:6])),
                          first_event = c(n1_repo_pla,1345),
                          first_sim = c(as.matrix(count_1st[2:6]),sum(count_1st[2:6])),
                          row.names = c("cvd", "mi", "revasc.","stroke", "hospital.","Total"))

compare_sub <- data.frame (report = n_event_repo_pla,
                           sim = n_events_sim,
                           row.names = c("subjets with 1 event", "subjets with 2 events", "subjets with 3+ event"))

