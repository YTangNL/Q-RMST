library(tidyverse)
library(readxl)
library(survival)
library(writexl)

#####
## First: check whether we can reconstruct the original (incorrect!) one minus Kaplan-Meiers
#####
cvdeath <- read_excel("cvd_ipd.xlsx")
mi <- read_excel("mi_ipd.xlsx")
revas <- read_excel("revasc_ipd.xlsx")
stroke <- read_excel("stroke_ipd.xlsx")
hosp <- read_excel("hosp_ipd.xlsx")

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
censseq <- punif(tseq, 18, 36, lower.tail = FALSE)
lines(tseq, censseq, type="l", col="blue", lty=3)

#####
## Third: generate data with multiple events under assumptions
#####

# Assume that for each of the events considered, the rate of first event is
# the estimated rate of event multiplied by an (as yet unknown) constant, called a

a1 <- c(0.616, 0.66885, 0.447, 0.82, 0.911)
a2 <- c(1.0, 0.7, 1.3, 0.3, 0.6)
a3 <- c(0.8, 0.6, 1.2, 0.2, 1.0)
a4 <- c(1.1, 0.8, 0.8, 0.18, 1.0)
a5 <- c(2.35, 3.1, 0.5, 0.25, 1.80)
  
# Also assume that for each event type, after an initial first event
# other than death, the event is double the first event rate
# maybe also need to assume another constant vector r for relapse chance 
r <- c(1, 10, 10, 20, 20) # discuss with Jan 

# Now we are going to generate multiple events for each individual, according
# to these assumptions, re-estimate cumulative events and check with the original
# data

#
# evolo
#

set.seed(2025)
n <- 13784
res1 <- matrix(NA, n, 11)
colnames(res1) <- c("id", "time1",'type1',"time2",'type2',"time3",'type3',
                   "time4",'type4','time5',"type5")
res1[, 1] <- 1:n
t1 <- rexp(n, rate = r[1]* a1[1] * tmp1$rate[2])
t2 <- rexp(n, rate = r[1]* a1[2] * tmp2$rate[2])
t3 <- rexp(n, rate = r[1]* a1[3] * tmp3$rate[2])
t4 <- rexp(n, rate = r[1]* a1[4] * tmp4$rate[2])
t5 <- rexp(n, rate = r[1]* a1[5] * tmp5$rate[2])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res1[, 2] <- tnext
res1[, 3] <- dnext
whnext <- which(dnext > 1)
nnext <- length(whnext)

# After first event, generate second event, except if dnext = 1
t1 <- rexp(nnext, rate = r[2] * a2[1] * tmp1$rate[2])
t2 <- rexp(nnext, rate = r[2] * a2[2] * tmp2$rate[2])
t3 <- rexp(nnext, rate = r[2] * a2[3] * tmp3$rate[2])
t4 <- rexp(nnext, rate = r[2] * a2[4] * tmp4$rate[2])
t5 <- rexp(nnext, rate = r[2] * a2[5] * tmp5$rate[2])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res1[whnext, 4] <- res1[whnext, 2] + tnext
res1[whnext, 5] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)


# After second event, generate third event, except if dnext = 1
t1 <- rexp(nnext, rate = r[3] * a3[1] * tmp1$rate[2])
t2 <- rexp(nnext, rate = r[3]* a3[2] * tmp2$rate[2])
t3 <- rexp(nnext, rate = r[3] * a3[3] * tmp3$rate[2])
t4 <- rexp(nnext, rate = r[3] * a3[4] * tmp4$rate[2])
t5 <- rexp(nnext, rate = r[3] * a3[5] * tmp5$rate[2])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res1[whnext, 6] <- res1[whnext, 4] + tnext
res1[whnext, 7] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After third event, generate fourth event, except if dnext = 1
t1 <- rexp(nnext, rate = r[4] * a4[1] * tmp1$rate[2])
t2 <- rexp(nnext, rate = r[4] * a4[2] * tmp2$rate[2])
t3 <- rexp(nnext, rate = r[4] * a4[3] * tmp3$rate[2])
t4 <- rexp(nnext, rate = r[4] * a4[4] * tmp4$rate[2])
t5 <- rexp(nnext, rate = r[4] * a4[5] * tmp5$rate[2])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res1[whnext, 8] <- res1[whnext, 6] + tnext
res1[whnext, 9] <- dnext
whnext <- whnext[dnext > 1]
nnext <- length(whnext)

# After fourth event, generate fifth (final) event, except if dnext = 1
t1 <- rexp(nnext, rate = r[5] * a5[1] * tmp1$rate[2])
t2 <- rexp(nnext, rate = r[5] * a5[2] * tmp2$rate[2])
t3 <- rexp(nnext, rate = r[5] * a5[3] * tmp3$rate[2])
t4 <- rexp(nnext, rate = r[5] * a5[4] * tmp4$rate[2])
t5 <- rexp(nnext, rate = r[5] * a5[5] * tmp5$rate[2])
tnext <- pmin(t1, t2, t3, t4, t5)
dnext <- 1
dnext[tnext == t1] <- 1
dnext[tnext == t2] <- 2
dnext[tnext == t3] <- 3
dnext[tnext == t4] <- 4
dnext[tnext == t5] <- 5
res1[whnext, 10] <- res1[whnext, 8] + tnext
res1[whnext, 11] <- dnext


#####
## Fourth: Censor data 
#####
set.seed(2025)
ctime <- runif(n, 18, 36)
res1 <- cbind(res1,ctime)
res1 <- data.frame(res1)

# censor first event
res1$type1 <- ifelse(res1$time1 < ctime, res1$type1, 0)
res1$time1 <- ifelse(res1$time1 > res1$ctime, res1$ctime, res1$time1)
res1[which(res1$type1 == 0), 4:(ncol(res1)-1)] <- NA
# censor second event
res1$type2 <- ifelse(res1$time2 < ctime, res1$type2, 0)
res1$time2 <- ifelse(res1$time2 > res1$ctime, res1$ctime, res1$time2)
res1[which(res1$type2 == 0), 6:(ncol(res1)-1)] <- NA
# censor third event
res1$type3 <- ifelse(res1$time3 < ctime, res1$type3, 0)
res1$time3 <- ifelse(res1$time3 > res1$ctime, res1$ctime, res1$time3)
res1[which(res1$type3 == 0), 8:(ncol(res1)-1)] <- NA
# censor fourth event
res1$type4 <- ifelse(res1$time4 < ctime, res1$type4, 0)
res1$time4 <- ifelse(res1$time4 > res1$ctime, res1$ctime, res1$time4)
res1[which(res1$type4 == 0), 10:(ncol(res1)-1)] <- NA
# censor fifth event
res1$type5 <- ifelse(res1$time5 < ctime, res1$type5, 0)
res1$time5 <- ifelse(res1$time5 > res1$ctime, res1$ctime, res1$time5)

# number of events each time
n_1st <- sum(res1$type1 > 0, na.rm = TRUE)
n_1st

n_2nd <- sum(res1$type2 > 0, na.rm = TRUE)
n_2nd

n_3rd <- sum(res1$type3 > 0, na.rm = TRUE)
n_3rd

n_4th <- sum(res1$type4 > 0, na.rm = TRUE)
n_4th

n_5th <- sum(res1$type5 > 0, na.rm = TRUE)
n_5th

n_total <- n_1st + n_2nd + n_3rd + n_4th + n_5th
n_total

# number of causes specific event frequencies
count_1st <- table(res1$type1)
count_1st

count_2nd <- table(res1$type2)
count_2nd

count_3rd <- table(res1$type3)
count_3rd

count_4th <- table(res1$type4)
count_4th

count_5th <- table(res1$type5)
count_5th

write_xlsx(res1,"censor_evo.xlsx")
#####
## Fifth: Extract Key counting numbers
#####

# Numbers from the clinical report
ncum_repo_evo <- c(251,468,759,207,236)
n1_repo_evo <- c(161,329,349,184,169)
n_event_repo_evo <- c(768,417,159)

# Numbers from simulation
n_multi <- c(1344-sum(count_1st[2:6]),576-sum(count_2nd[2:6]),159-sum(count_3rd[2:6]),178-sum(count_4th[2:6])-sum(count_5th[2:6]))
n_total_evo <- n_1st + n_2nd + n_3rd + n_4th + n_5th
n_1event_evo <- (n_1st + n_multi[1]) - (n_2nd + n_multi[2])
n_2event_evo <- (n_2nd + n_multi[2]) - (n_3rd + n_multi[3])
n_3pevent_evo <- n_3rd + n_multi[3]

n_events_sim = c(n_1event_evo,n_2event_evo,n_3pevent_evo)
count_sim <- rbind(count_1st,count_2nd, count_3rd,count_4th,count_5th)
count_cum_sim <- colSums(count_sim)

# output a comparison table
compare_evo <- data.frame(cumu_events = c(ncum_repo_evo,1921),
                          cumu_sim = c(count_cum_sim[2:6],sum(count_cum_sim[2:6])),
                          first_event = c(n1_repo_evo,1192),
                          first_sim = c(as.matrix(count_1st[2:6]),sum(count_1st[2:6])),
                          row.names = c("cvd", "mi", "revasc.","stroke", "hospital.","Total"))

# add multiple events same time term
compare_sub <- data.frame (report = n_event_repo_evo,
                sim = n_events_sim,
                row.names = c("subjets with 1 event", "subjets with 2 events", "subjets with 3+ event"))

  