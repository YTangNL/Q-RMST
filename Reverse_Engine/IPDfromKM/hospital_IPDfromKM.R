library(IPDfromKM)
library(haven)
library(survival)
library(survminer)
library(writexl)
library(readxl)
####
#Hospitalization
####
data0 <- read_excel("/Users/miyake/Library/CloudStorage/GoogleDrive-yixiaotang0626@gmail.com/My Drive/LUMC Master Thesis/Data/coordinate/hospital/placebo_hospital.xlsx")
data1 <- read_excel("/Users/miyake/Library/CloudStorage/GoogleDrive-yixiaotang0626@gmail.com/My Drive/LUMC Master Thesis/Data/coordinate/hospital/evolo_hospital.xlsx")

####
##placebo
####
coordinate <- data0[,c(1,5)]
trisk <- data0$trisk[1:8]
nrisk <- data0$nrisk[1:8]
pre_placebo <- preprocess(dat=coordinate, trisk=trisk, nrisk=nrisk, totalpts=13780, maxy=1) 
est_placebo <- getIPD(prep=pre_placebo,armID=0,tot.events=239)
hosp_ipd <- est_placebo$IPD
table(hosp_ipd$status)

####
##evolo
####
coordinate <- data1[,c(1,5)]
trisk <- data1$trisk[1:8]
nrisk <- data1$nrisk[1:8]
pre_evolo <- preprocess(dat=coordinate, trisk=trisk, nrisk=nrisk, totalpts=13784, maxy=1) 
est_evolo<- getIPD(prep=pre_evolo,armID=1,tot.events=236)
temp <- est_evolo$IPD
table(temp$status)
hosp_ipd <- rbind(hosp_ipd,temp)
# write_xlsx(hosp_ipd, "hosp_ipd.xlsx")

## check the plot
sf <- survfit(Surv(time, status) ~ treat, data = hosp_ipd)
hosp_summary <- summary(sf, times = trisk)

# Split at-risk numbers by group
at_risk_matrix <- matrix(mi_summary$n.risk, ncol = 2, byrow = FALSE)
colnames(at_risk_matrix) <- c("placebo","evolo")

# Plot cumulative incidence function
plot(sf, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of Hospitalization", xaxt = "n")


# Add x-axis labels at 6-month intervals
axis(1, at = trisk, labels = trisk)

# Add number at risk **above** the x-axis
text(x = trisk, y = -0.05 + (1 * 0.052), labels = at_risk_matrix[, 1],col = "red", cex = 0.6)
text(x = trisk, y = -0.05 + (1 * 0.05), labels = at_risk_matrix[, 2], col = "blue", cex = 0.6)
## Also: check the numbers of events in total 
table(hosp_ipd$status, hosp_ipd$treat)
at_risk_matrix
