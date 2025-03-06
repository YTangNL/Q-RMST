library(IPDfromKM)
library(haven)
library(survival)
library(survminer)
library(writexl)
library(readxl)
####
#Revascularization
####
data0 <- read_excel("/Users/miyake/Library/CloudStorage/GoogleDrive-yixiaotang0626@gmail.com/My Drive/LUMC Master Thesis/Data/coordinate/revasc/placebo_revasc.xlsx")
data1 <- read_excel("/Users/miyake/Library/CloudStorage/GoogleDrive-yixiaotang0626@gmail.com/My Drive/LUMC Master Thesis/Data/coordinate/revasc/evolo_revasc.xlsx")

####
##placebo
####
coordinate <- data0[,c(1,5)]
trisk <- data0$trisk[1:8]
nrisk <- data0$nrisk[1:8]
pre_placebo <- preprocess(dat=coordinate, trisk=trisk, nrisk=nrisk, totalpts=13780, maxy=1) 
est_placebo <- getIPD(prep=pre_placebo,armID=0,tot.events=965)
revasc_ipd <- est_placebo$IPD
table(revasc_ipd$status)

####
##evolo
####
coordinate <- data1[,c(1,5)]
trisk <- data1$trisk[1:8]
nrisk <- data1$nrisk[1:8]
pre_evolo <- preprocess(dat=coordinate, trisk=trisk, nrisk=nrisk, totalpts=13784, maxy=1) 
est_evolo<- getIPD(prep=pre_evolo,armID=1,tot.events=759)
temp <- est_evolo$IPD
table(temp$status)
revasc_ipd <- rbind(revasc_ipd,temp)
#write_xlsx(revasc_ipd, "revasc_ipd.xlsx")

## check the plot
sf <- survfit(Surv(time, status) ~ treat, data = revasc_ipd)
revasc_summary <- summary(sf, times = trisk)

# Split at-risk numbers by group
at_risk_matrix <- matrix(revasc_summary$n.risk, ncol = 2, byrow = FALSE)
colnames(at_risk_matrix) <- c("placebo","evolo")

# Plot cumulative incidence function
plot(sf, fun = "event", col = c("red", "blue"), lwd = 1,
     xlab = "Time (months)", ylab = "Cumulative Incidence",
     main = "Incident Event Curve of Revasc", xaxt = "n")


# Add x-axis labels at 6-month intervals
axis(1, at = trisk, labels = trisk)

# Add number at risk **above** the x-axis
text(x = trisk, y = -0.05 + (1 * 0.057), labels = at_risk_matrix[, 1],col = "red", cex = 0.6)
text(x = trisk, y = -0.05 + (1 * 0.05), labels = at_risk_matrix[, 2], col = "blue", cex = 0.6)
## Also: check the numbers of events in total 
table(stroke_ipd$status, stroke_ipd$treat)
at_risk_matrix
