setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New") # This is where the plots Output
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo")

#### Model Functions - Generation Time + Betas + ODEs ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Beta Functions - Added Compliance Scaling Factor which alters Betas TO and FROM the vulnerables
betacomp <- function(time, tstart1, tdur, comp) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (0.4+((1.7-0.4)*(1-comp)))*(gamma)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(betacomp(seq(0,730), 71, (6*7), 1), ylim = c(0,0.5))

beta1 <- function(time, tstart1, tdur) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*(gamma)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}
plot(beta1(seq(0,730), 71, (6*7)), ylim = c(0,0.5))

beta2 <- function(time, tstart1, tdur) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 1.85*(gamma)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}
plot(beta2(seq(0,730), 71, (6*7)))

beta3 <- function(time, tstart1, tdur) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 2.25*(gamma)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}
plot(beta3(seq(0,730), 71, (6*7)))

beta4 <- function(time,tstart1,tdur) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*gamma
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}
plot(beta4(seq(0,730), 71, (6*7)))

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    betacomp <- betacomp(time,tstart1,tdur,comp)
    beta1 <- beta1(time,tstart1,tdur)
    beta2 <- beta2(time,tstart1,tdur)
    beta3 <- beta3(time,tstart1,tdur)
    beta4 <- beta4(time,tstart1,tdur)
    
    dSv = - betacomp*Iv*Sv - betacomp*Is*Sv - betacomp*Ir1*Sv - betacomp*Ir2*Sv - betacomp*Ir3*Sv + zeta*Rv
    dSs = - betacomp*Iv*Ss - beta1*Is*Ss - beta2*Ir1*Ss - beta2*Ir2*Ss - beta2*Ir3*Ss + zeta*Rs
    dSr1 = - betacomp*Iv*Sr1 - beta2*Is*Sr1 - beta3*Ir1*Sr1 - beta3*Ir2*Sr1 - beta3*Ir3*Sr1 + zeta*Rr1
    dSr2 = - betacomp*Iv*Sr2 - beta2*Is*Sr2 - beta3*Ir1*Sr2 - beta3*Ir2*Sr2 - beta3*Ir3*Sr2 + zeta*Rr2
    dSr3 = - betacomp*Iv*Sr3 - beta2*Is*Sr3 - beta3*Ir1*Sr3 - beta3*Ir2*Sr3 - beta3*Ir3*Sr3 + zeta*Rr3
    
    dIv = betacomp*Iv*Sv + betacomp*Is*Sv + betacomp*Ir1*Sv + betacomp*Ir2*Sv + betacomp*Ir3*Sv - gamma*Iv
    dIs = betacomp*Iv*Ss + beta1*Is*Ss + beta2*Ir1*Ss + beta2*Ir2*Ss + beta2*Ir3*Ss - gamma*Is
    dIr1 = betacomp*Iv*Sr1 + beta2*Is*Sr1 + beta3*Ir1*Sr1 + beta3*Ir2*Sr1 + beta3*Ir3*Sr1 - gamma*Ir1
    dIr2 = betacomp*Iv*Sr2 + beta2*Is*Sr2 + beta3*Ir1*Sr2 + beta3*Ir2*Sr2 + beta3*Ir3*Sr2 - gamma*Ir2
    dIr3 = betacomp*Iv*Sr3 + beta2*Is*Sr3 + beta3*Ir1*Sr3 + beta3*Ir2*Sr3 + beta3*Ir3*Sr3 - gamma*Ir3
    
    dRv = gamma*Iv - zeta*Rv 
    dRs = gamma*Is - zeta*Rs
    dRr1 = gamma*Ir1 - zeta*Rr1
    dRr2 = gamma*Ir2 - zeta*Rr2
    dRr3 = gamma*Ir3 - zeta*Rr3

    return(list(c(dSv, dSs, dSr1, dSr2, dSr3,
                  dIv, dIs, dIr1, dIr2, dIr3,
                  dRv, dRs, dRr1, dRr2, dRr3)))
  })
}

#### Initial Common Parameters - Initial Conditions and Times ####
init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0)
times <- seq(0, 478, by = 1)

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

#### Running Simulations for 100, 75%, 50%, 25% and 0% Compliance

compliance <- seq(0,1,by = 0.25)

compdata <- data.frame(matrix(nrow = 0, ncol = 15))

for(i in 1:length(compliance)) {
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
                zeta = 1/365,
                tstart1 = 71, 
                tdur = 6*7,
                comp = compliance[i]) 
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
  out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
  out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
  out1 <- out1[, -grep("Sr|Ir|Rr", names(out1))]
  out1[,2:7] <- out1[,2:7]/0.2
  out1[,8:10] <- out1[,8:10]/0.6
  out1$Beta1 <- beta1(times, 71, (6*7)) 
  out1$Beta2 <- beta2(times, 71, (6*7)) 
  out1$Beta3 <- beta3(times, 71, (6*7))
  out1$Beta4 <- beta4(times, 71, (6*7))
  out1$compliance <- compliance[i]
  compdata <- rbind(out1, compdata)
}

colnames(compdata) <- c("Time", "Suscv", "Suscs", "Infected_Iv", "Infected_Is", "Recovv", "Recovs", "RemSusc", "RemInf", "RemRecov",
                        "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Compliance")

#Compliance 100%
statsinfecv <- melt(compdata[compdata$Compliance == 1,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf100 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 20) +
  geom_text(data = phase2, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 6) +
  geom_text(data = phase3, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 28) +
  geom_text(data = phase4, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 120) +
  geom_line(size = 1.02, stat = "identity")

#Compliance 75%
statsinfecv <- melt(compdata[compdata$Compliance == 0.75,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf75 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 20) +
  geom_text(data = phase2, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 6) +
  geom_text(data = phase3, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 28) +
  geom_text(data = phase4, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 120) +
  geom_line(size = 1.02, stat = "identity")

#Compliance 50%
statsinfecv <- melt(compdata[compdata$Compliance == 0.50,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf50 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 20) +
  geom_text(data = phase2, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 6) +
  geom_text(data = phase3, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 28) +
  geom_text(data = phase4, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 120) +
  geom_line(size = 1.02, stat = "identity")

#Compliance 25%
statsinfecv <- melt(compdata[compdata$Compliance == 0.25,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf25 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.3,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")


#Compliance 0%
statsinfecv <- melt(compdata[compdata$Compliance == 0,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf0 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.3,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 20) +
  geom_text(data = phase2, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 6) +
  geom_text(data = phase3, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 28) +
  geom_text(data = phase4, aes(x = xmin, y = 0.100, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 120) +
  geom_line(size = 1.02, stat = "identity")

#### Plotting of Compliance ####

plot <- ggarrange(pinf100, pinf75, pinf50, pinf25, pinf0, nrow = 3, ncol = 2, heights = c(0.5), font.label = c(size = 20),
                  labels = c("A", "B", "C", "D", "E"), common.legend = TRUE, legend = "bottom")

#### Compliance Range ####

compliance <- seq(0,1, by = 0.01)

output <- data.frame(matrix(ncol = 6, nrow = length(compliance)))
colnames(output) <- c("Compliance","TimeSecPeak","Second Peak","TimeFirPeak","First Peak","RelHeight1stvs2nd")

for (i in 1:length(compliance)) {
  temp <- numeric(6)
  parms1 = c(gamma = 1/(GenTime(3.3,2.8)), 
             zeta = 1/365,
             tstart1 = 71, 
             tdur = 6*7,
             comp = compliance[i])
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms1))
  out1$Iv <- out1$Iv/0.20
  temp[1] <- compliance[i]*100
  temp[2] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure Time of 2nd Peak
  temp[3] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure I(t)  of 2nd Peak
  temp[4] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure Time of 1st Peak
  temp[5] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure I(t) of 1st Peak
  temp[6] <- temp[3]/temp[5]
  output[i,] <- temp
  print(i/length(compliance))
}

statcomp <- melt(output, id.vars = c("Compliance"), measure.vars = c("RelHeight1stvs2nd"))

#### Plotting Compliance Range ####
p1 <- ggplot(statcomp, aes(x = Compliance, y = value)) + geom_line(size = 1.02, stat = "identity", col = "red") + theme_bw() +
  labs(x =expression("Compliance of the Vulnerable Population (%)"), y = "Relative Height of 2nd Peak vs 1st Peak") + scale_y_continuous(limits = c(0,3),expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 	1, lty = 2, size = 1.02, col = "black") +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))


plot1 <- ggarrange(p1, nrow = 1, ncol = 1, align = "v", heights = c(0.5), font.label = c(size = 20) )

#Anti-Aliasing in Plots - Can Ignore

ggsave(plot, filename = "Fig9.png", dpi = 300, type = "cairo",
       width = 14, height = 10, units = "in")

ggsave(plot1, filename = "Fig9_Range.png", dpi = 300, type = "cairo",
       width = 8, height = 7, units = "in")