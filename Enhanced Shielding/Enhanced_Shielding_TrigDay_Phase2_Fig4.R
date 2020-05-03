setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New") # This is where the plots Output
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("Cairo")

#This Script contains otuput for both the Trigger Day and Phase 2 Analysis - Part of Figure 4

#### Model Functions ####
#Generation Time Function - Uses Baseline R0 
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Beta Functions - Uses Two different Re - for beta1 and beta 4 + beta 2 and beta 3 
beta1 <- function(time, tstart1, tdur, r014) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (r014*(gamma))*0.5
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(r014*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         r014*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta1(seq(0,730), 71, (6*7), 0.6))

beta2 <- function(time, tstart1, tdur, r023) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - ((2.8*(gamma) - r023*(gamma))*0.5))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(r023*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         r023*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta2(seq(0,730), 71, (6*7), 0.7))

beta3 <- function(time, tstart1, tdur, r023) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - (2.8*(gamma) - 1.7*(gamma))*0.5)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(r023*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         r023*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta3(seq(0,730), 71, (6*7), 0.7))


beta4 <- function(time,tstart1,tdur, r014) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- r014*(gamma)*0.5
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(r014*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         r014*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta4(seq(0,730), 71, (6*7), 0.6))

#ODEs Enhanced Shielding Model - Adapted to use different phase 2 R_e
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- beta1(time,tstart1,tdur,r014)
    beta2 <- beta2(time,tstart1,tdur,r023)
    beta3 <- beta3(time,tstart1,tdur,r023)
    beta4 <- beta4(time,tstart1,tdur,r014)
    
    dSv = - beta1*Iv*Sv - beta1*Is*Sv - beta4*Ir1*Sv - beta4*Ir2*Sv - beta4*Ir3*Sv + zeta*Rv
    dSs = - beta1*Iv*Ss - beta1*Is*Ss - beta2*Ir1*Ss - beta2*Ir2*Ss - beta2*Ir3*Ss + zeta*Rs
    dSr1 = - beta4*Iv*Sr1 - beta2*Is*Sr1 - beta3*Ir1*Sr1 - beta3*Ir2*Sr1 - beta3*Ir3*Sr1 + zeta*Rr1
    dSr2 = - beta4*Iv*Sr2 - beta2*Is*Sr2 - beta3*Ir1*Sr2 - beta3*Ir2*Sr2 - beta3*Ir3*Sr2 + zeta*Rr2
    dSr3 = - beta4*Iv*Sr3 - beta2*Is*Sr3 - beta3*Ir1*Sr3 - beta3*Ir2*Sr3 - beta3*Ir3*Sr3 + zeta*Rr3
    
    dIv = beta1*Iv*Sv + beta1*Is*Sv + beta4*Ir1*Sv + beta4*Ir2*Sv + beta4*Ir3*Sv - gamma*Iv
    dIs = beta1*Iv*Ss + beta1*Is*Ss + beta2*Ir1*Ss + beta2*Ir2*Ss + beta2*Ir3*Ss - gamma*Is
    dIr1 = beta4*Iv*Sr1 + beta2*Is*Sr1 + beta3*Ir1*Sr1 + beta3*Ir2*Sr1 + beta3*Ir3*Sr1 - gamma*Ir1
    dIr2 = beta4*Iv*Sr2 + beta2*Is*Sr2 + beta3*Ir1*Sr2 + beta3*Ir2*Sr2 + beta3*Ir3*Sr2 - gamma*Ir2
    dIr3 = beta4*Iv*Sr3 + beta2*Is*Sr3 + beta3*Ir1*Sr3 + beta3*Ir2*Sr3 + beta3*Ir3*Sr3 - gamma*Ir3
    
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


#### Baseline Parameters - Initial Conditions, Times and Baseline Parameters ####

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0)

times <- seq(0, 478, by = 1)

parmsbase = c(gamma = 1/(GenTime(3.3,2.8)), 
              zeta = 1/365,
              tstart1 = 71, 
              tdur = 6*7) 

#### SET 1) TRIGGER DAY ####

trigday <- c(46, 71, 96)

trigdata <- data.frame(matrix(nrow = 0, ncol = 15))

for(i in 1:length(trigday)) {
  parms <- c(parmsbase, r014 = 0.8, r023 = 0.9)
  parms["tstart1"] <- trigday[i]
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
  out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
  out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
  out1 <- out1[, -grep("Sr|Ir|Rr", names(out1))]
  out1[,2:7] <- out1[,2:7]/0.2
  out1[,8:10] <- out1[,8:10]/0.6
  out1$Beta1 <- beta1(times, 71, (6*7), parms[["r014"]]) 
  out1$Beta2 <- beta2(times, 71, (6*7), parms[["r023"]]) 
  out1$Beta3 <- beta3(times, 71, (6*7), parms[["r023"]])
  out1$Beta4 <- beta4(times, 71, (6*7), parms[["r014"]])
  out1$trigday <- trigday[i]
  trigdata <- rbind(out1, trigdata)
}

colnames(trigdata) <- c("Time", "Suscv", "Suscs", "Infected_Iv", "Infected_Is", "Recovv", "Recovs", "RemSusc", "RemInf", "RemRecov",
                      "Beta_1", "Beta_2", "Beta_3", "Beta_4", "TriggerDay")

#Trigger Day 96

phase1 <- data.frame(xmin=0, xmax=96, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=96, xmax=96+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=96+(6*7), xmax=96+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=96+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

statsinfecv <- melt(trigdata[trigdata$TriggerDay == 96,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf96 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 36) +
  geom_text(data = phase2, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 7) +
  geom_text(data = phase3, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 27) +
  geom_text(data = phase4, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 110) +
  geom_line(size = 1.02, stat = "identity")

#Baseline Day 71

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

statsinfecv <- melt(trigdata[trigdata$TriggerDay == 71,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf71 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 11) +
  geom_text(data = phase2, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase4, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Trigger Day 46

phase1 <- data.frame(xmin=0, xmax=46, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=46, xmax=46+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=46+(6*7), xmax=46+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=46+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

statsinfecv <- melt(trigdata[trigdata$TriggerDay == 46,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf46 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 11) +
  geom_text(data = phase2, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase4, aes(x = xmin, y = 0.107, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#### SET 2) PHASE 2 ANALYSIS ####

retest <- data.frame(scenario = c("Lower", "Baseline", "Higher"),
                     Re14 = c(0.6, 0.8, 1),
                     Re23 = c(0.7, 0.9, 1.1))

redata <- data.frame(matrix(nrow = 0, ncol = 15))

for(i in 1:nrow(retest)) {
  parms <- c(parmsbase, r014 = retest[i,2], r023 = retest[i,3])
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
  out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
  out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
  out1 <- out1[, -grep("Sr|Ir|Rr", names(out1))]
  out1[,2:7] <- out1[,2:7]/0.2
  out1[,8:10] <- out1[,8:10]/0.6
  out1$Beta1 <- beta1(times, 71, (6*7), parms[["r014"]]) 
  out1$Beta2 <- beta2(times, 71, (6*7), parms[["r023"]]) 
  out1$Beta3 <- beta3(times, 71, (6*7), parms[["r023"]])
  out1$Beta4 <- beta4(times, 71, (6*7), parms[["r014"]])
  out1$r_e4 <- as.factor(retest[i,1])
  redata <- rbind(out1, redata)
}

colnames(redata) <- c("Time", "Suscv", "Suscs", "Infected_Iv", "Infected_Is", "Recovv", "Recovs", "RemSusc", "RemInf", "RemRecov",
                      "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Phase4_Re")

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

# R_e 0.6/0.7

statsinfecv <- melt(redata[redata$Phase4_Re == "Lower",], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf0607 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

# R_e 0.8/0.9

statsinfecv <- melt(redata[redata$Phase4_Re == "Baseline",], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf0809 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

# R_e 1.0/1.1

statsinfecv <- melt(redata[redata$Phase4_Re == "Higher",], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf1011 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.07,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#### Plotting - Phase 2 ####

plot <- ggarrange(pinf46,pinf0607,
                  
                  pinf71, pinf0809,
                  
                  pinf96, pinf1011,  
                  
                  ncol = 2, nrow = 3, 
                  
                  heights = c(1, 1, 
                              1.1), common.legend = TRUE, legend = "bottom",
                  labels = c("A", "B"), font.label = c(size = 20))

#Anti-Aliasing in Plots - Can Ignore 
ggsave(plot, filename = "Fig4.png", dpi = 300, type = "cairo",
       width = 14, height = 10, units = "in")
