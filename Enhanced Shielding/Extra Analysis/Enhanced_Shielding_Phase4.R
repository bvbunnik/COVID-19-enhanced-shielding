setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New") # This is where the plots Output
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("Cairo")

#### Model Functions - Generation Time + Betas + ODEs ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Betas for the Intervention - Phase 4 is a scaling factor which affects all Betas during Phase 4 equally 
beta1 <- function(time, tstart1, tdur, phase4) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*gamma*phase4
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta1(seq(0,730), 71, (6*7), 0.5), ylim = c(0,0.5))

beta2 <- function(time, tstart1, tdur, phase4) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 1.85*gamma*phase4
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta2(seq(0,730), 71, (6*7), 0.75))

beta3 <- function(time, tstart1, tdur, phase4) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 2.25*gamma*phase4
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta3(seq(0,730), 71, (6*7), 0.5))

beta4 <- function(time,tstart1,tdur, phase4) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*gamma*phase4
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

plot(beta4(seq(0,730), 71, (6*7), 0.75))

#ODEs integrating new betas in
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- beta1(time,tstart1,tdur, phase4)
    beta2 <- beta2(time,tstart1,tdur, phase4)
    beta3 <- beta3(time,tstart1,tdur, phase4)
    beta4 <- beta4(time,tstart1,tdur, phase4)
    
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

#### Baseline - Phase4 Scaling Factor 1 ####

#Initial Conditions and Times

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0)

times <- seq(0, 478, by = 1)

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

#### Running the Model For Phase 4 Analysis ####

phase4seq <- seq(0.75, 1.25, by = 0.25)

phase4data <- data.frame(matrix(nrow = 0, ncol = 15))

for(i in 1:length(phase4seq)) {
  parms = c(gamma = 1/(GenTime(3.3,2.8)), 
            zeta = 1/365,
            tstart1 = 71, 
            tdur = 6*7,
            phase4 = phase4seq[i]) 
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
  out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
  out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
  out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
  out1 <- out1[, -grep("Sr|Ir|Rr", names(out1))]
  out1[,2:7] <- out1[,2:7]/0.2
  out1[,8:10] <- out1[,8:10]/0.6
  out1$Beta1 <- beta1(times, 71, (6*7), parms[[5]]) 
  out1$Beta2 <- beta2(times, 71, (6*7), parms[[5]]) 
  out1$Beta3 <- beta3(times, 71, (6*7), parms[[5]])
  out1$Beta4 <- beta4(times, 71, (6*7), parms[[5]])
  out1$phase4 <- phase4seq[i]
  phase4data <- rbind(out1, phase4data)
}

colnames(phase4data) <- c("Time", "Suscv", "Suscs", "Infected_Iv", "Infected_Is", "Recovv", "Recovs", "RemSusc", "RemInf", "RemRecov",
                        "Beta_1", "Beta_2", "Beta_3", "Beta_4", "Phase4Mod")

# Baseline Phase 4 - 1.0
statsinfecv <- melt(phase4data[phase4data$Phase4Mod == 1,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinfbase <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
 geom_line(size = 1.02, stat = "identity")

# Phase 4 - 1.25

statsinfecv <- melt(phase4data[phase4data$Phase4Mod == 1.25,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf125 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

# Phase 4 - 0.75

statsinfecv <- melt(phase4data[phase4data$Phase4Mod == 0.75,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf075 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.110, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#### Plotting ####

plot <- ggarrange(pinf075, NULL, pinfbase,NULL, pinf125, nrow = 5, ncol =1, align = "v", heights = c(0.5, -0.03, 0.5,-0.05, 0.6), 
                  labels = c("A", "", "B", "", "C"),font.label = c(size = 20))

plotlim <- ggarrange(pinf075, NULL, pinf125, nrow = 3, ncol =1, align = "v", heights = c(0.5, -0.03, 0.6), 
                  labels = c("A", "", "B"),font.label = c(size = 20))

#Anti-Aliasing in Plots - Can Ignore 
ggsave(plot, filename = "Fig10.png", dpi = 300, type = "cairo",
       width = 8.5, height = 11, units = "in")

ggsave(plotlim, filename = "Fig10_2Limit.png", dpi = 300, type = "cairo",
       width = 8.5, height = 11, units = "in")

#### Ranging Phase 4 #### 

phase4scale <- seq(0.75, 1.25, by = 0.01)

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0)
times <- seq(0, 478, by = 1)

output <- data.frame(matrix(ncol = 10, nrow = length(phase4scale)))
colnames(output) <- c("Phase4BetaScale", "Beta1","Beta2", "Beta3", "Beta4","TimeSecPeak","Second Peak",
                      "TimeFirPeak","First Peak","RelHeight1stvs2nd")

for (i in 1:length(phase4scale)) {
  temp <- numeric(10)
  parms1 = c(gamma = 1/(GenTime(3.3,2.8)), 
             zeta = 1/365,
             tstart1 = 71, 
             tdur = 6*7,
             phase4 = phase4scale[i])
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms1))
  out1$Iv <- out1$Iv/0.20
  temp[1] <- phase4scale[i]
  temp[2] <- beta1(400, 71, (6*7), phase4scale[i])
  temp[3] <- beta2(400, 71, (6*7), phase4scale[i])
  temp[4] <- beta3(400, 71, (6*7), phase4scale[i])
  temp[5] <- beta4(400, 71, (6*7), phase4scale[i])
  temp[6] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure Time of 2nd Peak
  temp[7] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure I(t)  of 2nd Peak
  temp[8] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure Time of 1st Peak
  temp[9] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure I(t) of 1st Peak
  temp[10] <- temp[7]/temp[9]
  output[i,] <- temp
  print(i/length(phase4scale))
}

statr0 <- melt(output, id.vars = c("Phase4BetaScale"), measure.vars = c("RelHeight1stvs2nd"))

p1 <- ggplot(statr0, aes(x = Phase4BetaScale, y = value)) + geom_line(size = 1.02, stat = "identity", col = "red") + theme_bw() +
  labs(x =expression("Phase 4 Scaling"), y = "Relative Height of 2nd Peak vs 1st Peak") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 	1, lty = 2, size = 1.02, col = "black") +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Plotting for Fig 5 #### 

plotrange <- ggarrange(p1, nrow = 1, ncol = 1, align = "v", heights = c(0.5), font.label = c(size = 20) )

#Anti-Aliasing in Plots - Can Ignore 

ggsave(plotrange, filename = "Fig10_Phase4Range.png", dpi = 300, type = "cairo",
       width = 8, height = 7, units = "in")
