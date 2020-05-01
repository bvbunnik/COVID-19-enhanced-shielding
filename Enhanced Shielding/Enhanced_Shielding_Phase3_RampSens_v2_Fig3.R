rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer"); library("sensitivity");library("fast")

#### Model Functions ####

#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

gamma <- 1/(GenTime(3.3,2.8))

#Function to model intervention - currently set at baseline - added additional functionality to it

beta1 <- function(time, tstart1, tdur, tdur2_14, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (0.8*(gamma))*scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(tdur2_14*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(tdur2_14*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(tdur2_14*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

beta2 <- function(time, tstart1, tdur, tdur2_23, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - ((2.8*(gamma) - 0.9*(gamma))*scaling))
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(tdur2_23*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(tdur2_23*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(tdur2_23*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

beta3 <- function(time, tstart1, tdur,tdur2_23, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- (2.8*(gamma) - (2.8*(gamma) - 1.7*(gamma))*scaling)
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(tdur2_23*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(tdur2_23*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(tdur2_23*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

beta4 <- function(time,tstart1,tdur, tdur2_14, scaling) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.8*(gamma) *scaling
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(tdur2_14*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(tdur2_14*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(tdur2_14*7) & time <= 730),
                       beta1_2,
                       1.7*(gamma))))}

#Function for Shielded/non-Shielded Pop
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- beta1(time,tstart1,tdur,tdur2_14,scaling)
    beta2 <- beta2(time,tstart1,tdur,tdur2_23,scaling)
    beta3 <- beta3(time,tstart1,tdur,tdur2_23,scaling)
    beta4 <- beta4(time,tstart1,tdur,tdur2_14,scaling)
    
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

#### Baseline Parameters + Phase Setting ####

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

#### Ramp Up 12/6/18 ####
#12
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 12,
          tdur2_23 = 12) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))
colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf12 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
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

#6
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 12,
          tdur2_23 = 6) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))

colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf6 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
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

#8
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 12,
          tdur2_23 = 18) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))

colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf18 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#GGARRANGE
p1 <- ggarrange(pinf12, NULL, pinf6, NULL, pinf18, nrow = 5, ncol = 1, align = "v", heights = c(0.5, -0.03, 0.5, -0.03, 0.6), labels = c("A"), font.label = c(size = 20) )

#### Ramp DOWN 12/6/18 ####
#12
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 12,
          tdur2_23 = 12) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))

out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))
colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf12down <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
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

#6
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 6,
          tdur2_23 = 12) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))

colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf6down <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
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

#8
parms = c(gamma = 1/(GenTime(3.3,2.8)), 
          zeta = 1/365,
          tstart1 = 71, 
          tdur = 6*7,
          scaling = 0.5,
          tdur2_14 = 18,
          tdur2_23 = 12) #CAN VARY THIS - DEPENDING ON FACTOR EXPLORED

out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms))
out1$RemS <- out1$Sr1 + out1$Sr2 + out1$Sr3
out1$RemI <- out1$Ir1 + out1$Ir2 + out1$Ir3
out1$RemR <- out1$Rr1 + out1$Rr2 + out1$Rr3
out1$Sv <- out1$Sv/0.20; out1$Ss <- out1$Ss/0.20; out1$RemS <- out1$RemS/0.60
out1$Iv <- out1$Iv/0.20; out1$Is <- out1$Is/0.20; out1$RemI <- out1$RemI/0.60
out1$Rv <- out1$Rv/0.20; out1$Rs <- out1$Rs/0.20; out1$RemR <- out1$RemR/0.60
out1$Beta1 <- beta1(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5])) 
out1$Beta2 <- beta2(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5])) 
out1$Beta3 <- beta3(times, 71, (6*7), as.numeric(parms[7]), as.numeric(parms[5]))
out1$Beta4 <- beta4(times, 71, (6*7), as.numeric(parms[6]), as.numeric(parms[5]))

colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr1", "Suscr2","Suscr3", 
                    "Infected_Iv", "Infected_Is", "Infected_Ir1", "Infected_Ir2", "Infected_Ir3", 
                    "Recovv", "Recovs", "Recovr1", "Recovr2", "Recovr2",
                    "RemSusc", "RemInf", "RemRecov",
                    "Beta_1", "Beta_2", "Beta_3", "Beta_4")
statsinfecv <- melt(out1, id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "RemInf"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

pinf18down <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.07), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("General", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.06, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#GGARRANGE
#p2 <- ggarrange(pinf12down, NULL, pinf6down, NULL, pinf18down, nrow = 5, ncol = 1, align = "v", heights = c(0.5, -0.03, 0.5, -0.03, 0.6), labels = c("B"), font.label = c(size = 20) )
#p1 <- ggarrange(pinf12, NULL, pinf6, NULL, pinf18, nrow = 5, ncol = 1, align = "v", heights = c(0.5, -0.03, 0.5, -0.03, 0.6), labels = c("A"), font.label = c(size = 20) )

plot <- ggarrange(pinf12, pinf12down,
          
          pinf6,  pinf6down,
          
          pinf18, pinf18down, ncol = 2, nrow = 3, 
          heights = c(1, 1, 
                      1.25), common.legend = TRUE, legend = "bottom",
          labels = c("A", "B"), font.label = c(size = 20))

#Anti-Aliasing in Plots - Can Ignore 
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New")
ggsave(plot, filename = "Fig3.png", dpi = 300, type = "cairo",
       width = 14, height = 10, units = "in")

