setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New") # This is where the plots Output
rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("Cairo")

#### Model Functions  ####
#Function for the generation time/(1/gamma) parameter
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

#Function to model interventions - Beta1,2,3 and 4
beta1 <- function(time, tstart1, tdur, phase1R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*gamma
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       phase1R0*(gamma))))}

plot(beta1(seq(0,730), 71, (6*7), 0.5,1.7), ylim = c(0,0.5))

beta2 <- function(time, tstart1, tdur,phase1R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 1.85*gamma
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       phase1R0*(gamma))))}

plot(beta2(seq(0,730), 71, (6*7), 0.5,1.7))

beta3 <- function(time, tstart1, tdur,phase1R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 2.25*gamma
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.9*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.9*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       phase1R0*(gamma))))}

plot(beta3(seq(0,730), 71, (6*7), 0.5, 1.7))


beta4 <- function(time,tstart1,tdur, phase1R0) {
  gamma <- 1/(GenTime(3.3,2.8))
  beta1_2 <- 0.4*gamma
  betalin <- approxfun(x=c(tstart1+tdur, tstart1+tdur+(12*7)), y = c(0.8*(gamma), beta1_2), method="linear", rule  =2)
  ifelse((time >= tstart1 & time <= tstart1+tdur), #Phase 2
         0.8*(gamma),
         ifelse((time >= tstart1+tdur & time <= tstart1+tdur+(12*7)), #Phase 3
                betalin(time),
                ifelse((time >= tstart1+tdur+(12*7) & time <= 730),
                       beta1_2,
                       phase1R0*(gamma))))}

plot(beta4(seq(0,730), 71, (6*7), 0.5,1.7))

#Function for ODEs
SIRS <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta1 <- beta1(time,tstart1,tdur,phase1R0)
    beta2 <- beta2(time,tstart1,tdur,phase1R0)
    beta3 <- beta3(time,tstart1,tdur,phase1R0)
    beta4 <- beta4(time,tstart1,tdur,phase1R0)
    
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

#### Common Model Parameters - Initial Conditions and Timings ####

init <- c(Sv = 0.2 - 0.0001*0.2, Ss = 0.2 - 0.0001*0.2, 
          Sr1 = 0.2 - 0.0001*0.2, Sr2 = 0.2 - 0.0001*0.2, Sr3 = 0.2 - 0.0001*0.2,
          Iv = 0.0001*0.2, Is = 0.0001*0.2, Ir1 = 0.0001*0.2, Ir2 = 0.0001*0.2, Ir3 = 0.0001*0.2,   
          Rv= 0, Rs = 0, Rr1 = 0, Rr2 = 0, Rr3 = 0)
times <- seq(0, 478, by = 1)

#### Testing the Zeta Parameter #### 

#We are identifying the point at which the 2nd peak is higher than the 1st peak 
zetaseq <- seq(1,365)

output <- data.frame(matrix(ncol = 7, nrow = length(zetaseq)))
colnames(output) <- c("DayImmune","TimeSecPeak","HeightSecPeak","TimeFirPeak","HeightFirPeak","HigherFirPeak", "RelHeight1stvs2nd")

#Run the for model for different zetas
for (i in 1:length(zetaseq)) {
  temp <- numeric(7)
  parms1 = c(gamma = 1/(GenTime(3.3,2.8)), 
             zeta = 1/zetaseq[i],
             tstart1 = 71, 
             tdur = 6*7,
             phase1R0 = 1.7)
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms1))
  out1$Iv <- out1$Iv/0.20
  temp[1] <- zetaseq[i]
  temp[2] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][2]
  temp[3] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][2]
  temp[4] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][1]
  temp[5] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][1]
  temp[6] <- ifelse((out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)][2] > out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)][1]), 1, 0)
  temp[7] <- temp[3]/temp[5]
  output[i,] <- temp
  print(i/length(zetaseq))
}

#Plot comparing 2nd Peak vs 1st Peak height
statzeta <- melt(output, id.vars = c("DayImmune"), measure.vars = c("RelHeight1stvs2nd"))

p2 <- ggplot(statzeta, aes(x = DayImmune, y = value)) + geom_line(size = 1.02, stat = "identity", col = "red") + theme_bw() +
  labs(x ="Duration of Immunity (Days)", y = "Relative Height of 2nd Peak vs 1st Peak") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, lty = 2, size = 1.02, col = "black") +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Testing Phase 1 R0 ####

phase1R0seq <- seq(1.4, 2, by = 0.005)
output <- data.frame(matrix(ncol = 6, nrow = length(phase1R0seq)))
colnames(output) <- c("Phase1R0","TimeSecPeak","Second Peak","TimeFirPeak","First Peak","RelHeight1stvs2nd")

#Only assessing vulnerable population

for (i in 1:length(phase1R0seq)) {
  temp <- numeric(6)
  parms1 = c(gamma = 1/(GenTime(3.3,2.8)), 
             zeta = 1/365,
             tstart1 = 71, 
             tdur = 6*7,
             scaling = 0.5,
             phase1R0 = phase1R0seq[i])
  out1 <- data.frame(ode(y = init, func = SIRS, times = times, parms = parms1))
  out1$Iv <- out1$Iv/0.20
  temp[1] <- phase1R0seq[i]
  temp[2] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure Time of 2nd Peak
  temp[3] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][2] #Outcome Measure I(t)  of 2nd Peak
  temp[4] <- out1$time[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure Time of 1st Peak
  temp[5] <- out1$Iv[which(diff(sign(diff(out1$Iv)))==-2)+1][1] #Outcome Measure I(t) of 1st Peak
  temp[6] <- temp[3]/temp[5]
  output[i,] <- temp
  print(i/length(phase1R0seq))
}

statr0 <- melt(output, id.vars = c("Phase1R0"), measure.vars = c("RelHeight1stvs2nd"))

p1 <- ggplot(statr0, aes(x = Phase1R0, y = value)) + geom_line(size = 1.02, stat = "identity", col = "red") + theme_bw() +
  labs(x =expression("Phase 1 R"[e]), y = "Relative Height of 2nd Peak vs 1st Peak") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 	1, lty = 2, size = 1.02, col = "black") +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Plotting for Fig 5 #### 

plot <- ggarrange(p1, NULL, p2, nrow = 3, ncol = 1, align = "v", heights = c(0.5, -0.03, 0.5), 
                  labels = c("A","", "B"), font.label = c(size = 20) )

#Anti-Aliasing in Plots - Can Ignore 
ggsave(plot, filename = "Fig5.png", dpi = 300, type = "cairo",
       width = 8.5, height = 11, units = "in")