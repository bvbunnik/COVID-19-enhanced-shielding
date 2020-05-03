rm(list=ls())
library("deSolve"); library("ggplot2"); library("ggpubr"); library("reshape2"); library("dplyr"); library("RColorBrewer"); library("sensitivity");library("fast")
library("Cairo")

#### 2/2/96 ####
outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_2-2-96.csv")
outimp <- outimp[,1:(ncol(outimp)-3)]
out <- outimp[, -grep( "cum" , colnames(outimp) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = out[,grepl( "Sv" , names(out) )]/0.02,
                   "SuscS" = out[,grepl( "Sh" , names(out) )]/0.02,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.96,
                   "InfV" = out[,grepl( "Iv" , names(out) )]/0.02,
                   "InfS" = out[,grepl( "Ih" , names(out) )]/0.02,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.96,
                   "RecovV" = out[,grepl( "Rv" , names(out) )]/0.02,
                   "RecovS" = out[,grepl( "Rh" , names(out) )]/0.02,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.96)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf96<- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.05,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

pinf96solo <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.05,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New")
ggsave(pinf96solo, filename = "FigS2.png", dpi = 300, type = "cairo",
       width = 8, height = 7, units = "in")

#Cumulative Calc
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_2-2-96.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                   "cumV" = outc[,grepl( "Iv" , names(outc) )]/0.02,
                   "cumE" = (outc[,grepl( "Ih" , names(outc) )] + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.98)

out1c$cumV.cumIv[out1c$time == 113+(365)]-out1c$cumV.cumIv[out1c$time == 113]

out1c$cumE[out1c$time == 113+(365)]- out1c$cumE[out1c$time == 113]

#### 10/10/30 ####
outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_10-10-30_test.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = rowSums(out[,grepl( "Sv" , names(out) )])/0.2,
                   "SuscS" = rowSums(out[,grepl( "Sh" , names(out) )])/0.2,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.6,
                   "InfV" = rowSums(out[,grepl( "Iv" , names(out) )])/0.2,
                   "InfS" = rowSums(out[,grepl( "Ih" , names(out) )])/0.2,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.6,
                   "RecovV" = rowSums(out[,grepl( "Rv" , names(out) )])/0.2,
                   "RecovS" = rowSums(out[,grepl( "Rh" , names(out) )])/0.2,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.6)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf60 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Cumulative Calc
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_10-10-30_test.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                    "cumV" = rowSums(outc[,grepl( "Iv" , names(outc) )])/0.20,
                    "cumE" = (rowSums(outc[,grepl( "Ih" , names(outc) )]) + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.80)

out1c$cumV[out1c$time == 113+(365)]-out1c$cumV[out1c$time == 113]

out1c$cumE[out1c$time == 113+(365)]- out1c$cumE[out1c$time == 113]

#### 8/8/84 ####
outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_8-8-84.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = rowSums(out[,grepl( "Sv" , names(out) )])/0.08,
                   "SuscS" = rowSums(out[,grepl( "Sh" , names(out) )])/0.08,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.84,
                   "InfV" = rowSums(out[,grepl( "Iv" , names(out) )])/0.08,
                   "InfS" = rowSums(out[,grepl( "Ih" , names(out) )])/0.08,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.84,
                   "RecovV" = rowSums(out[,grepl( "Rv" , names(out) )])/0.08,
                   "RecovS" = rowSums(out[,grepl( "Rh" , names(out) )])/0.08,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.84)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf84 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,-0.1,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Cumulative Calculations
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_8-8-84.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                    "cumV" = rowSums(outc[,grepl( "Iv" , names(outc) )])/0.08,
                    "cumE" = (rowSums(outc[,grepl( "Ih" , names(outc) )]) + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.92)

out1c$cumV[out1c$time == 113+(365)]-out1c$cumV[out1c$time == 113]

out1c$cumE[out1c$time == 113+(365)]- out1c$cumE[out1c$time == 113]

#### 14/14/72 ####
outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_14-14-72.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = rowSums(out[,grepl( "Sv" , names(out) )])/0.14,
                   "SuscS" = rowSums(out[,grepl( "Sh" , names(out) )])/0.14,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.72,
                   "InfV" = rowSums(out[,grepl( "Iv" , names(out) )])/0.14,
                   "InfS" = rowSums(out[,grepl( "Ih" , names(out) )])/0.14,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.72,
                   "RecovV" = rowSums(out[,grepl( "Rv" , names(out) )])/0.14,
                   "RecovS" = rowSums(out[,grepl( "Rh" , names(out) )])/0.14,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.72)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf72 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.125), expand = c(0,0)) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.11, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Cumulative Calculations
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_14-14-72.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                    "cumV" = rowSums(outc[,grepl( "Iv" , names(outc) )])/0.14,
                    "cumE" = (rowSums(outc[,grepl( "Ih" , names(outc) )]) + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.86)

out1c$cumV[out1c$time == 113+(365)]-out1c$cumV[out1c$time == 113]

out1c$cumE[out1c$time == 113+(365)]- out1c$cumE[out1c$time == 113]


#### Plotting ####
#GGARRANGE

ggarrange(pinf96, NULL, pinf72, nrow = 3, ncol =1, align = "v", heights = c(0.5, -0.03, 0.6), labels = c("A","", "C",""), font.label = c(size = 20) )
ggarrange(pinf84, NULL, pinf60, nrow = 3, ncol =1, align = "v", heights = c(0.5, -0.03, 0.6), labels = c("B","", "D"), font.label = c(size = 20) )

plot <- ggarrange(pinf60, pinf72,
                  
                  pinf84,  pinf96,
                  
                  ncol = 2, nrow = 2, 
                  heights = c(1, 1), common.legend = TRUE, legend = "bottom",
                  labels = c("A", "B", "C", "D"), font.label = c(size = 20))

#Anti-Aliasing in Plots - Can Ignore 
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New")
ggsave(plot, filename = "FigS3_v1.png", dpi = 300, type = "cairo",
       width = 14, height = 10, units = "in")

#### 10/10/30 ####
outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_10-10-30_test.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = rowSums(out[,grepl( "Sv" , names(out) )])/0.2,
                   "SuscS" = rowSums(out[,grepl( "Sh" , names(out) )])/0.2,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.6,
                   "InfV" = rowSums(out[,grepl( "Iv" , names(out) )])/0.2,
                   "InfS" = rowSums(out[,grepl( "Ih" , names(out) )])/0.2,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.6,
                   "RecovV" = rowSums(out[,grepl( "Rv" , names(out) )])/0.2,
                   "RecovS" = rowSums(out[,grepl( "Rh" , names(out) )])/0.2,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.6)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf60 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.075), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#### 20-40-40 ####

outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_20-40-40.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]
out <- out[, -grep( "p" , colnames(out) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = out[,grepl( "Sv" , names(out) )]/0.2,
                   "SuscS" = rowSums(out[,grepl( "Sh" , names(out) )])/0.4,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.4,
                   "InfV" = out[,grepl( "Iv" , names(out) )]/0.2,
                   "InfS" = rowSums(out[,grepl( "Ih" , names(out) )])/0.4,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.4,
                   "RecovV" = out[,grepl( "Rv" , names(out) )]/0.2,
                   "RecovS" = rowSums(out[,grepl( "Rh" , names(out) )])/0.4,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.4)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf204040 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.075), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_blank(), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Cumulative Calculations
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_14-14-72.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                    "cumV" = rowSums(outc[,grepl( "Iv" , names(outc) )])/0.14,
                    "cumE" = (rowSums(outc[,grepl( "Ih" , names(outc) )]) + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.86)

out1c$cumV[out1c$time == 71+(365)]-out1c$cumV[out1c$time == 71]

out1c$cumE[out1c$time == 71+(365)]- out1c$cumE[out1c$time == 71]

#### 20-10-70 ####

outimp <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_20-10-70.csv")
out <- outimp[, -grep( "cum" , colnames(outimp) )]
out <- out[, -grep( "p" , colnames(out) )]

out1 <- data.frame("time" = out$t, 
                   "SuscV" = rowSums(out[,grepl( "Sv" , names(out) )])/0.2,
                   "SuscS" = out[,grepl( "Sh" , names(out) )]/0.1,
                   "SuscR" = rowSums(out[,grepl( "Sr" , names(out) )])/0.7,
                   "InfV" = rowSums(out[,grepl( "Iv" , names(out) )])/0.2,
                   "InfS" = out[,grepl( "Ih" , names(out) )]/0.1,
                   "InfR" = rowSums(out[,grepl( "Ir" , names(out) )])/0.7,
                   "RecovV" = rowSums(out[,grepl( "Rv" , names(out) )])/0.2,
                   "RecovS" = out[,grepl( "Rh" , names(out) )]/0.1,
                   "RecovR" = rowSums(out[,grepl( "Rr" , names(out) )])/0.7)


colnames(out1) <- c("Time", "Suscv", "Suscs", "Suscr", "Infected_Iv", "Infected_Is", "Infected_Ir", "Recovv", "Recovs", "Recovr")

statsinfecv <- melt(out1[1:483,], id.vars = c("Time"), measure.vars = c("Infected_Iv", "Infected_Is", "Infected_Ir"))
statsinfecv$variable <- factor(statsinfecv$variable, levels=rev(levels(statsinfecv$variable)))

phase1 <- data.frame(xmin=0, xmax=71, ymin=-Inf, ymax=Inf, name = "P1")
phase2 <- data.frame(xmin=71, xmax=71+(6*7), ymin=-Inf, ymax=Inf, name = "P2")
phase3 <- data.frame(xmin=71+(6*7), xmax=71+(6*7)+(12*7), ymin=-Inf, ymax=Inf, name = "P3")
phase4 <- data.frame(xmin=71+(6*7)+(12*7), xmax=Inf, ymin=-Inf, ymax=Inf, name = "P4")

pinf201070 <- ggplot(data = statsinfecv, aes(x = (Time), y = value, col = variable))  + theme_bw() +
  labs(x ="Time (Days)", y = "Proportion Infected", color = "Population") + scale_y_continuous(limits = c(0,0.075), expand = c(0,0)) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14), 
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) + scale_x_continuous(expand = c(0, 0)) + 
  scale_color_manual(values=c('blue','darkgreen','red'), labels= c("Remainders", "Shielders", "Vulnerable")) +
  geom_rect(data=phase2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_rect(data=phase3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.35, inherit.aes = FALSE) +
  geom_rect(data=phase4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.5, inherit.aes = FALSE) +
  geom_text(data = phase1, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 25) +
  geom_text(data = phase2, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 8) +
  geom_text(data = phase3, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 30) +
  geom_text(data = phase4, aes(x = xmin, y = 0.065, label = name),inherit.aes = FALSE, size = 8, vjust = 0, hjust = 0, nudge_x = 125) +
  geom_line(size = 1.02, stat = "identity")

#Cumulative Calculations
outcum <- read.csv("C:/Users/amorg/Documents/PhD/nCoV Work/Models/SIRS_MComp_14-14-72.csv")
outc <- outcum[, grep( "cum" , colnames(outcum) )]

out1c <- data.frame("time" = outcum$t, 
                    "cumV" = rowSums(outc[,grepl( "Iv" , names(outc) )])/0.14,
                    "cumE" = (rowSums(outc[,grepl( "Ih" , names(outc) )]) + rowSums(outc[,grepl( "Ir" , names(outc) )]))/0.86)

out1c$cumV[out1c$time == 71+(365)]-out1c$cumV[out1c$time == 71]

out1c$cumE[out1c$time == 71+(365)]- out1c$cumE[out1c$time == 71]

#### Plotting ####
plot <- ggarrange(pinf204040, NULL, pinf60,NULL, pinf201070, 
                  ncol = 1, nrow = 5, 
                  heights = c(0.5, -0.03, 0.5,-0.05, 0.6), common.legend = TRUE, legend = "bottom",
                  labels = c("A","", "B","", "C"), font.label = c(size = 20))

#Anti-Aliasing in Plots - Can Ignore 
setwd("C:/Users/amorg/Documents/PhD/nCoV Work/Figures/Enhanced Shielding/New")
ggsave(plot, filename = "PopStruct.png", dpi = 300, type = "cairo",
       width = 8.5, height = 11,units = "in")

