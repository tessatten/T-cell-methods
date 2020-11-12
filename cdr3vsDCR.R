

library("devtools")
library("ggplot2")
library("vegan")
library("ineq")
library("tidyverse")
library("praise")

#import main data file
tcrstats_new = read.csv("TCR_statistics_March.csv")  # read csv file
tcrstats_new$Month.after.transplant2 <- factor(tcrstats_new$Month.after.transplant, levels = c("0.5", "0.75", "0", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

alpha <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="alpha"))
beta <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="beta"))
alpha$Month.after.transplant2 <- factor(alpha$Month.after.transplant, levels = c("0.5", "0.75", "0", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))
beta$Month.after.transplant2 <- factor(beta$Month.after.transplant, levels = c("0.5", "0.75", "0", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

clinical = subset(tcrstats_new, subset=(tcrstats_new$Patient.or.control =="Patient" & tcrstats_new$CD4.count > 0.01))
clinical$Month.after.transplant3 <- factor(clinical$Month.after.transplant, levels = c("0.75", "1", "2", "3", "6", "12"))

alphapatients = subset(alpha, subset=(alpha$Patient.or.control=="Patient"))
betapatients = subset(beta, subset=(beta$Patient.or.control=="Patient"))

############################################3

#in this file am comparing the difference betwen cdr3 and dcr data

#########################################

#shannon in cdr3 and dcr
ggplot(tcrstats_new, aes(x=tcrstats_new$Shannon.entropy.CDR3, y=tcrstats_new$Shannon.entropy, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_point(size=16) + labs(title = " ", x = "Shannon index: CDR3", y = "Shannon index: DCR", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonCDR3vsDCR.pdf", width = 40, height = 30, units = c('cm'))
############################################
#shannon in cdr3 and dcr by CHAIN
p1 = ggplot(tcrstats_new, aes(x=tcrstats_new$Shannon.entropy.CDR3, y=tcrstats_new$Shannon.entropy, color=tcrstats_new$Chain, fill=tcrstats_new$Chain, shape = tcrstats_new$Chain)) +
  geom_point(size=8) + labs(title = " ", x = "Shannon index: CDR3", y = "Shannon index: DCR", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))
p1
ggsave("ShannonCDR3vsDCRCHAIN.pdf", width = 40, height = 30, units = c('cm'))
############################################
#gini in cdr3 and dcr by CHAIN
p2 = ggplot(tcrstats_new, aes(x=tcrstats_new$Gini.coefficient.CDR3, y=tcrstats_new$Gini.coefficient, color=tcrstats_new$Chain, fill=tcrstats_new$Chain, shape = tcrstats_new$Chain)) +
  geom_point(size=8) + labs(title = " ", x = "Gini coefficient: CDR3", y = "Gini coefficient: DCR", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))
p2
ggsave("GiniCDR3vsDCRCHAIN.pdf", width = 40, height = 30, units = c('cm'))
#
############################################
############################################
p1 + p2+
  plot_layout(guides = 'collect')
ggsave("diversityTogetherCDR3vDCR.pdf", width = 30, height = 20, units = c('cm'))
############################################
# Shapiro-Wilk normality test for wt
shapiro.test(tcrstats_new$Shannon.entropy) # => p = 0.09
shapiro.test(tcrstats_new$Shannon.entropy.CDR3) # => p = 0.09
spearman <-cor.test(tcrstats_new$Shannon.entropy, tcrstats_new$Shannon.entropy.CDR3,  method = "spearman")
spearman
#########################################

#gini in cdr3 and dcr
ggplot(tcrstats_new, aes(x=tcrstats_new$Gini.coefficient.CDR3, y=tcrstats_new$Gini.coefficient, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_point(size=16) + labs(title = " ", x = "Gini coefficient: CDR3", y = "Gini coefficient: DCR", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GiniCDR3vsDCR.pdf", width = 40, height = 30, units = c('cm'))
############################################
# Shapiro-Wilk normality test for wt
shapiro.test(tcrstats_new$Gini.coefficient) # => p = 0.09
shapiro.test(tcrstats_new$Gini.coefficient.CDR3) # => p = 0.09
spearman <-cor.test(tcrstats_new$Gini.coefficient, tcrstats_new$Gini.coefficient.CDR3,  method = "spearman")
spearman