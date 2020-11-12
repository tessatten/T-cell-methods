library("ape")
library("ggplot2")
library("plyr")
library("scales")
library("vegan")
library("knitr")
library("dplyr")
library("praise")
library("tidyverse")
library("phyloseq")
library("Bios2cor")
library("grid")
library("ggpubr")
library("patchwork")

statsSheet <- read.csv(file = "pentaFactSheet.csv") #import metadata for samples (samples that worked)
tcrstats_new_march <- read.csv(file = "TCR_statistics_March.csv")

#doing a t test to see if population means are differnt
alpha <- subset(statsSheet, subset=(statsSheet$Chain=="alpha"))
beta <- subset(statsSheet, subset=(statsSheet$Chain=="beta"))

noUCBs <- subset(statsSheet, subset=(statsSheet$Month.after.transplant !="UCB"))
noUCBalpha <- subset(noUCBs, subset=(noUCBs$Chain=="alpha"))
noUCBbeta <- subset(noUCBs, subset=(noUCBs$Chain=="beta"))


qqplot(alpha$Percentage.of.rearrangements.found)
t.test(alpha$Percentage.of.rearrangements.found, beta$Percentage.of.rearrangements.found)

#########################################

#reads vs vdj reads in alpha
p1 = ggplot(noUCBalpha, aes(x=noUCBalpha$Number.of.reads, y=noUCBalpha$Number.of.VJ.rearrangements, color=noUCBalpha$Individual, fill=noUCBalpha$Individual, shape = noUCBalpha$Patient.or.control)) +
  geom_point(size=8) + labs(title = "Alpha Chain", x = "Read count from MiSeq", y = "Number of identified reads", color = "Individual") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + scale_y_continuous(labels = comma)+ scale_x_continuous(labels = comma)
p1
ggsave("VJvsReadsALPHA.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs vdj reads beta
p2 = ggplot(noUCBbeta, aes(x=noUCBbeta$Number.of.reads, y=noUCBbeta$Number.of.VJ.rearrangements, color=noUCBbeta$Individual, fill=noUCBbeta$Individual, shape = noUCBbeta$Patient.or.control)) +
  geom_point(size=8) + labs(title = "Beta chain ", x = "Read count from MiSeq", y = "Number of identified reads", color = "Individual") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE) + scale_x_continuous(labels = comma)+ scale_y_continuous(labels = comma)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
p2
ggsave("VJvsReadsBETA.pdf", width = 40, height = 30, units = c('cm'))
############################################
p1 + p2+
  plot_layout(guides = 'collect')
ggsave("VJvsReadsBETAtogether.pdf", width = 65, height = 45, units = c('cm'))
############################################

############################################
#number of vj vs error correct
ggplot(noUCBs, aes(x=noUCBs$Number.of.VJ.rearrangements, y=noUCBs$Number.of.error.corrected.reads, color=noUCBs$Chain, fill=noUCBs$Chain, shape = noUCBs$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Number of reads (rearrangement identified)", y = "Number of reads (error-corrected)", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00")) + scale_y_continuous(labels = comma) + scale_x_continuous(labels = comma)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))

ggsave("VJvsErrorCorrectReadChain.pdf", width = 40, height = 30, units = c('cm'))
############################################
#number of vj vs error correct ALPHA
ggplot(noUCBalpha, aes(x=noUCBalpha$Number.of.error.corrected.reads, y=noUCBalpha$Number.of.VJ.rearrangements, color=noUCBalpha$Individual, fill=noUCBalpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error-corrected)", y = "Number of reads (VJ)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsErrorCorrectReadALPHA.pdf", width = 40, height = 30, units = c('cm'))
############################################
#number of vj vs error correct BETA
ggplot(noUCBbeta, aes(x=noUCBbeta$Number.of.error.corrected.reads, y=noUCBbeta$Number.of.VJ.rearrangements, color=noUCBbeta$Individual, fill=noUCBbeta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error-corrected)", y = "Number of reads (VJ)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsErrorCorrectReadBETA.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#vdj reads vs shannon
p1 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$Number.of.VJ.rearrangements, y=tcrstats_new_march$Shannon.entropy, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Number of identified reads", y = "Shannon entropy", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))
p1
ggsave("VJvsShannon.pdf", width = 50, height = 40, units = c('cm'))
############################################

#########################################
#error-corrected reads vs shannon
p2 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$Number.of.error.corrected.reads, y=tcrstats_new_march$Shannon.entropy, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Number of error-corrected reads", y = "Shannon entropy", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ scale_x_continuous(labels = comma)+ geom_vline(xintercept=1000)
p2
ggsave("error-correctdvsShannon.pdf", width = 50, height = 40, units = c('cm'))

############################################
p1 + p2+
  plot_layout(guides = 'collect')
ggsave("ReadsvsShannontogether.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#vdj reads vs gini
p1 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$Number.of.VJ.rearrangements, y=tcrstats_new_march$Gini.coefficient, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Number of identified reads", y = "Gini coefficient", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))
p1
ggsave("VJvsGini.pdf", width = 50, height = 40, units = c('cm'))
############################################

#########################################
#error-corrected reads vs gini
p2 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$Number.of.error.corrected.reads, y=tcrstats_new_march$Gini.coefficient, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Number of error-corrected reads", y = "Gini coefficient", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ scale_x_continuous(labels = comma)
p2
ggsave("error-correctdvsGini.pdf", width = 50, height = 40, units = c('cm'))

############################################
p1 + p2+
  plot_layout(guides = 'collect')
ggsave("ReadsvsGiniTogether.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################
#ng RNA vs e-corected reads
p2 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$ng.of.RNA, y=tcrstats_new_march$Number.of.error.corrected.reads, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Nanograms of RNA", y = "Number of error-corrected reads", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ scale_y_continuous(labels = comma)
p2
ggsave("error-correctdvsNGRNA.pdf", width = 40, height = 30, units = c('cm'))

############################################

#########################################
#ng RNA vs e-corected reads
p2 = ggplot(tcrstats_new_march, aes(x=tcrstats_new_march$ng.of.RNA, y=tcrstats_new_march$Shannon.entropy, color=tcrstats_new_march$Chain, fill=tcrstats_new_march$Chain, shape = tcrstats_new_march$Chain)) +
  geom_point(size=8) + labs(title = "", x = "Nanograms of RNA", y = "Shannon entropy", color = "Chain") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=24, face="bold"),
    axis.title.y = element_text(color="black", size=24, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE, shape=FALSE)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))  + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ scale_y_continuous(labels = comma)
p2
ggsave("ShannonvsNGRNA.pdf", width = 40, height = 30, units = c('cm'))

############################################