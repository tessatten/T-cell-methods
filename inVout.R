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

statsSheet <- read.csv(file = "TCR_statistics_Feb.csv") #import metadata for samples (samples that worked)

library(tidyr)
test =  statsSheet %>% drop_na(CD4.count)
test$Month.after.transplant2 <- factor(test$Month.after.transplant, levels = c("0.75", "1", "2", "3", "6", "12"))
statsSheet$Month.after.transplant2 <- factor(statsSheet$Month.after.transplant, levels = c("0.75", "1", "2", "3", "6", "12"))

#cell count vs read count
p = ggplot(test, aes(x=test$CD4.count, y=test$Number.of.error.corrected.reads, color=test$Month.after.transplant2, fill=test$Month.after.transplant)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "Clinical CD4 cell count vs. read count", x = "Clinical CD4 cell count ( x 10^9)", y = "Read count", color = "Month") +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ scale_y_continuous(labels = comma)+ geom_smooth(method = "lm", se = FALSE)
p
ggsave("CD4ellCountvsReadCountUCB.pdf", width = 27, height = 20, units = c('cm'))


#ginivShannon
p = ggplot(statsSheet, aes(x=statsSheet$Gini.coefficient, y=statsSheet$Shannon.entropy, color=statsSheet$Chain, fill=statsSheet$Chain, shape = statsSheet$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "", x = "Gini coefficient", y = "Shannon entropy", color = "Chain") +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ scale_y_continuous(labels = comma)+ geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))
p
ggsave("ShannonvGiniAll.pdf", width = 27, height = 20, units = c('cm'))


#ShannonvVJ
p1 = ggplot(statsSheet, aes(x=statsSheet$Number.of.VJ.rearrangements, y=statsSheet$Shannon.entropy, color=statsSheet$Chain, fill=statsSheet$Chain, shape = statsSheet$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "", x = "Number of identified reads", y = "Shannon entropy", color = "Chain") +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ scale_y_continuous(labels = comma)+ geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
p1
ggsave("ShannonvVJChain.pdf", width = 27, height = 20, units = c('cm'))


#ShannonvVJ
p2 = ggplot(statsSheet, aes(x=statsSheet$Number.of.error.corrected.reads, y=statsSheet$Shannon.entropy, color=statsSheet$Chain, fill=statsSheet$Chain, shape = statsSheet$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "", x = "Number of error-corrected reads", y = "Shannon entropy", color = "Chain") +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ scale_x_continuous(labels = comma)+ geom_smooth(method = "lm", se = FALSE) + guides(fill=FALSE, shape=FALSE)+
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#56B4E9", "#D55E00"))+ theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
p2
ggsave("ShannonvErrorCorrectChain.pdf", width = 27, height = 20, units = c('cm'))

p1 + p2+
  plot_layout(guides = 'collect')
ggsave("ShannonvReadCountChain.pdf", width = 30, height = 20, units = c('cm'))
