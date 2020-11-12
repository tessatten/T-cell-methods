library("devtools")
library("ggplot2")
library("arm")
library("vegan")
library("sets")
library("ineq")
library("tidyverse")
library("praise")

#import main data file
subsamplingIndividuals = read.csv("subsamplingIndividual.csv")  # read csv file

CD4_subsample<- subset(subsamplingIndividuals, subset=(subsamplingIndividuals$Cell_type=="CD4"))
CD8_subsample<- subset(subsamplingIndividuals, subset=(subsamplingIndividuals$Cell_type=="CD8"))

QUD_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUD"))
QUE_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUE"))
QUF_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUF"))
QUK_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUK"))
QUL_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUL"))
QUP_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUP"))
QUQ_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUQ"))
QUT_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUT"))
QUV_CD4_info<- subset(CD4_subsample, subset=(CD4_subsample$Patient=="QUV"))

QUE_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUE"))
QUF_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUF"))
QUK_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUK"))
QUL_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUL"))
QUQ_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUQ"))
QUT_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUT"))
QUV_CD8_info<- subset(CD8_subsample, subset=(CD8_subsample$Patient=="QUV"))


#QUD
#gini 
p = ggplot(QUD_CD4_info, aes(x=QUD_CD4_info$Week, y=QUD_CD4_info$Gini.coefficient, color=QUD_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUD CD4: Gini coefficient (library size 2,224)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUD_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUD_CD4_info, aes(x=QUD_CD4_info$Week, y=QUD_CD4_info$Shannon.entropy, color=QUD_CD4_info$Chain, fill=QUD_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUD CD4: Shannon index (library size 2,224)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUD_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUE
#gini 
p = ggplot(QUE_CD4_info, aes(x=QUE_CD4_info$Week, y=QUE_CD4_info$Gini.coefficient, color=QUE_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD4: Gini coefficient (library size 12,346)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUE_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUE_CD4_info, aes(x=QUE_CD4_info$Week, y=QUE_CD4_info$Shannon.entropy, color=QUE_CD4_info$Chain, fill=QUE_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD4: Shannon index (library size 12,346)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUE_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUF
#gini 
p = ggplot(QUF_CD4_info, aes(x=QUF_CD4_info$Week, y=QUF_CD4_info$Gini.coefficient, color=QUF_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUF CD4: Gini coefficient (library size 1,011)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUF_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUF_CD4_info, aes(x=QUF_CD4_info$Week, y=QUF_CD4_info$Shannon.entropy, color=QUF_CD4_info$Chain, fill=QUF_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUF CD4: Shannon index (library size 1,011)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUF_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUK
#gini 
p = ggplot(QUK_CD4_info, aes(x=QUK_CD4_info$Week, y=QUK_CD4_info$Gini.coefficient, color=QUK_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD4: Gini coefficient (library size 3,583)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUK_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUK_CD4_info, aes(x=QUK_CD4_info$Week, y=QUK_CD4_info$Shannon.entropy, color=QUK_CD4_info$Chain, fill=QUK_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD4: Shannon index (library size 3,583)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUK_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUL
#gini 
p = ggplot(QUL_CD4_info, aes(x=QUL_CD4_info$Week, y=QUL_CD4_info$Gini.coefficient, color=QUL_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD4: Gini coefficient (library size 2,319)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUL_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUL_CD4_info, aes(x=QUL_CD4_info$Week, y=QUL_CD4_info$Shannon.entropy, color=QUL_CD4_info$Chain, fill=QUL_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD4: Shannon index (library size 2,319)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUL_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUP
#gini 
p = ggplot(QUP_CD4_info, aes(x=QUP_CD4_info$Week, y=QUP_CD4_info$Gini.coefficient, color=QUP_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUP CD4: Gini coefficient (library size 3,083)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUP_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUP_CD4_info, aes(x=QUP_CD4_info$Week, y=QUP_CD4_info$Shannon.entropy, color=QUP_CD4_info$Chain, fill=QUP_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUP CD4: Shannon index (library size 3,083)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUP_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUQ
#gini 
p = ggplot(QUQ_CD4_info, aes(x=QUQ_CD4_info$Week, y=QUQ_CD4_info$Gini.coefficient, color=QUQ_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUQ CD4: Gini coefficient (library size 1,729)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUQ_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUQ_CD4_info, aes(x=QUQ_CD4_info$Week, y=QUQ_CD4_info$Shannon.entropy, color=QUQ_CD4_info$Chain, fill=QUQ_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUQ CD4: Shannon index (library size 1,729)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUQ_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUT
#gini 
p = ggplot(QUT_CD4_info, aes(x=QUT_CD4_info$Week, y=QUT_CD4_info$Gini.coefficient, color=QUT_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUT CD4: Gini coefficient (library size 2,387)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUT_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUT_CD4_info, aes(x=QUT_CD4_info$Week, y=QUT_CD4_info$Shannon.entropy, color=QUT_CD4_info$Chain, fill=QUT_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUT CD4: Shannon index (library size 2,387)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUT_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUV
#gini 
p = ggplot(QUV_CD4_info, aes(x=QUV_CD4_info$Week, y=QUV_CD4_info$Gini.coefficient, color=QUV_CD4_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUV CD4: Gini coefficient (library size 2,546)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUV_CD4_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUV_CD4_info, aes(x=QUV_CD4_info$Week, y=QUV_CD4_info$Shannon.entropy, color=QUV_CD4_info$Chain, fill=QUV_CD4_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUV CD4: Shannon index (library size 2,546)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUV_CD4_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))


#QUE
#gini 
p = ggplot(QUE_CD8_info, aes(x=QUE_CD8_info$Week, y=QUE_CD8_info$Gini.coefficient, color=QUE_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD8: Gini coefficient (library size 5,991)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUE_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUE_CD8_info, aes(x=QUE_CD8_info$Week, y=QUE_CD8_info$Shannon.entropy, color=QUE_CD8_info$Chain, fill=QUE_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD8: Shannon index (library size 5,991)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUE_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUF
#gini 
p = ggplot(QUF_CD8_info, aes(x=QUF_CD8_info$Week, y=QUF_CD8_info$Gini.coefficient, color=QUF_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUF CD8: Gini coefficient (library size 1,942)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUF_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUF_CD8_info, aes(x=QUF_CD8_info$Week, y=QUF_CD8_info$Shannon.entropy, color=QUF_CD8_info$Chain, fill=QUF_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUF CD8: Shannon index (library size 1,942)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUF_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUK
#gini 
p = ggplot(QUK_CD8_info, aes(x=QUK_CD8_info$Week, y=QUK_CD8_info$Gini.coefficient, color=QUK_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD8: Gini coefficient (library size 13,035)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUK_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUK_CD8_info, aes(x=QUK_CD8_info$Week, y=QUK_CD8_info$Shannon.entropy, color=QUK_CD8_info$Chain, fill=QUK_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD8: Shannon index (library size 13,035)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUK_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUL
#gini 
p = ggplot(QUL_CD8_info, aes(x=QUL_CD8_info$Week, y=QUL_CD8_info$Gini.coefficient, color=QUL_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD8: Gini coefficient (library size 5,614)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUL_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUL_CD8_info, aes(x=QUL_CD8_info$Week, y=QUL_CD8_info$Shannon.entropy, color=QUL_CD8_info$Chain, fill=QUL_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD8: Shannon index (library size 5,614)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUL_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUQ
#gini 
p = ggplot(QUQ_CD8_info, aes(x=QUQ_CD8_info$Week, y=QUQ_CD8_info$Gini.coefficient, color=QUQ_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUQ CD8: Gini coefficient (library size 2,762)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUQ_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUQ_CD8_info, aes(x=QUQ_CD8_info$Week, y=QUQ_CD8_info$Shannon.entropy, color=QUQ_CD8_info$Chain, fill=QUQ_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUQ CD8: Shannon index (library size 2,762)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUQ_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUT
#gini 
p = ggplot(QUT_CD8_info, aes(x=QUT_CD8_info$Week, y=QUT_CD8_info$Gini.coefficient, color=QUT_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUT CD8: Gini coefficient (library size 1,102)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUT_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUT_CD8_info, aes(x=QUT_CD8_info$Week, y=QUT_CD8_info$Shannon.entropy, color=QUT_CD8_info$Chain, fill=QUT_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUT CD8: Shannon index (library size 1,102)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUT_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUV
#gini 
p = ggplot(QUV_CD8_info, aes(x=QUV_CD8_info$Week, y=QUV_CD8_info$Gini.coefficient, color=QUV_CD8_info$Chain, fill = "black")) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUV CD8: Gini coefficient (library size 1,156)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2) + coord_cartesian(ylim = c(0, 1))
p
ggsave("QUV_CD8_individ_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUV_CD8_info, aes(x=QUV_CD8_info$Week, y=QUV_CD8_info$Shannon.entropy, color=QUV_CD8_info$Chain, fill=QUV_CD8_info$Chain)) +
  geom_point(size = 10) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUV CD8: Shannon index (library size 1,156)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+ geom_line(
    stat = "identity",
    position = "identity",
    inherit.aes = TRUE,size = 2)+ coord_cartesian(ylim = c(5, 16))
p
ggsave("QUV_CD8_individ_shannon.pdf", width = 27, height = 20, units = c('cm'))




