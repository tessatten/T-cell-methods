library("devtools")
library("ggplot2")
library("arm")
library("vegan")
library("sets")
library("ineq")
library("tidyverse")
library("praise")

#import main data file
subsampling5000 = read.csv("subsampling5000.csv")  # read csv file

QUK_5000_info<- subset(subsampling5000, subset=(subsampling5000$Patient=="QUK"))
QUL_5000_info<- subset(subsampling5000, subset=(subsampling5000$Patient=="QUL"))
QUE_5000_info<- subset(subsampling5000, subset=(subsampling5000$Patient=="QUE"))


#QUK
#gini 
p = ggplot(QUK_5000_info, aes(x=QUK_5000_info$Week, y=QUK_5000_info$Gini.coefficient, color=QUK_5000_info$Chain, fill=QUK_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD8: Gini coefficient (library size 5000)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUK_CD8_5000_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUK_5000_info, aes(x=QUK_5000_info$Week, y=QUK_5000_info$Shannon.entropy, color=QUK_5000_info$Chain, fill=QUK_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUK CD8: Shannon index (library size 5000)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUK_CD8_5000_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUL
#gini 
p = ggplot(QUL_5000_info, aes(x=QUL_5000_info$Week, y=QUL_5000_info$Gini.coefficient, color=QUL_5000_info$Chain, fill=QUL_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD8: Gini coefficient (library size 5000)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUL_CD8_5000_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUL_5000_info, aes(x=QUL_5000_info$Week, y=QUL_5000_info$Shannon.entropy, color=QUL_5000_info$Chain, fill=QUL_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUL CD8: Shannon index (library size 5000)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUL_CD8_5000_shannon.pdf", width = 27, height = 20, units = c('cm'))

#QUE
#gini 
p = ggplot(QUE_5000_info, aes(x=QUE_5000_info$Week, y=QUE_5000_info$Gini.coefficient, color=QUE_5000_info$Chain, fill=QUE_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD8: Gini coefficient (library size 5000)", x = "Week", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUE_CD8_5000_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(QUE_5000_info, aes(x=QUE_5000_info$Week, y=QUE_5000_info$Shannon.entropy, color=QUE_5000_info$Chain, fill=QUE_5000_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  labs(title = "QUE CD8: Shannon index (library size 5000)", x = "Week", y = "Shannon entropy", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("QUE_CD8_5000_shannon.pdf", width = 27, height = 20, units = c('cm'))

# Gini Index

for (i in unique(cdr3Dat$ID)){
  tempC <- c()
  tempM <- c()
  tempP <- c()
  tempSimp <- c()
  tempShan <- c()
  for (j in unique(cdr3Dat$month[cdr3Dat$ID==i])){
    for (k in unique(cdr3Dat$chain[cdr3Dat$ID==i & cdr3Dat$month==j])){
      temp <- cdr3Dat$V2[cdr3Dat$chain==k & cdr3Dat$month==j & cdr3Dat$ID==i]
      tempC <- append(tempC, k)
      tempM <- append(tempM, j)
      tempP <- append(tempP, ineq(temp))
      tempSimp <- c(tempSimp, diversity(temp, index = 'simpson'))
      tempShan <- c(tempShan, diversity(temp, index = 'shannon', base = 2))
    }
  }
  # Gini
  temp <- data.frame(tempP, tempM, tempC)
  ggp <- ggplot(data = temp)+
    geom_point(aes(x=tempM, y=tempP, group=tempC, colour=tempC))+
    geom_line(aes(x=tempM, y=tempP, group=tempC, colour=tempC))+
    scale_y_continuous(limits=c(0,1))+
    labs(x='Month', y='Gini Coefficient', colour='Chain', title=i)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 1))
  
  ggsave(paste(divGiniOutputPath, '/', i, '.pdf', sep = ''), 
         ggp, width = 9, height = 7, units = c('cm'))
  
  # Simpson
  temp <- data.frame(tempSimp, tempM, tempC)
  ggp <- ggplot(data = temp)+
    geom_point(aes(x=tempM, y=tempSimp, group=tempC, colour=tempC))+
    geom_line(aes(x=tempM, y=tempSimp, group=tempC, colour=tempC))+
    labs(x='Month', y='Simpson Index', colour='Chain', title=i)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 1))
  
  ggsave(paste(divSimpsonOutputPath, '/', i, '.pdf', sep = ''), 
         ggp, width = 9, height = 7, units = c('cm'))
  
  # Shannon
  temp <- data.frame(tempShan, tempM, tempC)
  ggp <- ggplot(data = temp)+
    geom_point(aes(x=tempM, y=tempShan, group=tempC, colour=tempC))+
    geom_line(aes(x=tempM, y=tempShan, group=tempC, colour=tempC))+
    labs(x='Month', y='Shannon Entropy', colour='Chain', title=i)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 1))
  
  ggsave(paste(divShannonOutputPath, '/', i, '.pdf', sep = ''), 
         ggp, width = 9, height = 7, units = c('cm'))
  
  
}

