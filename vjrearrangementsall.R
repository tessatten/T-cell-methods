tcrstats = read.csv("TCR basic statistics.csv")  # read csv file
tcrstats$Month.after.transplant2 <- factor(tcrstats$Month.after.transplant, levels = c("1", "2", "3", "6", "8", "12", "18", "22"))

# based on variable values
alpha <- subset(tcrstats, subset=(tcrstats$Chain=="alpha" & tcrstats$Month.after.transplant > 0.5))
beta <- subset(tcrstats, subset=(tcrstats$Chain=="beta"& tcrstats$Month.after.transplant > 0.5))

#by individual
RYF <- subset(tcrstats, subset=(tcrstats$Individual=="RYF"))
AAF <- subset(tcrstats, subset=(tcrstats$Individual=="AAF"))
DSF <- subset(tcrstats, subset=(tcrstats$Individual=="DSF"))
YOM <- subset(tcrstats, subset=(tcrstats$Individual=="YOM"))
HBF <- subset(tcrstats, subset=(tcrstats$Individual=="HBF"))
M <- subset(tcrstats, subset=(tcrstats$Individual=="M"))

A <- subset(tcrstats, subset=(tcrstats$Individual=="A"))
C <- subset(tcrstats, subset=(tcrstats$Individual=="C"))
E <- subset(tcrstats, subset=(tcrstats$Individual=="E"))
L <- subset(tcrstats, subset=(tcrstats$Individual=="L"))

#alpha
RYFalpha <- subset(alpha, subset=(alpha$Individual=="RYF"))
AAFalpha <- subset(alpha, subset=(alpha$Individual=="AAF"))
DSFalpha <- subset(alpha, subset=(alpha$Individual=="DSF"))
YOMalpha <- subset(alpha, subset=(alpha$Individual=="YOM"))
HBFalpha <- subset(alpha, subset=(alpha$Individual=="HBF"))
Malpha <- subset(alpha, subset=(alpha$Individual=="M"))

Aalpha <- subset(alpha, subset=(alpha$Individual=="A"))
Calpha <- subset(alpha, subset=(alpha$Individual=="C"))
Ealpha <- subset(alpha, subset=(alpha$Individual=="E"))
Lalpha <- subset(alpha, subset=(alpha$Individual=="L"))

#beta
RYFbeta <- subset(beta, subset=(beta$Individual=="RYF"))
AAFbeta <- subset(beta, subset=(beta$Individual=="AAF"))
DSFbeta <- subset(beta, subset=(beta$Individual=="DSF"))
YOMbeta <- subset(beta, subset=(beta$Individual=="YOM"))
HBFbeta <- subset(beta, subset=(beta$Individual=="HBF"))
Mbeta <- subset(beta, subset=(beta$Individual=="M"))

Abeta <- subset(beta, subset=(beta$Individual=="A"))
Cbeta <- subset(beta, subset=(beta$Individual=="C"))
Ebeta <- subset(beta, subset=(beta$Individual=="E"))
Lbeta <- subset(beta, subset=(beta$Individual=="L"))

# total and unique sequences

##########################################
#number of vj rearrangement reads in all samples
ggplot(tcrstats_new, aes(x=tcrstats_new$Sample.ID, y=tcrstats_new$Number.of.VJ.rearrangements, color=tcrstats_new$Individual, fill=tcrstats_new$Individual,las=2)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Reads (VJ rearrangements)") +
  guides(fill=FALSE) +theme_bw()

#number of vj rearrangement reads in alpha samples
ggplot(alpha, aes(x=alpha$Sample.ID, y=alpha$Number.of.VJ.rearrangements, color=alpha$Individual, fill=alpha$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Reads (VJ rearrangements)") +
  guides(fill=FALSE)

#number of vj rearrangement reads in beta samples
ggplot(beta, aes(x=beta$Sample.ID, y=beta$Number.of.VJ.rearrangements, color=beta$Individual, fill=beta$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Reads (VJ rearrangements)") +
  guides(fill=FALSE)

########################################
#number of clonotypes in all samples
ggplot(tcrstats_new, aes(x=tcrstats_new$Sample.ID, y=tcrstats_new$Number.of.sequences.observed, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Clonotypes") +
  guides(fill=FALSE)

#########################################
#number of unique clonotypes in all samples
ggplot(tcrstats_new, aes(x=tcrstats_new$Sample.ID, y=tcrstats_new$Number.of.unique.sequences, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Sample ID") + 
  ylab("Number of Unique Clonotypes") +
  guides(fill=FALSE)+ theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=8))

##########################################
#size of largest clonal expansion in all samples
ggplot(tcrstats, aes(x=tcrstats$Sample.ID, y=tcrstats$Largest.clonal.expansion, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Size of Largest Clonal Expansion") +
  guides(fill=FALSE)

#size of largest clonal expansion in all samples
ggplot(tcrstats, aes(x=tcrstats$Month.after.transplant2, y=tcrstats$Largest.clonal.expansion, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Size of Largest Clonal Expansion") +
  guides(fill=FALSE)

##########################################
# Manual month levels
tcrstats$Month.after.transplant2 <- factor(tcrstats$Month.after.transplant, levels = c("1", "2", "3", "6", "8", "12", "18", "22"))
############################################
#size of largest clonal expansion in all samples
ggplot(tcrstats, aes(x=tcrstats$Month.after.transplant2, y=tcrstats$Largest.clonal.expansion, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Size of Largest Clonal Expansion") +
  guides(fill=FALSE)
############################################
#size of largest clonal expansion in alpha samples ALPHA
ggplot(alpha, aes(x=alpha$Month.after.transplant2, y=alpha$Largest.clonal.expansion, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Size of Largest Clonal Expansion") +
  guides(fill=FALSE)
############################################
#size of largest clonal expansion in beta samples BETA
ggplot(beta, aes(x=alpha$Month.after.transplant2, y=beta$Largest.clonal.expansion, color=beta$Individual, fill=beta$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Size of Largest Clonal Expansion") +
  guides(fill=FALSE)

##############################################
#number of unique clonotypes in all samples - barchart
ggplot(tcrstats, aes(x=tcrstats$Sample.ID, y=tcrstats$Number.of.unique.sequences, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Unique Clonotypes") +
  guides(fill=FALSE)

#number of unique clonotypes in all samples - scatterplot by month
ggplot(tcrstats, aes(x=tcrstats$Month.after.transplant2, y=tcrstats$Number.of.unique.sequences, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Number of Unique Clonotypes") +
  guides(fill=FALSE)

##############################################
#number of clonotypes in all samples - barchart
ggplot(tcrstats, aes(x=tcrstats$Sample.ID, y=tcrstats$Number.of.different.sequences, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Number of Clonotypes") +
  guides(fill=FALSE)


#number of clonotypes in all samples - scatterplot by month
ggplot(tcrstats, aes(x=tcrstats$Month.after.transplant2, y=tcrstats$Number.of.different.clonotypes, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Number of Different Clonotypes") +
  guides(fill=FALSE)

tcrstats$Month.after.transplant <- as.character(tcrstats$Month.after.transplant)