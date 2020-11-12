library("devtools")
library("ggplot2")
library("arm")
library("vegan")
library("sets")
library("ineq")
library("tidyverse")
library("praise")


#####CD4
###QUD
QUDCD4a0total = read.csv("dcr_alpha_QUDCTW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUDCD4a0total, "names") <- c("CDR3", "Count")

QUDCD4a150total = read.csv("dcr_alpha_QUDCTW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUDCD4a150total, "names") <- c("CDR3", "Count")

QUDCD4b0total = read.csv("dcr_beta_QUDCTW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUDCD4b0total, "names") <- c("CDR3", "Count")

QUDCD4b150total = read.csv("dcr_beta_QUDCTW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUDCD4b150total, "names") <- c("CDR3", "Count")

###QUE

QUECD4a0total = read.csv("dcr_alpha_QUECTW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD4a0total, "names") <- c("CDR3", "Count")

QUECD4a150total = read.csv("dcr_alpha_QUECTW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD4a150total, "names") <- c("CDR3", "Count")

QUECD4b0total = read.csv("dcr_beta_QUECTW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD4b0total, "names") <- c("CDR3", "Count")

QUECD4b150total = read.csv("dcr_beta_QUECTW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD4b150total, "names") <- c("CDR3", "Count")

###QUF

QUFCD4a48total = read.csv("dcr_alpha_QUFCTW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4a48total, "names") <- c("CDR3", "Count")

QUFCD4b48total = read.csv("dcr_beta_QUFCTW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4b48total, "names") <- c("CDR3", "Count")

QUFCD4a0total = read.csv("dcr_alpha_QUFCTW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4a0total, "names") <- c("CDR3", "Count")

QUFCD4a150total = read.csv("dcr_alpha_QUFCTW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4a150total, "names") <- c("CDR3", "Count")

QUFCD4b0total = read.csv("dcr_beta_QUFCTW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4b0total, "names") <- c("CDR3", "Count")

QUFCD4b150total = read.csv("dcr_beta_QUFCTW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD4b150total, "names") <- c("CDR3", "Count")

###QUK

QUKCD4a48total = read.csv("dcr_alpha_QUKPTIW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4a48total, "names") <- c("CDR3", "Count")

QUKCD4b48total = read.csv("dcr_beta_QUKPTIW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4b48total, "names") <- c("CDR3", "Count")

QUKCD4a12total = read.csv("dcr_alpha_QUKPTIW12-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4a12total, "names") <- c("CDR3", "Count")

QUKCD4b12total = read.csv("dcr_beta_QUKPTIW12-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4b12total, "names") <- c("CDR3", "Count")

QUKCD4a0total = read.csv("dcr_alpha_QUKPTIW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4a0total, "names") <- c("CDR3", "Count")

QUKCD4a150total = read.csv("dcr_alpha_QUKPTIW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4a150total, "names") <- c("CDR3", "Count")

QUKCD4b0total = read.csv("dcr_beta_QUKPTIW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4b0total, "names") <- c("CDR3", "Count")

QUKCD4b150total = read.csv("dcr_beta_QUKPTIW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD4b150total, "names") <- c("CDR3", "Count")

######QUL

QULCD4a48total = read.csv("dcr_alpha_QULPTIW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4a48total, "names") <- c("CDR3", "Count")

QULCD4b48total = read.csv("dcr_beta_QULPTIW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4b48total, "names") <- c("CDR3", "Count")

QULCD4a12total = read.csv("dcr_alpha_QULPTIW12-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4a12total, "names") <- c("CDR3", "Count")

QULCD4b12total = read.csv("dcr_beta_QULPTIW12-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4b12total, "names") <- c("CDR3", "Count")

QULCD4a0total = read.csv("dcr_alpha_QULPTIW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4a0total, "names") <- c("CDR3", "Count")

QULCD4a150total = read.csv("dcr_alpha_QULPTIW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4a150total, "names") <- c("CDR3", "Count")

QULCD4b0total = read.csv("dcr_beta_QULPTIW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4b0total, "names") <- c("CDR3", "Count")

QULCD4b150total = read.csv("dcr_beta_QULPTIW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD4b150total, "names") <- c("CDR3", "Count")

######QUP

QUPCD4a0total = read.csv("dcr_alpha_QUPCTW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUPCD4a0total, "names") <- c("CDR3", "Count")

QUPCD4a150total = read.csv("dcr_alpha_QUPCTW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUPCD4a150total, "names") <- c("CDR3", "Count")

QUPCD4b0total = read.csv("dcr_beta_QUPCTW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUPCD4b0total, "names") <- c("CDR3", "Count")

QUPCD4b150total = read.csv("dcr_beta_QUPCTW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUPCD4b150total, "names") <- c("CDR3", "Count")


######QUQ

QUQCD4a48total = read.csv("dcr_alpha_QUQPTIW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4a48total, "names") <- c("CDR3", "Count")

QUQCD4b48total = read.csv("dcr_beta_QUQPTIW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4b48total, "names") <- c("CDR3", "Count")

QUQCD4a12total = read.csv("dcr_alpha_QUQPTIW12-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4a12total, "names") <- c("CDR3", "Count")

QUQCD4b12total = read.csv("dcr_beta_QUQPTIW12-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4b12total, "names") <- c("CDR3", "Count")

QUQCD4a0total = read.csv("dcr_alpha_QUQPTIW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4a0total, "names") <- c("CDR3", "Count")

QUQCD4a150total = read.csv("dcr_alpha_QUQPTIW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4a150total, "names") <- c("CDR3", "Count")

QUQCD4b0total = read.csv("dcr_beta_QUQPTIW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4b0total, "names") <- c("CDR3", "Count")

QUQCD4b150total = read.csv("dcr_beta_QUQPTIW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD4b150total, "names") <- c("CDR3", "Count")

######QUT

QUTCD4a48total = read.csv("dcr_alpha_QUTPTIW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4a48total, "names") <- c("CDR3", "Count")

QUTCD4b48total = read.csv("dcr_beta_QUTPTIW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4b48total, "names") <- c("CDR3", "Count")

QUTCD4a12total = read.csv("dcr_alpha_QUTPTIW12-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4a12total, "names") <- c("CDR3", "Count")

QUTCD4b12total = read.csv("dcr_beta_QUTPTIW12-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4b12total, "names") <- c("CDR3", "Count")

QUTCD4a0total = read.csv("dcr_alpha_QUTPTIW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4a0total, "names") <- c("CDR3", "Count")

QUTCD4a150total = read.csv("dcr_alpha_QUTPTIW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4a150total, "names") <- c("CDR3", "Count")

QUTCD4b0total = read.csv("dcr_beta_QUTPTIW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4b0total, "names") <- c("CDR3", "Count")

QUTCD4b150total = read.csv("dcr_beta_QUTPTIW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD4b150total, "names") <- c("CDR3", "Count")

######QUV

QUVCD4a48total = read.csv("dcr_alpha_QUVCTW48-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4a48total, "names") <- c("CDR3", "Count")

QUVCD4b48total = read.csv("dcr_beta_QUVCTW48-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4b48total, "names") <- c("CDR3", "Count")

QUVCD4a0total = read.csv("dcr_alpha_QUVCTW0-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4a0total, "names") <- c("CDR3", "Count")

QUVCD4a150total = read.csv("dcr_alpha_QUVCTW150-aCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4a150total, "names") <- c("CDR3", "Count")

QUVCD4b0total = read.csv("dcr_beta_QUVCTW0-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4b0total, "names") <- c("CDR3", "Count")

QUVCD4b150total = read.csv("dcr_beta_QUVCTW150-bCD4_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD4b150total, "names") <- c("CDR3", "Count")

#####CD8

###QUE

QUECD8a0total = read.csv("dcr_alpha_QUFCTW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD8a0total, "names") <- c("CDR3", "Count")

QUECD8a150total = read.csv("dcr_alpha_QUECTW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD8a150total, "names") <- c("CDR3", "Count")

QUECD8b0total = read.csv("dcr_beta_QUECTW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD8b0total, "names") <- c("CDR3", "Count")

QUECD8b150total = read.csv("dcr_beta_QUECTW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUECD8b150total, "names") <- c("CDR3", "Count")

###QUF

QUFCD8a48total = read.csv("dcr_alpha_QUFCTW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8a48total, "names") <- c("CDR3", "Count")

QUFCD8b48total = read.csv("dcr_beta_QUFCTW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8b48total, "names") <- c("CDR3", "Count")

QUFCD8a0total = read.csv("dcr_alpha_QUFCTW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8a0total, "names") <- c("CDR3", "Count")

QUFCD8a150total = read.csv("dcr_alpha_QUFCTW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8a150total, "names") <- c("CDR3", "Count")

QUFCD8b0total = read.csv("dcr_beta_QUFCTW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8b0total, "names") <- c("CDR3", "Count")

QUFCD8b150total = read.csv("dcr_beta_QUFCTW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUFCD8b150total, "names") <- c("CDR3", "Count")


###QUK

QUKCD8a48total = read.csv("dcr_alpha_QUKPTIW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8a48total, "names") <- c("CDR3", "Count")

QUKCD8b48total = read.csv("dcr_beta_QUKPTIW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8b48total, "names") <- c("CDR3", "Count")

QUKCD8a12total = read.csv("dcr_alpha_QUKPTIW12-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8a12total, "names") <- c("CDR3", "Count")

QUKCD8b12total = read.csv("dcr_beta_QUKPTIW12-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8b12total, "names") <- c("CDR3", "Count")

QUKCD8a0total = read.csv("dcr_alpha_QUKPTIW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8a0total, "names") <- c("CDR3", "Count")

QUKCD8a150total = read.csv("dcr_alpha_QUKPTIW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8a150total, "names") <- c("CDR3", "Count")

QUKCD8b0total = read.csv("dcr_beta_QUKPTIW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8b0total, "names") <- c("CDR3", "Count")

QUKCD8b150total = read.csv("dcr_beta_QUKPTIW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUKCD8b150total, "names") <- c("CDR3", "Count")

###QUL

QULCD8a48total = read.csv("dcr_alpha_QULPTIW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8a48total, "names") <- c("CDR3", "Count")

QULCD8b48total = read.csv("dcr_beta_QULPTIW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8b48total, "names") <- c("CDR3", "Count")

QULCD8a12total = read.csv("dcr_alpha_QULPTIW12-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8a12total, "names") <- c("CDR3", "Count")

QULCD8b12total = read.csv("dcr_beta_QULPTIW12-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8b12total, "names") <- c("CDR3", "Count")

QULCD8a0total = read.csv("dcr_alpha_QULPTIW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8a0total, "names") <- c("CDR3", "Count")

QULCD8a150total = read.csv("dcr_alpha_QULPTIW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8a150total, "names") <- c("CDR3", "Count")

QULCD8b0total = read.csv("dcr_beta_QULPTIW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8b0total, "names") <- c("CDR3", "Count")

QULCD8b150total = read.csv("dcr_beta_QULPTIW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QULCD8b150total, "names") <- c("CDR3", "Count")

###QUQ

QUQCD8a48total = read.csv("dcr_alpha_QUQPTIW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8a48total, "names") <- c("CDR3", "Count")

QUQCD8b48total = read.csv("dcr_beta_QUQPTIW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8b48total, "names") <- c("CDR3", "Count")

QUQCD8a12total = read.csv("dcr_alpha_QUQPTIW12-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8a12total, "names") <- c("CDR3", "Count")

QUQCD8b12total = read.csv("dcr_beta_QUQPTIW12-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8b12total, "names") <- c("CDR3", "Count")

QUQCD8a0total = read.csv("dcr_alpha_QUQPTIW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8a0total, "names") <- c("CDR3", "Count")

QUQCD8a150total = read.csv("dcr_alpha_QUQPTIW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8a150total, "names") <- c("CDR3", "Count")

QUQCD8b0total = read.csv("dcr_beta_QUQPTIW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8b0total, "names") <- c("CDR3", "Count")

QUQCD8b150total = read.csv("dcr_beta_QUQPTIW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUQCD8b150total, "names") <- c("CDR3", "Count")

###QUT

QUTCD8a48total = read.csv("dcr_alpha_QUTPTIW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8a48total, "names") <- c("CDR3", "Count")

QUTCD8b48total = read.csv("dcr_beta_QUTPTIW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8b48total, "names") <- c("CDR3", "Count")

QUTCD8a12total = read.csv("dcr_alpha_QUTPTIW12-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8a12total, "names") <- c("CDR3", "Count")

QUTCD8b12total = read.csv("dcr_beta_QUTPTIW12-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8b12total, "names") <- c("CDR3", "Count")

QUTCD8a0total = read.csv("dcr_alpha_QUTPTIW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8a0total, "names") <- c("CDR3", "Count")

QUTCD8a150total = read.csv("dcr_alpha_QUTPTIW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8a150total, "names") <- c("CDR3", "Count")

QUTCD8b0total = read.csv("dcr_beta_QUTPTIW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8b0total, "names") <- c("CDR3", "Count")

QUTCD8b150total = read.csv("dcr_beta_QUTPTIW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUTCD8b150total, "names") <- c("CDR3", "Count")

###QUV

QUVCD8a48total = read.csv("dcr_alpha_QUVCTW48-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8a48total, "names") <- c("CDR3", "Count")

QUVCD8b48total = read.csv("dcr_beta_QUVCTW48-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8b48total, "names") <- c("CDR3", "Count")

QUVCD8a0total = read.csv("dcr_alpha_QUVCTW0-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8a0total, "names") <- c("CDR3", "Count")

QUVCD8a150total = read.csv("dcr_alpha_QUVCTW150-a_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8a150total, "names") <- c("CDR3", "Count")

QUVCD8b0total = read.csv("dcr_beta_QUVCTW0-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8b0total, "names") <- c("CDR3", "Count")

QUVCD8b150total = read.csv("dcr_beta_QUVCTW150-b_CD8_cdr3.csv", header=FALSE)  # read csv file
attr(QUVCD8b150total, "names") <- c("CDR3", "Count")



################################################################
list_raw_CD4 <- list(QUDCD4a0total, QUDCD4b0total, QUDCD4a150total,
                     QUDCD4a150total)
diversity_table <- data.frame(shannon = 1:4,
                              gini = 1:4)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUECD4a0total, QUECD4b0total, QUECD4a150total,
                     QUECD4b150total)
diversity_table <- data.frame(shannon = 1:4,
                              gini = 1:4)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUECD4a0total, QUECD4b0total, QUECD4a150total,
                     QUECD4b150total)
diversity_table <- data.frame(shannon = 1:4,
                              gini = 1:4)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUFCD4a0total, QUFCD4b0total, QUFCD4a150total,
                     QUFCD4b150total,QUFCD4a48total, QUFCD4b48total)
diversity_table <- data.frame(shannon = 1:6,
                              gini = 1:6)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUKCD4a0total, QUKCD4b0total, QUKCD4a12total, QUKCD4b12total, QUKCD4a150total,
                     QUKCD4b150total,QUKCD4a48total, QUKCD4b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QULCD4a0total, QULCD4b0total, QULCD4a12total, QULCD4b12total, QULCD4a150total,
                     QULCD4b150total,QULCD4a48total, QULCD4b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUPCD4a0total, QUPCD4b0total, QUPCD4a150total,
                     QUPCD4b150total)
diversity_table <- data.frame(shannon = 1:4,
                              gini = 1:4)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUQCD4a0total, QUQCD4b0total, QUQCD4a12total, QUQCD4b12total, QUQCD4a150total,
                     QUQCD4b150total,QUQCD4a48total, QUQCD4b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUTCD4a0total, QUTCD4b0total, QUTCD4a12total, QUTCD4b12total, QUTCD4a150total,
                     QUTCD4b150total,QUTCD4a48total, QUTCD4b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD4 <- list(QUVCD4a0total, QUVCD4b0total, QUVCD4a150total,
                     QUVCD4b150total,QUVCD4a48total, QUVCD4b48total)
diversity_table <- data.frame(shannon = 1:6,
                              gini = 1:6)

diversity_table$shannon <- lapply(list_raw_CD4, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD4, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUECD8a0total, QUECD8b0total, QUECD8a150total,
                     QUECD8b150total)
diversity_table <- data.frame(shannon = 1:4,
                              gini = 1:4)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUFCD8a0total, QUFCD8b0total, QUFCD8a150total,
                     QUFCD8b150total,QUFCD8a48total, QUFCD8b48total)
diversity_table <- data.frame(shannon = 1:6,
                              gini = 1:6)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUKCD8a0total, QUKCD8b0total, QUKCD8a12total, QUKCD8b12total, QUKCD8a150total,
                     QUKCD8b150total,QUKCD8a48total, QUKCD8b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QULCD8a0total, QULCD8b0total, QULCD8a12total, QULCD8b12total, QULCD8a150total,
                     QULCD8b150total,QULCD8a48total, QULCD8b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUQCD8a0total, QUQCD8b0total, QUQCD8a12total, QUQCD8b12total, QUQCD8a150total,
                     QUQCD8b150total,QUQCD8a48total, QUQCD8b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUTCD8a0total, QUTCD8b0total, QUTCD8a12total, QUTCD8b12total, QUTCD8a150total,
                     QUTCD8b150total,QUTCD8a48total, QUTCD8b48total)
diversity_table <- data.frame(shannon = 1:8,
                              gini = 1:8)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
################################################################
list_raw_CD8 <- list(QUVCD8a0total, QUVCD8b0total, QUVCD8a150total,
                     QUVCD8b150total,QUVCD8a48total, QUVCD8b48total)
diversity_table <- data.frame(shannon = 1:6,
                              gini = 1:6)

diversity_table$shannon <- lapply(list_raw_CD8, function(x) diversity(x[,2], index = "shannon", base = 2))

diversity_table$gini <- lapply(list_raw_CD8, function(x) ineq(x[,2], type="Gini"))

diversity_table$shannon <- as.numeric(diversity_table$shannon)
diversity_table$gini <- as.numeric(diversity_table$gini)

write.csv(diversity_table, "diversity_table.csv")

############################
