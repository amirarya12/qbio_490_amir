
# sets working directory
setwd("/Users/amirarya/Desktop/qbio_490_amir/analysis_data")
# initializes all the necessary packages required view the data set
library(BiocManager)
library(TCGAbiolinks)
library(survival)
library(survminer)


clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")

clinical <- GDCprepare_clinic(clinical_query, clinical.info = "patient")

clinical <- read.csv("/Users/amirarya/Desktop/qbio_490_amir/brca_clinical_data.csv")

clinical.drug <- GDCprepare_clinic(query = clinical_query,
                                   clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query,
                                  clinical.info = "radiation")
# adds in the clinical data set in to the script 

is.na(clinical$lymph_node_examined_count)
sum(is.na(clinical$lymph_node_examined_count))
# sum of n/a values = 139

sum(!is.na(clinical$lymph_node_examined_count))
# sum of non n/a values = 1035

sum(is.na(clinical$race_list))
# Number of N/A: 0

sum(!is.na(clinical$race_list))
# Number of valid values: 1174



filter_mask1 = ifelse(is.numeric(clinical$lymph_node_examined_count), TRUE, FALSE)
clinical_filtered1 = clinical[filter_mask1, ]
filter_mask2 = ifelse(!is.na(clinical_filtered1$race_list), TRUE, FALSE)
clinical_filtered2 = clinical_filtered1[filter_mask2, ]
filter_mask3 = ifelse(clinical_filtered2$race_list != "", TRUE, FALSE)
clinical_filtered3 = clinical_filtered2[filter_mask3, ]
# creates masks to filter out na values, and blank spaces that are in the lymph node and race colums.

plot(factor(clinical_filtered3$race_list), clinical_filtered3$lymph_node_examined_count, xlab = "Race", ylab = "Lymph Node Examined Count")
# plots the lymph node count data against race list and labels the x axis as race and y axis as lymph node examined count



#Survival plots
clinical_filtered3$survival_time = ifelse(is.na(clinical_filtered3$days_to_death), clinical_filtered3$days_to_last_followup, clinical_filtered3$days_to_death)
clinical_filtered3$death_event = ifelse(clinical_filtered3$vital_status == "Dead", TRUE, FALSE)


surv_object_age <- Surv(time = clinical_filtered3$survival_time,
                        event = clinical_filtered3$death_event)


race_fit <- surv_fit(surv_object_age ~ clinical_filtered3$race_list,
                     data = clinical_filtered3 )

survplot_race = ggsurvplot(race_fit,
                           pval=TRUE,
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                           legend = "right")


KM_plot_race = survplot_race$plot +
  theme_bw() +  #
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

KM_plot_race
# this whole section creates a plot of race against survival rates 


#Stratification of lymph node variable
clinical_filtered3$lymphgroups = ifelse(clinical_filtered3$lymph_node_examined_count < 9, "0-8 Lymph Nodes Examined",
                                        ifelse(clinical_filtered3$lymph_node_examined_count < 18, "9-17 Lymph Nodes Examined",
                                               ifelse(clinical_filtered3$lymph_node_examined_count < 27, "18-26 Lymph Nodes Examined",
                                                      ifelse(clinical_filtered3$lymph_node_examined_count < 36, "27-35 Lymph Nodes Examined", "36-44 Lymph Nodes Examined"
                                                      ))))

#Survival plot for lymph node
clinical_filtered3$survival_time = ifelse(is.na(clinical_filtered3$days_to_death), clinical_filtered3$days_to_last_followup, clinical_filtered3$days_to_death)
clinical_filtered3$death_event = ifelse(clinical_filtered3$vital_status == "Dead", TRUE, FALSE)

surv_object_age <- Surv(time = clinical_filtered3$survival_time,
                        event = clinical_filtered3$death_event)

lymph_fit <- surv_fit(surv_object_age ~ clinical_filtered3$lymphgroups,
                      data = clinical_filtered3 )

survplot_lymph = ggsurvplot(lymph_fit,
                            pval=TRUE,
                            ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")),
                            legend = "right")

KM_plot_lymph = survplot_lymph$plot +
  theme_bw() +  #
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

KM_plot_lymph
# this sections creates a survival plot for amount of lymph nodes and survival rate

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager :: install("maftools")

