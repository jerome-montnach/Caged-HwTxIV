## Description

## Functions

## Libraries ----
#install.packages("parallel")
#install.packages("tidyverse", dependencies = T)
#install.packages("data.table", dependencies = T)
#install.packages("FedData", dependencies = T)
#install.packages("ggplot2", dependencies = T)
#install.packages("cowplot",dependencies = T)
#install.packages("qqpubr", dependencies = T)
#install.packages("settings", dependencies = T)
#install.packages("nls2", dependencies = T)
#install.packages("backports")
#install.packages("xlsx")


## load librairies ----
library("settings")
library("data.table")
library("FedData")
library("dplyr")
library("ggplot2")
library("cowplot")
library("tidyverse")
library("spectral")
library("FinCal")
library("stringr")
library("pracma")
library("pspline")
library("nls2")
library("xlsx")

## set variables ----
# Dossier dans lequel se trouvent les fichiers
# dossier racine
maindir = "D:/AnalysePatch/Nav1.6_HwtxIVNvoc" 

# dossier du jour de manip
dir = "201201-K36/raw"

# fichier avec le plan de plaque
plaque_file = "20W68562_Plaque" 

## Script ----
# Plan de plaque ----
setwd(paste0(maindir,"/",dir))


df_plaque = read.csv(paste0(plaque_file,".csv"), sep=";",header = F, dec = ",", stringsAsFactors=FALSE)
row.names(df_plaque) = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
colnames(df_plaque) = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24")

well_id = NULL
cellules_id = NULL
condition_id = NULL
for (i in 1:nrow(df_plaque)){
  well = paste0(row.names(df_plaque)[i],colnames(df_plaque))
  well_id =  c(well_id,well)
  condition_id0 = t(df_plaque[i,])
  colnames(condition_id0) = "Compose"
  condition_id = rbind(condition_id,condition_id0)
}
well_id = as.data.frame(well_id)
colnames(well_id) = "Well"
plaque_id = cbind(well_id,condition_id)
plaque_id [1] = as.character(plaque_id[,1])
plaque_id [2] = as.character(plaque_id[,2])
for (i in 1:nrow(plaque_id)){
  plaque_id[i,3] = strsplit(plaque_id[i,2],"_")[[1]][1]
  plaque_id[i,4] = strsplit(plaque_id[i,2],"_")[[1]][2]
}
plaque_id = plaque_id[,-2]
colnames(plaque_id)[2] = "Compose"
colnames(plaque_id)[3] = "Concentration_nM"
plaque_id [3] = str_replace(plaque_id [,3],",",".")
plaque_id [3] = as.numeric(plaque_id [,3])

fwrite(plaque_id,paste0(plaque_file,".txt"),sep="\t")

#plan plaque patch auto
plaque = rep(1,nrow(plaque_id))
plaque_id2 = cbind(plaque,plaque_id)
colnames(plaque_id2) = c("PLATE","POSITION","COMPOUND","CONCENTRATION")
fwrite(plaque_id2,paste0(plaque_file,"_patchauto.txt"),sep="\t")
