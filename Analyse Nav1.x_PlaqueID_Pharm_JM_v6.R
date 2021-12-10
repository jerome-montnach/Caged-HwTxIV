## Description ----
# Objectif: Analyse d'un protocole single step dépolarisant avec une période Controle et une période avec 1 compose.
# Ne fonctionne pas pour plusieurs composés ni pour un protocole dans lequel il y a photoactivation

# Si extraction des donnees depuis l'ordinateur du patch auto, le separateur decimal est un "." et non une "," dans le point 2.2


## Functions ----

## Libraries ----
#install.packages("doParallel")
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
#install.packages("formattable")



## load librairies ----
library(doParallel)
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


# 0 - set variables ----
  # 0.1 - Dossier dans lequel se trouvent les fichiers ----
# dossier racine
maindir = "D:/AnalysePatch/Nav1.6_HwtxIVNvoc" 

# dossier du jour de manip
dir = "201201-K36/raw"

# dossier des races non corrigees 
dir.file = "20W68562_Pharm_non corr"

# dossier des races  corrigees 
dir.file.corr = "20W68562_Pharm_corr"

# 0.2 - Fichiers de parametres de patch et plaque de composes ----
# fichier avec les parametres de Seal, Capa et Rs pour chaque sweep
param_file = "20W68562_Pharm_Seal_capa_Rs.csv" 

# fichier avec le plan de plaque
plaque_file = "20W68562_Plaque.txt"  
  
  # 0.3 - Parametres d'acquisition ----  
    # Potentiels (mV)
    holding_potential = -100
    test_potential = 10
    
    # Temps debut depolarisation (us)
    Time_start_depol1 = 10000
    Time_end_depol1 = 60000
    
    # Sample periode d'acquisition (us)
    sample_period = 50
  
  # 0.4 - Curseurs et zone d'analyse ----
    # Curseurs de detection du pic (ms)
    curseur_debut_pharm_pic1  = 10.3
    curseur_fin_pharm_pic1  = 20
    curseur_debut_late_pic1  = 58
    curseur_fin_late_pic1  = 59.90
  
    # Nombre de points autour du pic pour faire la moyenne
    nb_points_pic = 1 
    
    # Decalage fit inactivation (nombre de points)
    shift_inact = 3 
    
    # 0.5 - Sweeps a utiliser pour la dose reponse ----
    # Numero du premier sweep avec compose
    sweep_comp = 53
    
    # Nombre de sweeps a moyenner pour le calcul du % d'inhibition
    n_sweeps_mean = 5
    
    # Sweep max avec compose pour le calcul du % d'inhibition
    sweep_max_comp = 200
  
  # 0.6 - QCs pour la selection des puits ----
    RSeal_min = 300 # MOhm
    courant_min = -200 # pA
    courant_max = -10000 # pA
    fuite_max = 0.4 # ratio du pic
    Rs_max = 25 # MOhm
    late_max = 50 # pA
    pourcentage_sweeps_valides = 70

## script ----
registerDoParallel(4)

# 1 - Creation du dossier d'analyse format jour et heure ----
  dir_ana = strsplit(dir,"/")[[1]][1]
  setwd(paste0(maindir,"/",dir_ana))
  date_id = format(Sys.time(), "%d-%m-%Y_%X")
  date_id = str_replace(date_id,":","-")
  date_id = str_replace(date_id,":","-")
  analyse_dir = paste0(strsplit(dir.file, "_")[[1]][2],"_","Analyse_",strsplit(dir.file, "_")[[1]][1],"_",date_id)
  dir.create(analyse_dir)      
  setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))

# 2 - Chargement des fichiers de parametres et de plan de plaque ----   
  # 2.1 - Chargement du plan de plaque ----
  setwd(paste0(maindir,"/",dir))
  plaque_id = fread(plaque_file, sep="\t",header = T, dec = ",", stringsAsFactors=FALSE)
  plaque_id = as.data.frame(plaque_id)

  # 2.2 - Chargement du fichier de parametres de patch ----
  setwd(paste0(maindir,"/",dir))
  df_param = read.csv(param_file, sep="\t",header = T, dec = ",", stringsAsFactors=FALSE, skip=1)
  colnames(df_param)[1] = "Well"

  # 2.3 - Filtre de la table de Seal_Rs_Capa ----
    # Table de Rs
    df_Rs = df_param[,grepl( "Series.Resistance", names(df_param))]
    row.names(df_Rs) = df_param[,1]
    df_Rs = df_Rs/10^6
    df_Rs[df_Rs> Rs_max] = NA
    df_Rs = df_Rs[-which(rowMeans(is.na(df_Rs)) > (100-pourcentage_sweeps_valides)/100), ]
    df_Rs = cbind(rownames(df_Rs),df_Rs)
    colnames(df_Rs)[1] = "Well"

    # Table de Seal
    df_Seal = df_param[,grepl( "Seal.Resistance", names(df_param))]
    row.names(df_Seal) = df_param[,1]
    df_Seal = df_Seal/10^6
    df_Seal[df_Seal< RSeal_min] = NA
    df_Seal = df_Seal[-which(rowMeans(is.na(df_Seal)) > (100-pourcentage_sweeps_valides)/100), ]
    df_Seal = cbind(rownames(df_Seal),df_Seal)
    colnames(df_Seal)[1] = "Well"

    # Table de Capacitance
    df_Capa = df_param[,grepl( "Capacitance", names(df_param))]
    df_Capa = cbind(df_param[,1],df_Capa)
    colnames(df_Capa)[1] = "Well"

    # Merging des tables
    df_param1 = inner_join(df_Rs,df_Seal,by="Well")
    df_param1 = inner_join(df_param1,df_Capa,by="Well")

  # 2.4 - Liste des cellules analysables ----
  wells_to_analyse = df_param1$Well
  print(paste0(length(wells_to_analyse)," puits analysables"))

# 3 - Analyse des traces de courant pour le suivi du patch auto (leak non corrige) ----
  # 3.1 - Nombre de fichiers disponibles ----
  setwd(paste0(maindir,"/",dir,"/",dir.file))
  list_files = list.files()
  print(paste0(length(list_files)," fichiers disponibles"))
  
  # 3.2 - Analyse des 10 premiers sweeps de chaque fichier ----
  options(warn=-1)
  res.df.suvi = NULL
  for (i in 1:length(list_files)){ 
    setwd(paste0(maindir,"/",dir,"/",dir.file))
    file = list_files[i]
    
    # ouverture et formatage du fichier source
    df_data0 = read.csv(file, sep=";",header = T, dec = ".", stringsAsFactors=FALSE)
    df_data0 = df_data0[,-c(1,2)]
    df_names = colnames(df_data0)
    rm(df_data0)
    df_data = read.csv(file, sep=";",header = T, dec = ",", stringsAsFactors=FALSE, skip = 2)
    df_data = df_data[,-1]
    colnames(df_data) = c("Time",df_names)
    if(is.na(df_data[1,ncol(df_data)])){df_data = df_data[,-ncol(df_data)]}
    dim(df_data)
    head(df_data[,1:10])
    
    # Extraction des parametres de patch 
    if(nrow(df_param[df_param$Well == tail(strsplit(file, "_")[[1]],n=2)[1],]) != 0){
      df_param0 = df_param[df_param$Well == tail(strsplit(file, "_")[[1]],n=2)[1],]
      # Table de Rs
      df_Rs = df_param0[,grepl( "Series.Resistance", names(df_param0))]
      # Table de Seal
      df_Seal = df_param0[,grepl( "Seal.Resistance", names(df_param0))]
      # Table de Capacitance
      df_Capa = df_param0[,grepl( "Capacitance", names(df_param0))]
    } 
    
    # Analyse de la fuite et de l'amplitude sur les 10 premiers sweeps
    res.df = NULL
    df1 = NULL
    for (sweep in 1:10){
      df = df_data[,c(1,sweep+1)]
      chip_id = paste0(head(strsplit(file, "_")[[1]],n=3)[2])
      sweep_id = as.numeric(strsplit(colnames(df)[2],"X")[[1]][2])
        
      if(!is.nan(df[1,2]) & !is.na(df[1,2])){
          df[2] = as.numeric(df[,2])
          
          # Calcul de la fuite moyenne (pA) sur les 10 premieres ms
          df_leak = df[df$Time/1000<10,]
          leak_mean = mean(df_leak[,2])/10e-13
          
          # Amplitude (pA) du pic et time (ms) to peak
          df_pic= df[df$Time/1000>curseur_debut_pharm_pic1 & df$Time/1000<curseur_fin_pharm_pic1,]
          time_pic1 = df_pic[df_pic[2]== min(df_pic[2]),1][1]
          df_data_pic = df_pic[df_pic$Time>time_pic1-nb_points_pic*sample_period & df_pic$Time<time_pic1+nb_points_pic*sample_period,]
          max_df_data1 = mean(df_data_pic[,2])/10e-13
          Tp1 = (time_pic1 - Time_start_depol1)/1000
          
          # Amplitude (pA) du courant residuel en fin de depolarisation
          df_late = df[df$Time/1000>curseur_debut_late_pic1 & df$Time/1000<curseur_fin_late_pic1,]
          mean_late1 = mean(df_late[,2])/10e-13
          
          # Formatage de la table de sortie
          res.df.suvi1 = data.frame(chip_id, tail(strsplit(file, "_")[[1]],n=2)[1],sweep_id,
                                    leak_mean,max_df_data1,df_Seal[,sweep]/10^6,df_Rs[,sweep]/10^6,df_Capa[,sweep]/10^-12)
          colnames(res.df.suvi1) = c("Chip","Well","Sweep","Leak(pA)","Max_current(pA)","Seal(MOhm)","Rs(MOhm)","Capacitance(pF)")
          res.df.suvi = rbind(res.df.suvi,res.df.suvi1)
        }
      sweep=sweep+1
    }  
    i=i+1
  }  
  options(warn=0)

  # 3.3 - Enregistrement de la table de suivi du patch-auto
  setwd(paste0(strsplit(maindir,"/")[[1]][1],"/",strsplit(maindir,"/")[[1]][2],"/SuiviPatch/",strsplit(file, "_")[[1]][1]))
  fwrite(res.df.suvi,paste0("Suivi_pharm_",strsplit(dir,"/")[[1]][1],"_",strsplit(file, "_")[[1]][2],".csv"), sep=";", dec=",")  

# 4 - Selection des cellules basee sur le ratio fuite/courant (leak non corrigee)----
  # a partir des cellules deja selectionnees a l'etape 2
  valid.well.list = NULL
  for (i in 1:length(list_files)){
    setwd(paste0(maindir,"/",dir,"/",dir.file))
    file = list_files[i]
    
    if(tail(strsplit(file, "_")[[1]],n=2)[1] %in% wells_to_analyse){
      
      # Ouverture et formatage du fichier source
      df_data0 = read.csv(file, sep=";",header = T, dec = ".", stringsAsFactors=FALSE)
      df_data0 = df_data0[,-c(1,2)]
      df_names = colnames(df_data0)
      rm(df_data0)
      df_data = read.csv(file, sep=";",header = T, dec = ",", stringsAsFactors=FALSE, skip = 2)
      df_data = df_data[,-1]
      colnames(df_data) = c("Time",df_names)
      if(is.na(df_data[1,ncol(df_data)])){df_data = df_data[,-ncol(df_data)]}
      dim(df_data)
      head(df_data[,1:10])
      
      # Analyse pour chaque sweep
      valid.well = NULL
      for (sweep in 1:(ncol(df_data)-1)){
        df = df_data[,c(1,sweep+1)]
        sweep_id = as.numeric(strsplit(colnames(df)[2],"X")[[1]][2])
        
        if(!is.nan(df[1,2]) & !is.na(df[1,2])){
          df[2] = as.numeric(df[,2])
          
          # Calcul de la fuite moyenne sur les 10 premieres ms
          df_leak = df[df$Time/1000<10,]
          leak_mean = mean(df_leak[,2])/10e-13
          
          # Amplitude du pic 
          df_pic= df[df$Time/1000>curseur_debut_pharm_pic1 & df$Time/1000<curseur_fin_pharm_pic1,]
          time_pic1 = df_pic[df_pic[2]== min(df_pic[2]),1][1]
          df_data_pic = df_pic[df_pic$Time>time_pic1-nb_points_pic*sample_period & df_pic$Time<time_pic1+nb_points_pic*sample_period,]
          max_df_data1 = mean(df_data_pic[,2])/10e-13
          
          # Amplitude max de la fuite
          if(sweep<sweep_comp){fuite_max0 = fuite_max*max_df_data1}
          if(sweep>=sweep_comp){fuite_max1 = fuite_max0}
          
          # Selection de la cellule
          if(leak_mean>fuite_max0){valid.well0 = data.frame(Well = tail(strsplit(file, "_")[[1]],n=2)[1], Sweep = sweep_id , Leak_validation = "YES")}
          if(leak_mean<=fuite_max0){valid.well0 = data.frame(Well = tail(strsplit(file, "_")[[1]],n=2)[1], Sweep = sweep_id , Leak_validation = "NO")}
          
          if(sweep<sweep_comp){
            if(max_df_data1<courant_min & max_df_data1>courant_max){valid.well0 = cbind(valid.well0,data.frame(Peak_validation = "YES"))}
            if(max_df_data1>courant_min | max_df_data1<courant_max){valid.well0 = cbind(valid.well0,data.frame(Peak_validation = "NO"))}
          }
          if(sweep>=sweep_comp){valid.well0 = cbind(valid.well0,data.frame(Peak_validation = "YES"))}
          
          valid.well = rbind(valid.well,valid.well0)
        }
        sweep=sweep+1
      }
      
      # selection des puits avec leak ok
      if(nrow(valid.well[valid.well$Leak_validation == "YES",])/(ncol(df_data)-1)>pourcentage_sweeps_valides/100){
        valid.well.list0 = valid.well[,-2]
        valid.well.list0 = valid.well.list0[valid.well.list0$Leak_validation == "YES" & valid.well.list0$Peak_validation == "YES",]
        if(nrow(valid.well.list0[valid.well.list0$Leak_validation == "YES",])/(ncol(df_data)-1)>pourcentage_sweeps_valides/100){
          valid.well.list0 = valid.well0[1,-2]
          valid.well.list = rbind(valid.well.list,valid.well.list0)
        }
      }
    }
    i=i+1
  }
  valid.well.list = valid.well.list$Well
  length(valid.well.list)

# 5 - Analyse des fichiers corriges de la fuite  ----
  setwd(paste0(maindir,"/",dir,"/",dir.file.corr))
  list_files_corr = list.files()
  
  # 5.1 - Processing des fichiers ----
  options(warn=-1)
  res.pharm = NULL
  for (i in 1:length(list_files_corr)){
    setwd(paste0(maindir,"/",dir,"/",dir.file.corr))
    file = list_files_corr[i]
    file
    
    # 5.1 - Processing de chaque fichier valide ----
    if(tail(strsplit(file, "_")[[1]],n=2)[1] %in% valid.well.list){

      df_data0 = read.csv(file, sep=";",header = T, dec = ".", stringsAsFactors=FALSE)
      df_data0 = df_data0[,-c(1,2)]
      df_names = colnames(df_data0)
      rm(df_data0)
      df_data = read.csv(file, sep=";",header = T, dec = ",", stringsAsFactors=FALSE, skip = 2)
      df_data = df_data[,-1]
      colnames(df_data) = c("Time",df_names)
      if(is.na(df_data[1,ncol(df_data)])){df_data = df_data[,-ncol(df_data)]}
      dim(df_data)
      head(df_data[,1:10])
      
      # Extraction des parametres de patch
      if(nrow(df_param[df_param1$Well == tail(strsplit(file, "_")[[1]],n=2)[1],]) != 0){
        df_param0 = df_param1[df_param1$Well == tail(strsplit(file, "_")[[1]],n=2)[1],]
        # Table de Rs
        df_Rs = t(df_param0[,grepl( "Series.Resistance", names(df_param0))])
        # Table de Seal
        df_Seal = t(df_param0[,grepl( "Seal.Resistance", names(df_param0))])
        # Table de Capacitance
        df_Capa = t(df_param0[,grepl( "Capacitance", names(df_param0))])
      } 
      
      # Analyse pour chaque sweep
      res.df = NULL
      plot_list = list()
      for (sweep in 1:(ncol(df_data)-1)){
        df = df_data[,c(1,sweep+1)]
        chip_id = paste0(head(strsplit(file, "_")[[1]],n=3)[2])
        sweep_id = as.numeric(strsplit(colnames(df)[2],"X")[[1]][2])
        
        if(!is.nan(df[1,2]) & !is.na(df[1,2])){
          df[2] = as.numeric(df[,2])
          
          # Trace du courant du sweep
          p1 = ggplot(df, aes(x=Time/1000, y=df[,2]/10e-13)) + 
            geom_line() +
            ggtitle(paste0("Courant_",tail(strsplit(file, "_")[[1]],n=2)[1],"_Sweep_",sweep_id)) + 
            scale_x_continuous(name="Time (ms)",limits = c((Time_start_depol1-3000)/1000,Time_end_depol1/1000)) +
            scale_y_continuous(name="Current (pA)", limits=c((min(df_data,na.rm = T)/10e-13), 300))
          theme(legend.position="none")
          
          # Amplitude du pic et time to peak
          df_pic= df[df$Time/1000>curseur_debut_pharm_pic1 & df$Time/1000<curseur_fin_pharm_pic1,]
          time_pic1 = df_pic[df_pic[2]== min(df_pic[2]),1][1]
          df_data_pic = df_pic[df_pic$Time>time_pic1-nb_points_pic*sample_period & df_pic$Time<time_pic1+nb_points_pic*sample_period,]
          max_df_data1 = mean(df_data_pic[,2])/10e-13
          Tp1 = (time_pic1 - Time_start_depol1)/1000
          
          # Amplitude du courant residuel (late) en fin de depolarisation
          df_late = df[df$Time/1000>curseur_debut_late_pic1 & df$Time/1000<curseur_fin_late_pic1,]
          mean_late1 = mean(df_late[,2])/10e-13
          
          # Half-time of inactivation
          half_amplitude = (max_df_data1-mean_late1)/2
          df_half = df[df$Time>=time_pic1+sample_period*shift_inact & df$Time<=Time_end_depol1-+sample_period*shift_inact,]
          half_time = (df_half[which.min(abs(half_amplitude-(df_half[,2]/10e-13))),1]-time_pic1)/1000
          
          # ratio Late/Peak
          ratio_late_peak = mean_late1/max_df_data1
          
          # Aire sous la courbe (AUC) 
          df_area = df[df$Time > Time_start_depol1 & df$Time < Time_end_depol1,]
          AUC_current = trapz(df_area$Time,df_area[,2])
          df_AUC = data.frame(Area_Under_Curve = AUC_current)
          
          # Extraction des parametres de patch
          if(length(df_Seal)!=0 & length(df_Rs)!=0){
            # Parametres du sweep
            Seal_sweep = df_Seal[sweep,1]
            if (length(df_Seal)==0){Seal_sweep = 0}
            Rs_sweep = df_Rs[sweep,1]
            if (length(df_Rs)==0){Rs_sweep = 0}
            Capa_sweep = df_Capa[sweep,1]
            if (length(df_Capa)==0){Capa_sweep = 0}
            
            # Densite de courant
            density = max_df_data1 / (Capa_sweep/10e-13)
            
            # Fit de l'inactivation du courant
            df_inact = df[df$Time>=time_pic1+sample_period*shift_inact & df$Time<=Time_end_depol1-+sample_period*shift_inact,]
            df_inact$Time = df_inact$Time/1000
            df_inact$Time = df_inact$Time-(Time_start_depol1/1000)
            colnames(df_inact) = c("Time","Current")
            df_inact$Time = as.numeric(df_inact$Time)
            
            # Fit avec 2 exponentielles
            fit_inact = tryCatch(nls(Current ~ A1*exp(-Time/Tau1) + A2*exp(-Time/Tau2) + c,
                                     start=list(Tau1=2, A1 = -5000, Tau2=10, A2 = -200, c= -1),
                                     data = df_inact, control = nls.control(maxiter = 1000, minFactor = 1e-10)),
                                 warning = function(w) {print("Warning error")},
                                 error = function(e) {print("Erreur dans les paramètres initiaux")})
            
            # Fit avec 1 exponentielle en cas d'échec avec 2 exponentielles
            if(fit_inact =="Erreur dans les paramètres initiaux"){
              fit_inact = tryCatch(nls(Current ~ A1*exp(-Time/Tau1) + c,
                                       start=list(Tau1=2, A1 = -5000, c= -1),
                                       data = df_inact, control = nls.control(maxiter = 1000, minFactor = 1e-10)),
                                   warning = function(w) {print("Warning error")},
                                   error = function(e) {print("Erreur dans les paramètres initiaux")})    
            }
            
            # Sauvegarde des données brutes du model
            #setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir,"/Fits"))
            #res.fit_inact = paste0("Fit_Inactivation_",tail(strsplit(file, "_")[[1]],n=2)[1],"_sweep_",sweep,".txt")
            #capture.output(summary(fit_inact),file = res.fit_inact)
            
            # Formatage de la table de sortie du model d'inactivation
            if(fit_inact !="Erreur dans les paramètres initiaux"){
              fit_coeff = as.data.frame(t(coefficients(fit_inact)))
              if(ncol(fit_coeff)==3){
                fit_coeff = cbind(fit_coeff[,1:2],fit_coeff)
                fit_coeff[3] = NA
                fit_coeff[4] = NA
              }
              colnames(fit_coeff)[1] = "Tau_fast"
              colnames(fit_coeff)[2] = "Af"
              fit_coeff[2] = fit_coeff[,2]/10e-13
              colnames(fit_coeff)[3] = "Tau_slow"
              colnames(fit_coeff)[4] = "As"
              fit_coeff[4] = fit_coeff[,4]/10e-13
              fit_coeff[5] = fit_coeff[,5]/10e-13
              
              # Prediction et plot de la trace de courant avec le fit
              smoothed_fit_inact = predict(fit_inact)
              df_inact = cbind(df_inact,smoothed_fit_inact)
              colnames(df_inact)[3] = "Smoothed"
              p2 = ggplot(df_inact, aes(x=Time, y=df_inact[,2]/10e-13)) + 
                geom_point(color= "black",size=2) +
                ggtitle(paste0("Fit inactivation_",tail(strsplit(file, "_")[[1]],n=2)[1],"_Sweep_",sweep_id)) + 
                geom_line(aes(Time, Smoothed/10e-13, color="red"), df_inact,size=1) +
                theme(legend.position="none")
            } else (fit_coeff = data.frame(Tau_fast = NA, Af = NA, Tau_slow = NA, As = NA, c = NA))
            
            # Formatage de la table de sortie
            res.df0 = data.frame(chip_id, tail(strsplit(file, "_")[[1]],n=2)[1],sweep_id,
                                 max_df_data1,density,
                                 Tp1,mean_late1,AUC_current, half_time,ratio_late_peak)
            colnames(res.df0) = c("Chip","Well","Sweep",
                                  "Max_current(pA)","Densite(pA/pF)",
                                  "TimeToPeak(ms)", 
                                  "Late Step1(pA)","AUC","Time half inact (ms)", "Ratio Late_Peak")
            res.df0 = cbind(res.df0,fit_coeff)
            res.df = rbind(res.df,res.df0)     
  
        } else (res.df=res.df)
          
          # Formatage des graphiques pour l'export en PDF
          if(exists("p1") & exists("p2")){
            p = plot_grid(p1, p2, labels = c('A', 'B'), ncol=1,align = "v")
            plot_list[[sweep]] = p
          }
          if(exists("p1") & !exists("p2")){
            p = plot_grid(p1, labels = c('A', 'B'), ncol=1,align = "v")
            plot_list[[sweep]] = p
          }
          if(!exists("p1") & !exists("p2")){
            plot_list[[sweep]] = p
          }
        }
        rm(p1,p2)
        sweep=sweep+1
      }
      
      # Sauvegarde des traces en PDF
      setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))
      res.plot = paste0("TracesCourants_",tail(strsplit(file, "_")[[1]],n=2)[1],".pdf")
      pdf(res.plot,width = 10, height = 10)
      for (sweep in 1:length(plot_list)) {
        print(plot_list[[sweep]])
      }
      dev.off()
      
      # Selection des cellules ayabt des sweeps en condition control et avec compose et Condition(pourcentage_sweeps_valides) valide
      res.df.ctl = res.df[res.df$Sweep < sweep_comp,]
      res.df.comp = res.df[res.df$Sweep >= sweep_comp,]
      if(!is.null(nrow(res.df.ctl)) & !is.null(nrow(res.df.comp))){
        if(nrow(res.df.ctl) != 0 & nrow(res.df.comp) != 0){
          res.df.ctl$Group = "Control"
          res.df.comp$Group = "Compose"
          res.df = rbind(res.df.ctl,res.df.comp)
          if(length(unique(res.df$Group)) == 2 & max(res.df$Sweep)/(ncol(df_data)-1) >= pourcentage_sweeps_valides/100){
            res.df = res.df
          } else (res.df = NULL)
        } else (res.df = NULL)
      } else (res.df = NULL)
      
      # Formatage de la table de sortie
      res.pharm = rbind(res.pharm,res.df)
      
      i=i+1
    }
  }
  options(warn=0)

  # 5.2 - Sauvegarde de l'environnement ----
  setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))
  save.image("Analyse_Step5.RData")
  rm(plot_list,list_files,list_files_corr, smoothed_fit_inact)

# 6 - Lien avec le plan de plaque ----
  res.df = inner_join(plaque_id,res.pharm, by="Well")
  res.df = as.data.frame(res.df)
  colnames(res.df)[3] = "Concentration(nM)"
  res.df_Condition = do.call(paste, c(res.df[c("Compose", "Concentration(nM)")], sep = "_")) 
  res.df = cbind(res.df[,1:3],res.df_Condition,res.df[,-c(1:3)])
  colnames(res.df)[4] = "Condition"
  head(res.df)
  tail(res.df)

  # Sauvegarde de la table
  setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))
  fwrite (res.df, "param_timecourse_all.csv", sep=";", dec=",")
  rm(res.pharm)
  rm(res.df_Condition, df_data, df_inact,df,df_area, df_param1)
  rm(res.df.comp,res.df.ctl, df_half,df_leak,df_pic)
  rm(valid.well,df_Capa,df_area,df_param,df_Rs,df_Seal)

# 7 - Analyse Time Course et DR ----
# reprendre au point 7.2 pour aggrager des donnees issues de diferentes plaques  
  # 7.1 - Time course de l'effet des composes ----
    # Normalisation de chaque cellule / aux sweeps avant compose
      res.df.norm = NULL
      for (well in 1:length(unique(res.df$Well))){
        res.df0 = res.df[res.df$Well == unique(res.df$Well)[well],]
        # Normalisation / aux x sweeps avant le compose
        # Moyenne des x swweps CONTROL
        res.df.ctl = res.df0[res.df0$Group == "Control",]
        res.df.ctl = res.df.ctl[res.df.ctl$Sweep > (max(res.df.ctl$Sweep)-n_sweeps_mean) & res.df.ctl$Sweep <= max(res.df.ctl$Sweep),]
        res.df.ctl = res.df.ctl[,-c(1:3,5:6,19)]
        # cas ou on a 1 expo et non pas 2 exp
        if(is.na(mean(res.df.ctl[,11])) & !is.na(mean(res.df.ctl[,9]))){
          res.df.ctl = res.df.ctl[,-c(11:12)]
          res.df.ctl = aggregate(. ~ Condition, data = res.df.ctl, FUN = mean)
          res.df.ctl = data.frame(res.df.ctl[,1:10],Tau_slow=NA,As=NA,c=res.df.ctl[,11])
        } else if(is.na(mean(res.df.ctl[,11])) & is.na(mean(res.df.ctl[,9]))){ # cas ou on a pas de fit de l'inactivation
          res.df.ctl = res.df.ctl[,-c(9:13)]
          res.df.ctl = aggregate(. ~ Condition, data = res.df.ctl, FUN = mean)
          res.df.ctl = data.frame(res.df.ctl[,1:8],Tau_fast=NA,Af=NA,Tau_slow=NA,As=NA,c=NA)
        } else (res.df.ctl = aggregate(. ~ Condition, data = res.df.ctl, FUN = mean))
        
        # Normalisation / a la moyennedes x sweeps CONTROL  
        res.df.norm1 = NULL
        for (i in 1:nrow(res.df0)){
          res.df.norm0 = res.df0[i,-c(1:3,5:6,19)]
          res.df.norm0 = res.df.norm0[,-1]/res.df.ctl[,-1]
          res.df.norm1 = rbind(res.df.norm1,res.df.norm0)
        }
        res.df.norm1 = cbind(res.df0[,1:6],res.df.norm1)  
          
        res.df.norm = rbind(res.df.norm,res.df.norm1)
        rm(res.df0,res.df.norm0,res.df.norm1)
      }
      head(res.df.norm)
      rm(res.df, plaque_id2,valid.well,well_id)
      
      # Sauvegarde de la table
      fwrite (res.df.norm, "param_timecourse_norm_all.csv", sep=";", dec=",")
    
      # save plot time course individuels
      res.plot = "TimeCourse_wells.pdf"
      pdf(res.plot, width = 20, height = 20)
      ggplot(res.df.norm, aes(x=Sweep, y=`Max_current(pA)`)) + 
        geom_point() +
        geom_vline(xintercept=sweep_comp, linetype="dashed", 
                   color = "red") +
        ggtitle("Time course Pic") + 
        facet_wrap(Condition~Well, ncol=14) +
        theme(legend.position="none")
      ggplot(res.df.norm, aes(x=Sweep, y=Tau_fast)) + 
        geom_point() +
        geom_vline(xintercept=sweep_comp, linetype="dashed", 
                   color = "red") +
        ggtitle("Time course Tau_fast") + 
        facet_wrap(Condition~Well, ncol=14) +
        theme(legend.position="none")
      ggplot(res.df.norm, aes(x=Sweep, y=Tau_slow)) + 
        geom_point() +
        geom_vline(xintercept=sweep_comp, linetype="dashed", 
                   color = "red") +
        ggtitle("Time course Tau_slow") + 
        facet_wrap(Condition~Well, ncol=14) +
        theme(legend.position="none")
      dev.off()

  
  # 7.2 - Aggregation en fonction de chaque condition (Compose + Concentration) ----
    res.df.norm1 = res.df.norm[,-c(1,5)]
    res.df.norm.courant = res.df.norm1[,1:11]
    res.df.norm.inact = res.df.norm1[,-c(5:11)]
    agg.res.df.norm.courant = do.call(data.frame, aggregate(. ~  Condition + Compose + `Concentration(nM)` + Sweep, data = res.df.norm.courant,
                                                    FUN = function(x) c(mean = mean(x), se = std.error(x), n = length(x))))
    agg.res.df.norm.inact = do.call(data.frame, aggregate(. ~  Condition + Compose + `Concentration(nM)` + Sweep, data = res.df.norm.inact,
                                                            FUN = function(x) c(mean = mean(x), se = std.error(x), n = length(x))))
    
    # Sauvegarde de la table
    fwrite (agg.res.df.norm.courant, "param_Courant_timecourse_norm_agg.csv", sep=";", dec=",")
    fwrite (agg.res.df.norm.inact, "param_Inact_timecourse_norm_agg.csv", sep=";", dec=",")
   
    rm(res.df.norm.courant,res.df.norm.inact,res.df.norm1)
    
    # save plot time course aggreges
    res.plot = "TimeCourse_agg.pdf"
    pdf(res.plot, width = 20, height = 10)
    ggplot(agg.res.df.norm.courant, aes(x=Sweep, y=Max_current.pA..mean)) + 
      geom_point() +
      geom_vline(xintercept=sweep_comp, linetype="dashed", 
                 color = "red") +
      ggtitle("Time course Pic") + 
      facet_wrap(~Condition, ncol=14) +
      theme(legend.position="none")
    ggplot(agg.res.df.norm.inact, aes(x=Sweep, y=Tau_fast.mean)) + 
      geom_point() +
      geom_vline(xintercept=sweep_comp, linetype="dashed", 
                 color = "red") +
      ggtitle("Time course Tau_fast") + 
      facet_wrap(~Condition, ncol=14) +
      theme(legend.position="none")
    ggplot(agg.res.df.norm.inact, aes(x=Sweep, y=Tau_slow.mean)) + 
      geom_point() +
      geom_vline(xintercept=sweep_comp, linetype="dashed", 
                 color = "red") +
      ggtitle("Time course Tau_slow") + 
      facet_wrap(~Condition, ncol=14) +
      theme(legend.position="none")
    dev.off()
    
    rm(agg.res.df.norm.courant,agg.res.df.norm.inact)
    rm(valid.well.list0,valid.well0,plaque_id,df_AUC,df_Capa,df_data_pic,df_late,df_Rs,df_Seal,condition_id,condition_id0,df_param,df_param0,df_plaque,res.df.ctl)
    rm(p,fit_coeff,fit_inact)
    gc()
    
  # 7.3 - Normalisation de la table normalisee / au run-up de la condition sans compose ----
    # (moyenne de toutes les cellules sans composes, si cette condition existe)
    if(nrow(res.df.norm[res.df.norm$`Concentration(nM)` == 0,]) != 0){
      res.df.norm.norm = NULL
      for (compose in 1:length(unique(res.df.norm$Compose))){
        for (sweep in 1:length(unique(res.df.norm$Sweep))){
          # moyenne des concentrations 0
          res.df.norm.ctl1 = res.df.norm[res.df.norm$Sweep == unique(res.df.norm$Sweep)[sweep] & res.df.norm$`Concentration(nM)` == 0,-c(1,2,4,5)]
          res.df.norm.ctl1 = aggregate(. ~ `Concentration(nM)`, data = res.df.norm.ctl1, FUN = mean, na.rm=TRUE, na.action=NULL)
          res.df.norm.ctl1 = res.df.norm.ctl1[,-c(1,2)]
          
          # correction sur le compose cible
          res.df.norm.ctl = res.df.norm[res.df.norm$Compose == unique(res.df.norm$Compose)[compose],]
          res.df.norm.sweep = res.df.norm.ctl[res.df.norm.ctl$Sweep == unique(res.df.norm.ctl$Sweep)[sweep],-c(1:6)]
          res.df.norm.norm1 = NULL
          for (i in 1:nrow(res.df.norm.sweep)){
            res.df.norm.norm0 = res.df.norm.sweep[i,] / res.df.norm.ctl1
            res.df.norm.norm1 = rbind(res.df.norm.norm1,res.df.norm.norm0)
          }
          res.df.norm.norm1 = cbind(res.df.norm.ctl[res.df.norm.ctl$Sweep == unique(res.df.norm.ctl$Sweep)[sweep],c(1:4,6)],res.df.norm.norm1)
          res.df.norm.norm = rbind(res.df.norm.norm,res.df.norm.norm1) 
        }
        rm(res.df.norm.ctl,res.df.norm.ctl1,res.df.norm.norm0,res.df.norm.norm1,res.df.norm.sweep)
      }
      rm(res.df.norm)
      res.df.norm.norm.courant = res.df.norm.norm[,2:12]
      res.df.norm.norm.inact = res.df.norm.norm[,-c(1,6:12)]
      agg.res.df.norm.norm.courant = do.call(data.frame, aggregate(. ~  Condition + Compose + `Concentration(nM)` + Sweep, data = res.df.norm.norm.courant,
                                                              FUN = function(x) c(mean = mean(x), se = std.error(x), n = length(x))))
      agg.res.df.norm.norm.inact = do.call(data.frame, aggregate(. ~  Condition + Compose + `Concentration(nM)` + Sweep, data = res.df.norm.norm.inact,
                                                            FUN = function(x) c(mean = mean(x), se = std.error(x), n = length(x))))
      
      # Sauvegarde de la table
      fwrite (agg.res.df.norm.norm.courant, "param_Courant_timecourse_doublenorm_agg.csv", sep=";", dec=",")
      fwrite (agg.res.df.norm.norm.inact, "param_Inact_timecourse_doublenorm_agg.csv", sep=";", dec=",")
    } else (res.df.norm.norm = res.df.norm)
    rm(res.df.norm.ctl,res.df.norm.ctl1,res.df.norm.norm.inact,res.df.norm.norm.courant,res.df.norm.norm0,res.df.norm.norm1,res.df.norm.sweep,res.df.norm)
    rm(agg.res.df.norm.norm.courant,agg.res.df.norm.norm.inact)
  
  # 7.4 - Formatage pour ratio Compose / Control ----
  # concentration réponse à partir de res.df.norm.norm  
    res.df.DR = res.df.norm.norm[,-c(4)]
    for (i in 1:nrow(res.df.DR)){
      if(res.df.DR[i,4]< sweep_comp){res.df.DR[i,17] = "Control"}
      if(res.df.DR[i,4]> sweep_comp){res.df.DR[i,17] = "Compose"}
    }
    colnames(res.df.DR)[17] = "Group"
    res.df.DR_control = res.df.DR[res.df.DR$Group == "Control",]
    res.df.DR_control = res.df.DR_control[!is.na(res.df.DR_control$Sweep),]
    res.df.DR_control = res.df.DR_control[res.df.DR_control$Sweep > (max(res.df.DR_control$Sweep)-n_sweeps_mean),]
    res.df.DR_control.courant = res.df.DR_control[,c(1:11,17)]
    res.df.DR_control.inact = res.df.DR_control[,-c(5:11)]
    agg.res.df.DR_control.courant = aggregate(. ~ Well + Compose + `Concentration(nM)` + Group, data = res.df.DR_control.courant, FUN = mean)
    agg.res.df2_control.inact = aggregate(. ~ Well + Compose + `Concentration(nM)` + Group, data = res.df.DR_control.inact, FUN = mean)
  
    res.df.DR_compose = res.df.DR[res.df.DR$Group == "Compose",]
    res.df.DR_compose = res.df.DR_compose[res.df.DR_compose$Sweep > (sweep_max_comp-n_sweeps_mean) & res.df.DR_compose$Sweep <= sweep_max_comp,]
    res.df.DR_compose.courant = res.df.DR_compose[,c(1:11,17)]
    res.df.DR_compose.inact = res.df.DR_compose[,-c(5:11)]
    agg.res.df.DR_compose.courant = aggregate(. ~ Well + Compose + `Concentration(nM)` + Group, data = res.df.DR_compose.courant, FUN = mean)
    agg.res.df.DR_compose.inact = aggregate(. ~ Well + Compose + `Concentration(nM)` + Group, data = res.df.DR_compose.inact, FUN = mean)
    
    res.df.DR.courant = rbind(agg.res.df.DR_control.courant,agg.res.df.DR_compose.courant)
    res.df.DR.inact = rbind(agg.res.df2_control.inact,agg.res.df.DR_compose.inact)
    
    # Sauvegarde des tables
    fwrite(res.df.DR.courant,"param_Courant_DR_agg.csv", sep=";", dec=",")
    fwrite(res.df.DR.inact,"param_Inact_DR_agg.csv", sep=";", dec=",")
  
  # 7.5 - Ratio Compose / Control ----
    # ratio courant  
    res.df.DR.courant.ratio = NULL
    for (well in unique(res.df.DR.courant$Well)){
      res.df.DR.courant.ratio0 = res.df.DR.courant[res.df.DR.courant$Well == well,]
      id = res.df.DR.courant.ratio0[1,1:3]
      res.df.DR.courant.ratio0 = res.df.DR.courant.ratio0[,-c(1:5)]
      res.df.DR.courant.ratio0 = res.df.DR.courant.ratio0[2,]/ res.df.DR.courant.ratio0[1,]
      colnames(res.df.DR.courant.ratio0) = c("Ratio_Peak_Corr","Ratio_Density","Ratio_TimeToPeak", "Ratio_Late_Corr","Ratio_AUC","Ratio_T0.5inact","Ratio_Late_Peak")
      res.df.DR.courant.ratio0 = cbind(id,res.df.DR.courant.ratio0)
      res.df.DR.courant.ratio = rbind(res.df.DR.courant.ratio,res.df.DR.courant.ratio0)
    }
    res.df.DR.courant.ratio[3] = as.numeric(res.df.DR.courant.ratio[,3])
    
    # ratio fit inactivation
    res.df.DR.inact.ratio = NULL
    for (well in unique(res.df.DR.inact$Well)){
      res.df.DR.inact.ratio0 = res.df.DR.inact[res.df.DR.inact$Well == well,]
      id = res.df.DR.inact.ratio0[1,1:3]
      res.df.DR.inact.ratio0 = res.df.DR.inact.ratio0[,-c(1:5)]
      res.df.DR.inact.ratio0 = res.df.DR.inact.ratio0[2,]/ res.df.DR.inact.ratio0[1,]
      colnames(res.df.DR.inact.ratio0) = c("Ratio_Tau_fast","Ratio_Af","Ratio_Tau_slow","Ratio_As","Ratio_c")
      res.df.DR.inact.ratio0 = cbind(id,res.df.DR.inact.ratio0)
      res.df.DR.inact.ratio = rbind(res.df.DR.inact.ratio,res.df.DR.inact.ratio0)
    }
    res.df.DR.inact.ratio[3] = as.numeric(res.df.DR.inact.ratio[,3])
    
    # Sauvegarde des tables de ratio
    fwrite(res.df.DR.courant.ratio,"ratio_Courant_DR.csv", sep=";", dec=",")
    fwrite(res.df.DR.inact.ratio,"ratio_Inact_DR.csv", sep=";", dec=",")
  
  rm(list = ls())
  gc()
  
  