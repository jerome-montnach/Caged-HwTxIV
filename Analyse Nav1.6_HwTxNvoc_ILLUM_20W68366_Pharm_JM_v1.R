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
    maindir = "Z:/AnalysePatch/Nav1.6_HwTxIV-Nvoc" 
    
    # dossier du jour de manip
    dir = "200930-Nav1.6-ST/raw"
    
    # dossier des races non corrigees 
    dir.file = "20W68366_Pharm_non corr"
    
    # dossier des races  corrigees 
    dir.file.corr = "20W68366_Pharm_corr"
  
  # 0.2 - Fichiers de parametres de patch et plaque de composes ----
    # fichier avec les parametres de Seal, Capa et Rs pour chaque sweep
    param_file = "20W68366_Pharm_Seal_capa_Rs.csv" 
    
    # fichier avec le plan de plaque
    plaque_file = "20W68366_Plaque.txt" 
  
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
    sweep_comp = 42
    
    # Numero du premier sweep avec illumination
    sweep_light = 75
    
    # Nombre de sweeps a moyenner pour le calcul du % d'inhibition
    n_sweeps_mean = 5
    
    # Sweep max avec compose pour le calcul du % d'inhibition
    sweep_max_comp = 260
  
  # 0.6 - QCs pour la selection des puits ----
    RSeal_min = 100 # MOhm
    courant_min = -200 # pA
    courant_max = -10000 # pA
    fuite_max = 0.4 # ratio du pic
    Rs_max = 25 # MOhm
    late_max = 50 # pA
    pourcentage_sweeps_valides = 70

## script ----
registerDoParallel(4)
memory.limit(18000)

# 1 - Creation du dossier d'analyse format jour et heure ----
  dir_ana = strsplit(dir,"/")[[1]][1]
  setwd(paste0(maindir,"/",dir_ana))
  date_id = format(Sys.time(), "%d-%m-%Y_%X")
  date_id = str_replace(date_id,":","-")
  date_id = str_replace(date_id,":","-")
  analyse_dir = paste0(strsplit(dir.file, "_")[[1]][2],"_","Analyse_",strsplit(dir.file, "_")[[1]][1],"_",date_id)
  dir.create(analyse_dir)      
  setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))
  dir.create("Fits")

# 2 - Chargement des fichiers de parametres et de plan de plaque ----   
  # 2.1 - Chargement du plan de plaque ----
  setwd(paste0(maindir,"/",dir))
  plaque_id = fread(plaque_file, sep="\t",header = T, dec = ",", stringsAsFactors=FALSE)
  plaque_id = as.data.frame(plaque_id)

  # 2.2 - Chargement du fichier de parametres de patch ----
  setwd(paste0(maindir,"/",dir))
  df_param = read.csv(param_file, sep="\t",header = T, dec = ",", stringsAsFactors=FALSE, skip=1)
  colnames(df_param)[1] = "Well"
  df_param = df_param[df_param$Well != "",]

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
  

# 4 - Selection des cellules basee sur le ratio fuite/courant (leak non corrigee)----
  # a partir des cellules deja selectionnees a l'etape 2
  valid.well.list = NULL
  for (i in 1:length(list_files)){
    setwd(paste0(maindir,"/",dir,"/",dir.file))
    file = list_files[i]
    
    if(tail(strsplit(file, "_")[[1]],n=2)[1] %in% wells_to_analyse){
      
      # Ouverture et formatage du fichier source
      df_data0 = read.csv(file, sep=";",header = T, dec = ",", stringsAsFactors=FALSE)
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
    
    # 5.1 - Processing de chaque fichier valide ----
    if(tail(strsplit(file, "_")[[1]],n=2)[1] %in% valid.well.list){

      df_data0 = read.csv(file, sep=";",header = T, dec = ",", stringsAsFactors=FALSE)
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
            setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir,"/Fits"))
            res.fit_inact = paste0("Fit_Inactivation_",tail(strsplit(file, "_")[[1]],n=2)[1],"_sweep_",sweep,".txt")
            capture.output(summary(fit_inact),file = res.fit_inact)
            
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
                                 Tp1,mean_late1,AUC_current)
            colnames(res.df0) = c("Chip","Well","Sweep",
                                  "Max_current(pA)","Densite(pA/pF)",
                                  "TimeToPeak(ms)", 
                                  "Late Step1(pA)","AUC")
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
  rm(plot_list)
  save.image("Analyse_Step5.RData")

# 6 - Lien avec le plan de plaque ----
  res.df = inner_join(plaque_id,res.pharm, by="Well")
  res.df = as.data.frame(res.df)
  colnames(res.df)[3] = "Concentration(nM)"
  res.df_Condition = do.call(paste, c(res.df[c("Compose", "Concentration(nM)")], sep = "_")) 
  res.df = cbind(res.df[,1:3],res.df_Condition,res.df[,-c(1:3)])
  colnames(res.df)[4] = "Condition"
  head(res.df)
  tail(res.df)
  unique(res.df$Well)
  
  res.df = res.df[res.df$Well !="J18",]


  

  # Sauvegarde de la table
  setwd(paste0(maindir,"/",dir_ana,"/",analyse_dir))
  fwrite (res.df, "param_timecourse_all.csv")
  
  rm(list = ls())
  