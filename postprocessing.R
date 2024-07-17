library("readxl")
library("writexl")
library("dplyr")
library("rxode2")
library("ggplot2")
library("bayesplot")
library("deSolve")
library("nimble")
library("reshape2")

#install.packages("installr")
#library(installr)
#install.rtools()
#pkgbuild::has_build_tools(debug = TRUE)
getwd()
#setwd("C:/Users/Юлия/Desktop/Summer Internship 2024 M&S Decisions")

df_for_plots<- function(ds_path, db_path){
  data_eff1 <- read_excel(ds_path)
  data_eff1<- data_eff1[2:31]
  #delete 102
  data_eff1<-data_eff1[data_eff1$ID!='102_1',]
  ##for labels with dose add Dose column
  data_dose <- read_excel(db_path,  sheet = "Arms")
  #View(data_dose)
  arms<- unique(data_eff1$ID)
  vect<- c()
  for (arm in arms){
    n_items <-nrow(filter(data_eff1, ID==arm))
    vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Dose,n_items)))}
  data_eff1$Dose<-vect
  
  col_n<-c("TIME", "DV", "CMT", "ADM", "AMT", "EVID","MDV","Sex_M", "Sex_F", "BMI", "BMI_SD", "Weight", "Weight_SD", "Age", "Age_SD", 
           "Duration_of_T1D", "Duration_of_T1D_SD", "Meal_kcal", "Meal_carbohydrates","Meal_protein", "Meal_fat","Nsub", "Dose") 
  
  for (col in col_n){
    
    data_eff1[col]<-as.numeric(unlist(data_eff1[col]))} 
  
  plot_data2 <- data_eff1[data_eff1$DVNAME=='LISPRO_PK',]
  
  #data with dose
  plot_data2 %>%
    select(TIME, DV, Dose, ID) %>%
    mutate(Arm_dose= paste(paste(ID, as.character(Dose)), 'U')) -> plot_data2_dose
  return(plot_data2_dose)
}

df_for_param_plots<- function(file_paths, num_chains, num_params){
  if (num_chains==1){
    data_out1<- read.csv(file_paths[1])
    data_out_tar1<-data_out1[,8:(7+num_params)]
    data_out_tar_time1 <- stack(data_out_tar1)
    names(data_out_tar_time1)<- c("values", "param")
    data_out_tar_time1$Chain<- rep("Chain 1", 1000)
    data_out_tar_merged<-data_out_tar_time1
    data_list<-list(data_out1, data_out_tar1, data_out_tar_merged)
    return(data_list)
  }
  else{
    data_out1<- read.csv(file_paths[1])
    data_out2<- read.csv(file_paths[2])
    data_out3<- read.csv(file_paths[3])
    data_out4<- read.csv(file_paths[4])
    #datas<-c(data_out1, data_out2, data_out3, data_out4)
    data_out_tar1<-data_out1[,8:(7+num_params)]
    data_out_tar2<-data_out2[,8:(7+num_params)]
    data_out_tar3<-data_out3[,8:(7+num_params)]
    data_out_tar4<-data_out4[,8:(7+num_params)]
    
    data_out_tar_time1 <- stack(data_out_tar1)
    data_out_tar_time2 <- stack(data_out_tar2)
    data_out_tar_time3 <- stack(data_out_tar3)
    data_out_tar_time4 <- stack(data_out_tar4)
    
    names(data_out_tar_time1)<- c("values", "param")
    names(data_out_tar_time2)<- c("values", "param")
    names(data_out_tar_time3)<- c("values", "param")
    names(data_out_tar_time4)<- c("values", "param")
    
    data_out_tar_time1$Chain<- rep("Chain 1", 1000)
    data_out_tar_time2$Chain<- rep("Chain 2", 1000)
    data_out_tar_time3$Chain<- rep("Chain 3", 1000)
    data_out_tar_time4$Chain<- rep("Chain 4", 1000)
    
    data_out_tar_merged<-rbind(data_out_tar_time1, data_out_tar_time2, data_out_tar_time3, data_out_tar_time4)
    data_list<-list(data_out1, data_out_tar1, data_out_tar_merged)
    return(data_list)
  }
}

lispro_1cmp <-rxode({
  cmt(Ad)
  cmt(Ac)
  adlag <- 1
  frac1 <- 1
  alag(Ad) <- adlag
  #f(Ad) <- frac1
  d/dt(Ad) <- -ka * Ad
  d/dt(Ac) <- ka * Ad - Cl * Ac / Vd
  Cc <- Ac / Vd;
})

df_plot_diff_eq <- function(doses_df, tim, Frac0,tinf,theta, func=lispro_1cmp){
  Dose <-doses_df[1,1]
  ev <- et(timeUnits="hr") %>%
    et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
    et(amt=Dose*(1 - Frac0), cmt = 1) %>%
    et(tim)
  Cc<-as.numeric(unlist(rxSolve(lispro_1cmp, theta, ev)['Cc']))
  Arm_dose<-rep(doses_df[1,2],100) 
  time<-tim
  df<-data.frame(Cc,Arm_dose,time)
  
  for (i in 2:8) {
    Dose <-doses_df[i,1]
    ev <- et(timeUnits="hr") %>%
      #et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
      #et(amt=Dose*(1 - Frac0), cmt = 1) %>%
      et(amt=Dose*(1-Frac0), dur = tinf, cmt = 2) %>%
      et(amt=Dose*Frac0, cmt = 1) %>%
      
      et(tim)
    Cc<-as.numeric(unlist(rxSolve(func, theta, ev)['Cc']))
    Arm_dose<-rep(doses_df[i,2],100)
    time<-tim
    df_new<-data.frame(Cc,Arm_dose, time)
    df<- rbind(df,df_new)
  }
  return(df)
}


#1. param distr plot
#2. plot of params in time
#The same but also for 4 chains
#3. plot mean, median, quantile solutions
#4. select some parameters from mean
#5. solve ode with params and compare with data  


##Analyze results
ds_path<-"dataset_lispro_with86_Dosepmol_v16.xlsx"
db_path<-"DB_Insulin_Lispro_longBL.xlsx"
plot_data2_dose<-df_for_plots(ds_path,db_path)

num_chains <- 1
num_params <- 7
one_stud<-FALSE
file_paths<-c("output_v11_mod.csv")
#file_paths<-c("output_1_v8_mod.csv","output_2_v8_mod.csv","output_3_v8_mod.csv","output_4_v8_mod.csv")
dfs<-df_for_param_plots(file_paths, num_chains, num_params)
data_out1<-dfs[1][[1]]
data_out_tar1<-dfs[2][[1]]
data_plot_param<-dfs[3][[1]]#df_for_param_plots(file_paths, num_chains)

#params in time
ggplot()+geom_line(data = data_plot_param,
                   mapping = aes(x=rep(seq(1, 1000),num_params*num_chains), y = values,
                                 color=Chain)) +facet_wrap(~param, scales="free")+
  labs(x = "iter", y = "values")

#ggsave("param_change_v9_merged.png",dpi=350,height=5,width=7.5)

#plot param distribution 
ggplot(data = data_plot_param,
       mapping = aes(x= values,
                     color = Chain)) + geom_density()+facet_wrap(~param, scales="free")

#ggsave("param_distr_v9_merged2.png",dpi=350,height=5,width=7.5)





#x_stan<-data_out1[,47:84] #93
#x_stan<-data_out1[,48:89] #103
#x2_stan <- x_stan[,seq(2, 38, by=2)] #93
#x2_stan <- x_stan[,seq(2, 42, by=2)] #103
#x_stan<-data_out[,189:510] #all
x_stan<-data_out1[,175:468] #all without 102 with frac
#x_stan<-data_out1[,174:467] #all without 102
#cHat_stan<-data_out1[,468:614] #all without 102
cHat_stan<-data_out1[,469:615] #all without 102 with frac
#cHatObs_stan<-data_out1[,615:753] # all without 102
#View(x_stan)
#x2_stan <- x_stan[,seq(2, 322, by=2)] #all
x2_stan <- x_stan[,seq(2, 294, by=2)] #all without 102
x2_stan_normV <- x2_stan / data_out_tar1$V_lis_hat



mean_param<- apply(data_out_tar1,2,mean)
median_param<-apply(data_out_tar1,2,median)
print(mean_param)
print(median_param)

mean_sol<- apply(x2_stan_normV,2,mean)
median_sol<-apply(x2_stan_normV,2,median)

t_mean <- plot_data2_dose$TIME
Arm_dose<-plot_data2_dose$Arm_dose
plot_mean <- data.frame(t_mean,mean_sol, Arm_dose)
names(plot_mean)<-c('t','sol','Arm_dose')
plot_median <- data.frame(t_mean, median_sol, Arm_dose)
names(plot_median)<-c('t','sol','Arm_dose')
iter_ind<-100
iter_sol<-as.numeric(x2_stan_normV[iter_ind,])
plot_iter <- data.frame(t_mean,iter_sol, Arm_dose)
names(plot_iter)<-c('t','sol','Arm_dose')
#View(plot_mean)
doses <- unique((plot_data2_dose[plot_data2_dose$TIME==0,]$Dose)*3.47*10^7/5813.68)
Arm_dose<- unique(plot_data2_dose$Arm_dose)
doses_df<- data.frame(doses, Arm_dose)

tim<-seq(0, 8, length.out=100)

doses_df$doses<-as.numeric(doses_df$doses)

#define sets of params for diff eq solution
mean_p <- as.numeric(mean_param)
theta_mean <- c('ka' = mean_p[1], 'Cl' = mean_p[2], 'Vd' = mean_p[3], 'adlag' = mean_p[4])
Frac0_mean <- mean_p[6]#0.271#mean_p[6]#1-0.729=0.271
tinf_mean <- mean_p[5]

median_p <- as.numeric(median_param)
theta_med <- c('ka' = median_p[1], 'Cl' = median_p[2], 'Vd' = median_p[3], 'adlag' = median_p[4])
Frac0_med <- median_p[6]#0.271#mean_p[6]#1-0.729=0.271
tinf_med <- median_p[5]

theta_hand <- c('ka' = mean_p[1], 'Cl' = 96, 'Vd' = 70, 'adlag' = mean_p[4])

df<-df_plot_diff_eq(doses_df, tim, Frac0_mean,tinf_mean,theta_mean)
df_med<-df_plot_diff_eq(doses_df, tim, Frac0_med,tinf_med,theta_med)
df_hand<-df_plot_diff_eq(doses_df, tim, Frac0_mean,tinf_mean,theta_hand)
#View(df)

Arm_dose_one<-'103_1 10 U'
Id_one<-'103_1'
if (one_stud==TRUE){
df<-df[df$Arm_dose=='103_1 10 U',]
plot_data2_dose<-plot_data2_dose[plot_data2_dose$ID=='103_1',]
}
#View(df)

#plot 
ggplot()+geom_point(data = plot_data2_dose,
                    mapping = aes(x = TIME, 
                                  y = DV, color="data")) +
  facet_wrap(~Arm_dose)+
  geom_point(data = plot_mean, mapping = aes(x = t, y = sol, color="Torsten"))+
#  facet_wrap(~Arm_dose)+ 
  geom_line(data = df, mapping = aes(x = time, y = Cc, color="diff eq mean"))+
  facet_wrap(~Arm_dose)+ 
  geom_line(data = df_med, mapping = aes(x = time, y = Cc, color="diff eq med"))+
  #facet_wrap(~Arm_dose)+ 
labs(x = " time, h",
     y = "Serum Insulin, pmol/L")

#ggsave("datapk_torsten_v11.png",dpi=350,height=5,width=7.5)

