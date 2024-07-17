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

getwd()
#setwd("C:/Users/Юлия/Desktop/Summer Internship 2024 M&S Decisions")


data_eff1 <- read_excel("dataset_lispro_with86_Dosepmol_v16.xlsx")
data_eff1<- data_eff1[2:31]

data_eff1<-data_eff1[data_eff1$ID!='102_1',]

##for labels with dose add Dose column

data_dose <- read_excel("DB_Insulin_Lispro_longBL.xlsx",  sheet = "Arms")
View(data_dose)


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


#View(plot_data2_dose)

#1. param distr plot
#2. plot of params in time
#The same but also for 4 chains
#3. plot mean, median, quantile solutions
#4. select some parameters from mean
#5. solve ode with params and compare with data  


##Analyze results

num_chains <-4

data_out1<- read.csv("output_1_v9_mod.csv")
data_out2<- read.csv("output_2_v8_mod.csv")
data_out3<- read.csv("output_3_v8_mod.csv")
data_out4<- read.csv("output_4_v8_mod.csv")
View(data_out1)

#x_stan<-data_out[,189:510]
x_stan<-data_out1[,174:467]
cHat_stan<-data_out1[,468:614]
#cHatObs_stan<-data_out1[,615:753]
#View(cHat_stan)
#x2_stan <- x_stan[,seq(2, 322, by=2)]
x2_stan <- x_stan[,seq(2, 294, by=2)]

data_out_tar1<-data_out1[,8:13]
data_out_tar2<-data_out2[,8:13]
data_out_tar3<-data_out3[,8:13]
data_out_tar4<-data_out4[,8:13]

mean_param<- apply(data_out_tar1,2,mean)

median_param<-apply(data_out_tar1,2,median)


print(mean_param)
print(median_param)
#plot parameters changing in time
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

View(data_out_tar_merged)
#params in time
num_chains <-4
ggplot()+geom_line(data = data_out_tar_merged,
                    mapping = aes(x=rep(seq(1, 1000),6*num_chains), y = values,
                                  color = Chain)) +facet_wrap(~param, scales="free")+
labs(x = "iter", y = "values")

ggsave("param_change_v9_merged.png",dpi=350,height=5,width=7.5)

#plot param distribution 

ggplot(data = data_out_tar_merged,
                   mapping = aes(x= values,
                                 color = Chain)) + geom_density()+facet_wrap(~param, scales="free")

ggsave("param_distr_v9_merged2.png",dpi=350,height=5,width=7.5)

View(data)
#individual plots:
num_chains<- 1
ggplot()+geom_line(data = data_out_tar_time1,
                   mapping = aes(x=rep(seq(1, 1000),6*num_chains), y = values,
                                 color = param)) +facet_wrap(~param, scales="free")+
  labs(x = "iter", y = "values")

ggsave("param_change_103_frac02.png",dpi=350,height=5,width=7.5)

ggplot(data = data_out_tar_time1,
       mapping = aes(x= values,
                     color = param)) + geom_density()+facet_wrap(~param, scales="free")

ggsave("param_distr_103_frac02.png",dpi=350,height=5,width=7.5)





#plot mean solution from Torsten


x2_stan_normV <- x2_stan / data_out_tar1$V_lis_hat

View(x2_stan_normV)

mean_sol<- apply(x2_stan_normV,2,mean)

median_sol<-apply(x2_stan_normV,2,median)


t_mean <- plot_data2_dose$TIME
Arm_dose<-plot_data2_dose$Arm_dose
plot_mean <- data.frame(t_mean,mean_sol, Arm_dose)

View(plot_data2_dose)
#solve diff eq with new params
#ode_rhs <-function(t,y,p=param_arr) {

#  Ka <-p[1]
#  CL<-p[2]
#  V<-p[3]
#  tlag<-p[4]
#  tinf<-p[5]
#  frac<-p[6]
#  dose<-p[7]
#  CL_V <- CL / V
  
#  k0 <- (1-frac)*dose/tinf
 
#  if (t>= tinf) {
#    k0 <- 0} 
  
#  k0_cond - CL_V*y}
#install.packages("pracma")
#library("pracma")
#sol <- pracma::ode45(ode_rhs, t0, tend, y0)

ggplot()+geom_point(data = plot_data2_dose,
       mapping = aes(x = TIME, 
                     y = DV,
                     color = Arm_dose)) +
  facet_wrap(~Arm_dose)+
  #geom_point(alpha = .7, size = 2) +
  geom_point(data = plot_mean, mapping = aes(x = t_mean, y = mean_sol))+
  facet_wrap(~Arm_dose) #+
  #labs(title = "Mean insulin concentration versus time after single s.c. dose of insulin lispro",
  #     x = " time, h",
  #     y = "Serum Insulin, pmol/L")
ggsave("pkdata_torstensolch1_v9.png",dpi=350,height=5,width=7.5)

#pkgbuild::has_build_tools(debug = TRUE)


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

#View(plot_data2_dose)

#select optimized params ans solve diff eq with them
mean_p <- as.numeric(mean_param)
theta <- c('ka' = mean_p[1], 'Cl' = mean_p[2], 'Vd' = mean_p[3], 'adlag' = mean_p[4])
Frac0 <- mean_p[6]#0.729
tinf <- mean_p[5]
doses <- unique((plot_data2_dose[plot_data2_dose$TIME==0,]$Dose)*3.47*10^7/5813.68)

Arm_dose<- unique(plot_data2_dose$Arm_dose)
doses_df<- data.frame(doses, Arm_dose)

tim<-seq(0, 8, length.out=100)

doses_df$doses<-as.numeric(doses_df$doses)

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
  et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
  et(amt=Dose*(1 - Frac0), cmt = 1) %>%
  et(tim)
Cc<-as.numeric(unlist(rxSolve(lispro_1cmp, theta, ev)['Cc']))
Arm_dose<-rep(doses_df[i,2],100)
time<-tim
df_new<-data.frame(Cc,Arm_dose, time)
df<- rbind(df,df_new)
}

View(df)
#plot it 

ggplot()+geom_point(data = plot_data2_dose,
                    mapping = aes(x = TIME, 
                                  y = DV)) +
  facet_wrap(~Arm_dose)+
  #geom_point(alpha = .7, size = 2) +
  geom_point(data = plot_mean, mapping = aes(x = t_mean, y = mean_sol))+
  facet_wrap(~Arm_dose)+ 
  geom_line(data = df, mapping = aes(x = time, y =Cc ))+
  facet_wrap(~Arm_dose)+ 
  
labs(x = " time, h",
     y = "Serum Insulin, pmol/L")
#ggsave("pkdata_torstensolch1_v9.png",dpi=350,height=5,width=7.5)




#Explore solution with different sets of parameters

#define sol to plot : plot_mean and params to use: mean_param


x_stan<-data_out1[,189:510]
x2_stan <- x_stan[,seq(2, 322, by=2)]
data_out_tar1<-data_out1[,8:13]
print(data_out_tar1)
#x_stan<-data_out1[,174:467]
#x2_stan <- x_stan[,seq(2, 294, by=2)]
x2_stan_normV <- x2_stan / data_out_tar1$V_lis_hat
iter_ind<- 100

mean_sol <- as.numeric(cHat_stan[iter_ind,])#apply(cHat_stan,2,mean)
#mean_sol<-as.numeric(x2_stan_normV[iter_ind,])
t_mean <- plot_data2_dose$TIME
Arm_dose<-plot_data2_dose$Arm_dose
#View(mean_sol)
plot_mean <- data.frame(t_mean,mean_sol, Arm_dose)
#View(plot_mean)


mean_param<-data_out_tar1[iter_ind,]#apply(data_out_tar1,2,median)

print(mean_param)
mean_p <- as.numeric(mean_param)
theta <- c('ka' = mean_p[1], 'Cl' = mean_p[2], 'Vd' = mean_p[3], 'adlag' = mean_p[4])
Frac0 <- 0.729
tinf <- mean_p[5]
doses <- unique((plot_data2_dose[plot_data2_dose$TIME==0,]$Dose)*3.47*10^7/5813.68)

Arm_dose<- unique(plot_data2_dose$Arm_dose)
doses_df<- data.frame(doses, Arm_dose)

tim<-seq(0, 8, length.out=100)

doses_df$doses<-as.numeric(doses_df$doses)

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
    et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
    et(amt=Dose*(1 - Frac0), cmt = 1) %>%
    et(tim)
  Cc<-as.numeric(unlist(rxSolve(lispro_1cmp, theta, ev)['Cc']))
  Arm_dose<-rep(doses_df[i,2],100)
  time<-tim
  df_new<-data.frame(Cc,Arm_dose, time)
  df<- rbind(df,df_new)
}

#View(df)
#plot it 

ggplot()+geom_point(data = plot_data2_dose,
                    mapping = aes(x = TIME, 
                                  y = DV, color="data")) +
  facet_wrap(~Arm_dose)+
  geom_point(data = plot_mean, mapping = aes(x = t_mean, y = mean_sol, color="Torsten"))+
  facet_wrap(~Arm_dose)+ 
  geom_line(data = df, mapping = aes(x = time, y = Cc, color="diff eq"))+
  facet_wrap(~Arm_dose)+ 
labs(x = " time, h",
     y = "Serum Insulin, pmol/L")



#for one article




data_out1<- read.csv("output_2_103frac02_mod.csv")
#View(x_stan)
#x_stan<-data_out1[,47:84] #93
x_stan<-data_out1[,48:89] #103
#x2_stan <- x_stan[,seq(2, 38, by=2)] #93
x2_stan <- x_stan[,seq(2, 42, by=2)] #103
data_out_tar1<-data_out1[,8:14]
data_out_tar1
#x_stan<-data_out1[,174:467]
#x2_stan <- x_stan[,seq(2, 294, by=2)]
x2_stan_normV <- x2_stan / data_out_tar1$V_lis_hat
iter_ind<- 100
mean_sol<-apply(x2_stan_normV,2,median)#as.numeric(x2_stan_normV[iter_ind,])
#t_mean <- plot_data2_dose[plot_data2_dose$ID=='93_1', 'TIME']
#Arm_dose<-plot_data2_dose[plot_data2_dose$ID=='93_1', 'Arm_dose']
t_mean <- plot_data2_dose[plot_data2_dose$ID=='103_1', 'TIME']
Arm_dose<-plot_data2_dose[plot_data2_dose$ID=='103_1', 'Arm_dose']

#View(mean_sol)
plot_mean <- data.frame(t_mean,mean_sol, Arm_dose)
#View(plot_mean)


mean_param<-apply(data_out_tar1,2,mean)#data_out_tar1[iter_ind,]#apply(data_out_tar1,2,median)

print(mean_param)
mean_p <- as.numeric(mean_param)
theta <- c('ka' = mean_p[1], 'Cl' = mean_p[2], 'Vd' = mean_p[3], 'adlag' = 0)#mean_p[4])
Frac0 <- 0#0#0.729#mean_p[6]
tinf <- mean_p[5]
doses <- unique((plot_data2_dose[plot_data2_dose$TIME==0,]$Dose)*3.47*10^7/5813.68)

Arm_dose<- unique(plot_data2_dose$Arm_dose)
doses_df<- data.frame(doses, Arm_dose)

tim<-seq(0, 8, length.out=100)

doses_df$doses<-as.numeric(doses_df$doses)

Dose <-doses_df[1,1]
ev <- et(timeUnits="hr") %>%
  et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
  et(amt=Dose*(1 - Frac0), cmt = 1) %>%
  et(tim)
Cc<-as.numeric(unlist(rxSolve(lispro_1cmp, theta, ev)['Cc']))
Arm_dose<-rep(doses_df[1,2],100) 
time<-tim
df<-data.frame(Cc,Arm_dose,time)

#frac and 1-frac confused
for (i in 2:8) {
  Dose <-doses_df[i,1]
  ev <- et(timeUnits="hr") %>%
    et(amt=Dose*Frac0, dur = tinf, cmt = 2) %>%
    et(amt=Dose*(1 - Frac0), cmt = 1) %>%
    et(tim)
  Cc<-as.numeric(unlist(rxSolve(lispro_1cmp, theta, ev)['Cc']))
  Arm_dose<-rep(doses_df[i,2],100)
  time<-tim
  df_new<-data.frame(Cc,Arm_dose, time)
  df<- rbind(df,df_new)
}
#df
#View(df)
#plot it 
#df1<-df[df$Arm_dose=='93_1 6 U',]
#plot_data2_dose1<-plot_data2_dose[plot_data2_dose$ID=='93_1',]
df1<-df[df$Arm_dose=='103_1 10 U',]
plot_data2_dose1<-plot_data2_dose[plot_data2_dose$ID=='103_1',]

#View(df[df$Arm_dose=='93_1 10 U',])
ggplot()+geom_point(data = plot_data2_dose1,
                    mapping = aes(x = TIME, 
                                  y = DV, color="data")) +
  #facet_wrap(~Arm_dose)+
  geom_point(data = plot_mean, mapping = aes(x = TIME, y = mean_sol, color="Torsten"))+
  #facet_wrap(~Arm_dose)+ 
  geom_line(data = df1, mapping = aes(x = time, y = Cc, color="diff eq"))+
  #facet_wrap(~Arm_dose)+ 
  labs(x = " time, h",
       y = "Serum Insulin, pmol/L")

ggsave("insulin103_torsten_frac02.png",dpi=350,height=5,width=7.5)

View(plot_mean)




