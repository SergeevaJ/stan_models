##Script with all preprocessing steps to make convenient for R dataset with pk and pd data for insulin lispro from database  


#install.packages("readxl")
#install.packages("writexl")
#install.packages("dplyr")

library("readxl")
library("writexl")
library("dplyr")
library("ggplot2")

#getwd()
#setwd("C:/Users/Юлия/Desktop/Summer Internship 2024 M&S Decisions")


#aggregate dose, efficacy and BL data in one dataset

data_eff <- read_excel("DB_Insulin_Lispro_longBL.xlsx",  sheet = "Efficacy")

#typeof(data_eff)
#class(data_eff)
#names(data_eff)
#View(data_eff1)

data_eff1 <- data_eff[1:8]

data_bl <- read_excel("DB_Insulin_Lispro_longBL.xlsx",  sheet = "BL_unstructured")

#typeof(data_bl)
#class(data_bl)
#View(data_bl)
#names(data_bl)


arms<- unique(data_eff1$Arm_ID)
arms

col_names<-names(data_bl)
col_names<-as.list(col_names)[2:17]

for (col in col_names){
  vect<- c()
  for (arm in arms){
  n_items <-nrow(filter(data_eff1, Arm_ID==arm))
  vect <-append(vect, c(rep(filter(data_bl, Arm_ID==arm)[col],n_items)))}
  data_eff1$col<-vect
  names(data_eff1)[names(data_eff1) == 'col'] <- col
}

#View(data_eff1)


#add number of patients
vect<- c()
for (arm in arms){
  n_items <-nrow(filter(data_eff1, Arm_ID==arm))
  vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Nsub,n_items)))}
data_eff1$Nsub<-vect

View(data_eff1)

names(data_eff1)



names(data_eff1)<- c("Arm_ID","TIME","TIME_Unit", "DV","NAME", "UNIT",                   
 "TYPE", "MEAS",                     
 "Sex_M","Sex_F",                    
 "BMI","BMI_SD",           
 "Weight","Weight_SD",            
 "Age","Age_SD",            
 "Duration_of_T1D", "Duration_of_T1D_SD",
 "Meal_description", "Meal_kcal",                
 "Meal_carbohydrates", "Meal_protein",             
 "Meal_fat", "Comments_B", "Nsub")

col_names3<-c("Sex_M","Sex_F",                    
              "BMI","BMI_SD",           
              "Weight","Weight_SD",            
              "Age","Age_SD",            
              "Duration_of_T1D", "Duration_of_T1D_SD",
              "Meal_kcal",                
              "Meal_carbohydrates", "Meal_protein",             
              "Meal_fat", "Nsub")

for (col in col_names3){
  data_eff1[col]<-as.numeric(unlist(data_eff1[col]))}




#make all data in the same units

data_eff1[data_eff1$TIME_Unit=='min', 'TIME']<- data_eff1[data_eff1$TIME_Unit=='min', 'TIME']/60

data_eff1[data_eff1$TIME_Unit=='min', 'TIME_Unit']<- c(rep('h', nrow(data_eff1[data_eff1$TIME_Unit=='min', 'TIME']) ))


data_eff1[data_eff1$UNIT=='mg/dL', 'DV']<-data_eff1[data_eff1$UNIT=='mg/dL', 'DV']*10000/180156

data_eff1[data_eff1$UNIT=='mg/dL', 'UNIT']<- c(rep('mmol/L', nrow(data_eff1[data_eff1$UNIT=='mg/dL', 'UNIT'])))


data_eff1[data_eff1$UNIT=='mU/L', 'DV']<-data_eff1[data_eff1$UNIT=='mU/L', 'DV']*34700/5813.68

data_eff1[data_eff1$UNIT=='mU/L', 'UNIT']<- c(rep('pmol/L', nrow(data_eff1[data_eff1$UNIT=='mU/L', 'UNIT'])))


#unique(data_eff1["UNIT"])


#make every start time point = 0 h 

data_eff1[data_eff1$Arm_ID=='93_1', 'TIME'] <- data_eff1[data_eff1$Arm_ID=='93_1', 'TIME']-7

data_eff1[data_eff1$Arm_ID=='103_1', 'TIME'] <- data_eff1[data_eff1$Arm_ID=='103_1', 'TIME']-7.5

data_eff1[data_eff1$Arm_ID=='105_1', 'TIME'] <- data_eff1[data_eff1$Arm_ID=='105_1', 'TIME']-7


View(data_eff1[data_eff1$Arm_ID=='105_1', 'DV'])

#delete 1 one point corresponding to pre-bolus insulin and glucose in arm 93_1

data_eff1 <- data_eff1[data_eff1$Arm_ID!='93_1' | (data_eff1$Arm_ID=='93_1'& data_eff1$TIME>0),]

View(data_eff1[data_eff1$Arm_ID=='93_1', 'TIME'])

#shift time 

zero_time <- as.numeric(data_eff1[data_eff1$Arm_ID=='93_1', 'TIME'][1,])
data_eff1[data_eff1$Arm_ID=='93_1',"TIME"] <- data_eff1[data_eff1$Arm_ID=='93_1',"TIME"]-zero_time


#start time 1 h earlier for 105_1

data_eff1 <- data_eff1[data_eff1$Arm_ID!='105_1' | (data_eff1$Arm_ID=='105_1'& data_eff1$TIME>=1),]

data_eff1[data_eff1$Arm_ID=='105_1', 'TIME'] <- data_eff1[data_eff1$Arm_ID=='105_1', 'TIME']-1


#set all insulin pk curves to start at zero value

#25, 86, 93, 102, 103, 105
as.numeric(data_eff1[data_eff1$Arm_ID=='25_1' & data_eff1$NAME=='LISPRO_PK', 'DV'][1,])


arms_corr_ins <-c('25_1', '93_1', '102_1', '103_1', '105_1')

for (arm in arms_corr_ins){
   zero_ins <- as.numeric(data_eff1[data_eff1$Arm_ID==arm & data_eff1$NAME=='LISPRO_PK', 'DV'][1,])
   data_eff1[data_eff1$Arm_ID==arm & data_eff1$NAME=='LISPRO_PK', 'DV'] <- data_eff1[data_eff1$Arm_ID==arm& data_eff1$NAME=='LISPRO_PK', 'DV']-zero_ins
}

#for 86_1 the minimal insulin concentration is subtracted
zero_ins <- as.numeric(data_eff1[data_eff1$Arm_ID=='86_1' & data_eff1$NAME=='LISPRO_PK', 'DV'][12,])
data_eff1[data_eff1$Arm_ID=='86_1' & data_eff1$NAME=='LISPRO_PK', 'DV'] <- data_eff1[data_eff1$Arm_ID=='86_1'& data_eff1$NAME=='LISPRO_PK', 'DV']-zero_ins


View(data_eff1[data_eff1$Arm_ID=='105_1', 'DV'])

# set all glucose pd curves to start at non-zero value
data_eff1[data_eff1$Arm_ID=='19_1' & data_eff1$NAME=='LISPRO_PD', 'DV'] <- data_eff1[data_eff1$Arm_ID=='19_1'& data_eff1$NAME=='LISPRO_PD', 'DV']+7.5



#remove data after lunch and dinner for 105_1

View(data_eff1[data_eff1$Arm_ID!='105_1' | (data_eff1$Arm_ID=='105_1'& data_eff1$TIME<=4),])

data_eff1 <- data_eff1[data_eff1$Arm_ID!='105_1' | (data_eff1$Arm_ID=='105_1'& data_eff1$TIME<=4),]

#remove negative points from 25 curve from consideration

data_eff1 <- data_eff1[data_eff1$Arm_ID!='25_1' | (data_eff1$Arm_ID=='25_1'& data_eff1$DV>=0),]

#remove 86_1 and 25_1 arms from dataset

#data_eff1 <- data_eff1[data_eff1$Arm_ID!='25_1' & data_eff1$Arm_ID!='86_1',]

#delete the last point from 103_1 because it is negative 
data_eff1<-data_eff1[(data_eff1$Arm_ID!='103_1')|(data_eff1$Arm_ID=='103_1')&(data_eff1$TIME!=5.5),]




data_dose <- read_excel("DB_Insulin_Lispro_longBL.xlsx",  sheet = "Arms")
View(data_dose)

#vect<- c()
#for (arm in arms){
#  n_items <-nrow(filter(data_eff1, Arm_ID==arm))
#  vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Dose,n_items)))}
#data_eff1$Dose<-vect

#vect<- c()
#for (arm in arms){
#  n_items <-nrow(filter(data_eff1, Arm_ID==arm))
#  vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Dose_SD,n_items)))}
#data_eff1$Dose_SD<-vect

#make a column which distinguishes dose administration events and measurements(pk,pd)
#ADM: 0 - no administration, 1 - administration
data_eff1$ADM <- c(rep(0, nrow(data_eff1)))

#View(data_eff1)

#make rows with dose

# Function to create new dataframe 
insertRow <- function(data, new_row, r) { 
  if (r==0) {data_new<- rbind(new_row, data)}
  else
  {data_new <- rbind(data[1:r, ],             
                     new_row,                 
                     data[- (1:r), ])}
  
  return(data_new) 
} 


#View(data_eff1[1,6:26])

for (arm in arms){
  dose_row<- c(arm, 0, 'h', 0,'LISPRO_PK', as.character(data_eff1[data_eff1$Arm_ID==arm,][1,6:25]), 1)
  ind <- which((data_eff1$Arm_ID==arm&data_eff1$NAME=='LISPRO_PK'& data_eff1$TIME==0),arr.ind=T)
  print(ind)
  data_eff1 <- insertRow(data_eff1,dose_row,ind-1)}
View(data_eff1)


#add dose amount

vect<- c()
for (arm in arms){
  n_items <-nrow(filter(data_eff1, Arm_ID==arm))
  vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Dose,n_items)))}
data_eff1$Dose<-vect

data_eff1$Dose<-as.numeric(data_eff1$Dose)
data_eff1$ADM<-as.numeric(data_eff1$ADM)

data_eff1 %>%
  select(names(data_eff1)) %>%
  mutate(AMT= Dose*ADM) -> data_eff1



data_eff1 <- subset(data_eff1, select = -Dose)

#View(data_eff1)



##Read existing dataset for plots

data_eff1 <- read_excel("dataset_lispro_with86_Dosepmol_v16.xlsx")
data_eff1<- data_eff1[2:31]
View(data_eff1)



##for labels with dose add Dose column


data_dose <- read_excel("DB_Insulin_Lispro_longBL.xlsx",  sheet = "Arms")
View(data_dose)


arms<- unique(data_eff1$ID)
vect<- c()
for (arm in arms){
  n_items <-nrow(filter(data_eff1, ID==arm))
  vect <-append(vect, c(rep(filter(data_dose, Arm_ID==arm)$Dose,n_items)))}
data_eff1$Dose<-vect


#View(data_eff1)



#plot all glucose curves

#library("ggplot2")
plot_data1 <- filter(data_eff1, DVNAME=='LISPRO_PD')
#plot_data1 <- filter(data_eff1, NAME=='LISPRO_PD')

plot_data1 %>%
  select(TIME, DV, Dose, ID) %>%
  mutate(Arm_dose= paste(paste(ID, as.character(Dose)), 'U')) -> plot_data1_dose

#plot_data1 %>%
#  select(TIME, DV, Dose, Arm_ID) %>%
#  mutate(Arm_dose= paste(paste(Arm_ID, as.character(Dose)), 'U')) -> plot_data1_dose

View(plot_data1_dose)


ggplot(data = plot_data1_dose,
       mapping = aes(x = TIME, 
                     y = DV,
                     color = Arm_dose)) +
  geom_point(alpha = .7,
             size = 2) +
labs(title = "Mean glucose concentration versus time after single s.c. dose of insulin lispro",
     x = " time, h",
     y = "Blood Glucose, mmol/L") #+theme_minimal()

ggsave("blood_glucose_v6_add86_dose.png",dpi=350,height=5,width=7.5)

#plot all insulin curves

plot_data2 <- data_eff1[data_eff1$DVNAME=='LISPRO_PK',]

#plot_data2 <- data_eff1[data_eff1$NAME=='LISPRO_PK',]

#ggplot(data = plot_data2,
#       mapping = aes(x = TIME, 
#                     y = DV,
#                     color = Arm_ID)) +
#  geom_point(alpha = .7,
#             size = 2) +
#  labs(title = "Mean insulin concentration versus time after single s.c. dose of insulin lispro",
#       x = " time, h",
#       y = "Serum Insulin, pmol/L")#+theme_minimal()

ggplot(data = plot_data2,
       mapping = aes(x = TIME, 
                     y = DV,
                     color = ID)) +
  geom_point(alpha = .7,
             size = 2) +
  labs(title = "Mean insulin concentration versus time after single s.c. dose of insulin lispro",
       x = " time, h",
       y = "Serum Insulin, pmol/L")#+theme_minimal()

ggsave("serum_insulin_with86_v10.png",dpi=350,height=5,width=7.5)

#plot with dose

plot_data2 %>%
  select(TIME, DV, Dose, ID) %>%
  mutate(Arm_dose= paste(paste(ID, as.character(Dose)), 'U')) -> plot_data2_dose


ggplot(data = plot_data2_dose,
       mapping = aes(x = TIME, 
                     y = DV,
                     color = Arm_dose)) +
  geom_point(alpha = .7,
             size = 2) +
  labs(title = "Mean insulin concentration versus time after single s.c. dose of insulin lispro",
       x = " time, h",
       y = "Serum Insulin, pmol/L")
ggsave("serum_insulin_with86_v10_dose.png",dpi=350,height=5,width=7.5)



#plot dose normalized curves
plot_data2$Dose <- as.numeric(plot_data2$Dose)
plot_data2_dose$Dose <- as.numeric(plot_data2_dose$Dose)

#plot_data2 %>%
#  select(TIME, DV, Dose, Arm_ID) %>%
#  mutate(Norm_DV = DV/Dose) -> plot_data2_norm


plot_data2_dose %>%
    select(TIME, DV, Dose, Arm_dose) %>%
    mutate(Norm_DV = DV/Dose) -> plot_data2_norm
  
ggplot(data = plot_data2_norm,
       mapping = aes(x = TIME, 
                     y = Norm_DV,
                     color = Arm_dose)) +
  geom_point(alpha = .7,
             size = 2) +
  labs(title = "Mean dose-nomalized insulin concentration versus time after single s.c. dose of insulin lispro",
       x = " time, h",
       y = "Serum Insulin, pmol/(L*U)")#+theme_minimal()


#ggplot(data = plot_data2_norm,
#       mapping = aes(x = TIME, 
#                     y = Norm_DV,
#                     color = Arm_ID)) +
#  geom_point(alpha = .7,
#             size = 2) +
#  labs(title = "Mean dose-nomalized insulin concentration versus time after single s.c. dose of insulin lispro",
#       x = " time, h",
#       y = "Serum Insulin, pmol/(L*U)")#+theme_minimal()

ggsave("serum_insulin_norm_with86_v5_dose.png",dpi=350,height=5,width=7.5)




View(data_eff1)


#if you used your dataset, better delete Dose column
#data_eff1 <- subset(data_eff1, select = -Dose)


#save our dataset

#reorder and rename columns

#Needed columns: ID, TIME, DV, DVID, DVNAME, CMT, ADM, AMT, EVID, MDV
#Add: DVID=0 or 1(== DVNAME), CMT==1, EVID==MDV==ADM
#Rename ArmID->ID, NAME -> DVNAME, 


data_eff2<- data_eff1
data_eff2$DVID <- as.numeric(data_eff2$NAME=='LISPRO_PK')

data_eff2$CMT <- c(rep(1, nrow(data_eff2)))

data_eff2$EVID <- data_eff2$ADM
data_eff2$MDV <- data_eff2$ADM


names(data_eff2)


new_cols <- c("Arm_ID","TIME","TIME_Unit", "DV", "DVID", "NAME",  "CMT", "ADM", "AMT", "EVID", "MDV", "UNIT","TYPE", "MEAS", "Sex_M",             
 "Sex_F","BMI","BMI_SD","Weight","Weight_SD","Age",               
 "Age_SD", "Duration_of_T1D", "Duration_of_T1D_SD",
 "Meal_description","Meal_kcal", "Meal_carbohydrates", "Meal_protein", "Meal_fat", "Comments_B", "Nsub")


data_eff2 <- data_eff2[,new_cols]

col_names4 <-c("DVID", "CMT", "EVID", "MDV")
for (col in col_names4){
  
  data_eff2[col]<-as.numeric(unlist(data_eff2[col]))}


names(data_eff2) <- c("ID","TIME","TIME_Unit", "DV", "DVID", "DVNAME",  "CMT", "ADM", "AMT", "EVID", "MDV", "UNIT","TYPE", "MEAS", "Sex_M",             
                      "Sex_F","BMI","BMI_SD","Weight","Weight_SD","Age",               
                      "Age_SD", "Duration_of_T1D", "Duration_of_T1D_SD",
                      "Meal_description","Meal_kcal", "Meal_carbohydrates", "Meal_protein", "Meal_fat", "Comments_B", "Nsub")
#View(data_eff2)


data_eff2<- subset(data_eff2, select= -Comments_B)

data_eff2$Meal_description <- as.character(data_eff2$Meal_description)


#writexl::write_xlsx(data_eff2, path= 'dataset_lispro_with86_v15.xlsx')
#write.csv(data_eff2, "dataset_lispro_with86_v15.csv")
