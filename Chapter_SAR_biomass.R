#Chapater SAR biomass mapping
#Xiaoxuan Li 
library(lidR)
library(tools) 
library(raster)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(ggcorrplot)
library(ggpmisc)
library(dplyr)
library(tidyverse) 
library(rsq)
library(readxl)
library(grid)
library(mltools)
library(Metrics)
library(devtools)
library(ggsignif)
library(igraph)
library(ff)
library(corrplot)
library(psych)
library(xlsx)
library(car)
library(MASS)
library(ggrepel)
library(rlist)
library(data.table)
library(lemon)
library(lmtest)
library(reshape2)
library(ggpubr)
library(tidyr)
library(stats)
library(igraph)
library(caret)
library(vegan)
library(rasterdiv)
library(snow)

# ------------------------------------------------------------------------------------------------ #
list_ini <- 0
list_ini2 <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table"
sar_date <- "_20170706_75m.xlsx"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")
#list_dir <- c("Welverdiendt")

dir_hcc_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_25m.csv"
data_hcc_sar = read.csv(dir_hcc_sar)

for (l in list_dir){
  print(l)
  dir_index <- file.path(dir,paste0(l,"_index_25m.csv"))
  Data_index = read.csv(dir_index)
  colnames(Data_index) <- c("xx","ID_25m","ID_75m")
  
  dir_chm <- file.path(dir,paste0(l,"_chm_25m.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  #Data_chm <- rbind(list_ini_chm, Data_chm)

  dir_cover <- file.path(dir,paste0(l,"_cc_25m.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  #Data_cover <- rbind(list_ini_cover, Data_cover)

  dir_sar <- file.path(dir,paste0(l,"_sar_25m.csv"))
  Data = read.csv(dir_sar)
  Data_sar = subset(Data, select=c(2,5))
  #Data_sar <- rbind(list_ini_sar, Data_sar)
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_index, Data_chm, by.x = "ID_25m", by.y = "FID")
  Data_index <- merge(Data_index, Data_cover, by.x = "ID_25m", by.y = "FID")
  Data_index <- merge(Data_index, Data_sar, by.x = "ID_25m", by.y = "FID")
  colnames(Data_index) <- c("xx","ID_25m","ID_75m","CHM","Cover","SAR")
  Data_index <- na.omit(Data_index)
  
  #data conversion, 75m agbd cal
  Data_index$Volume <- Data_index$CHM*625
  Data_index$Cover <- Data_index$Cover/625
  Data_index$AGBD <- 9.0665*Data_index$CHM*Data_index$Cover*0.0625/0.5625
  #Data_index$AGBD <- 4.852459*exp(0.3274118*Data_index$CHM)/0.5625
  Data_index$HCC <- Data_index$CHM*Data_index$Cover
  
  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
  Data_index <- Data_index[Data_index$SAR > -25 & Data_index$SAR < -10,]
  
  #75m aggregation
  df_agbd <- aggregate(AGBD ~ ID_75m, data = Data_index, FUN = sum)
  df_sar <- aggregate(SAR ~ ID_75m, data = Data_index, FUN = mean)
  df_agbd_sar <- cbind(df_agbd, df_sar)

  #merge by site
  list_ini <- rbind(list_ini, df_agbd_sar)

  df_hcc <- data.frame(Data_index$HCC)
  df_sar <- data.frame(Data_index$SAR)
  df_hcc_sar <- cbind(df_hcc,df_sar)
  list_ini2 <- rbind(list_ini2, df_hcc_sar)
}

#file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_25m.csv")
#write.csv(list_ini2, file = file_path)

Data <- list_ini[list_ini$AGBD > 1,]

#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR))
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)

summary(Data$AGBD)
Data = subset(Data, select=c(2,4))
ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 90),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  labs(x="ALS AGBD (Mg/ha)", 
       y="Mean of ScanSAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


#Cover
reg2 <- lm(log(Cover) ~ SAR, data = Data)
pred <- exp((coef(reg2)[1] + coef(reg2)[2]*Data$SAR))
coef(reg2)[2]
coef(reg2)[1]
round(summary(reg2)$adj.r.squared,2)
mean(pred - Data$Cover)
sqrt(mean((pred - Data$Cover)^2))
100*mean(pred - Data$Cover)/mean(Data$Cover)
100*sqrt(mean((pred - Data$Cover)^2))/mean(Data$Cover)




# ------------------------------------------------------------------------------------------------ #
#output ALS volume, cover, and SAR backscatter in one csv file based on SAR date 
#update version 20230907
#fine beam only
list_ini_chm <- 0
list_ini_cover <- 0
list_ini_sar <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_0625\\site_75m"
sar_date <- "_20170706_75m.xlsx"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
#fine beam 4 sites only
list_dir <- c("Agincourt","Ireagh","Justicia","Welverdiendt")
for (l in list_dir){
  print(l)
  dir_index <- file.path(dir,l,"index.csv")
  Data_index = read.csv(dir_index)
  colnames(Data_index)[1] <- "ID"
  
  dir_chm <- file.path(dir,l,paste0("CHM_",l,"_75m.xlsx"))
  Data = read_excel(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  Data_chm <- rbind(list_ini_chm, Data_chm)
  Data_index <- merge(Data_index, Data_chm, by.x = "ID", by.y = "FID")
  
  dir_cover <- file.path(dir,l,paste0("Cover_",l,"_75m.xlsx"))
  Data = read_excel(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  Data_cover <- rbind(list_ini_cover, Data_cover)
  Data_index <- merge(Data_index, Data_cover, by.x = "ID", by.y = "FID")
  
  dir_sar <- file.path(dir,l,paste0("FineBeam_",l,sar_date))
  Data = read_excel(dir_sar)
  Data_sar = subset(Data, select=c(2,5))
  Data_sar <- rbind(list_ini_sar, Data_sar)
  Data_index <- merge(Data_index, Data_sar, by.x = "ID", by.y = "FID")

}


Data <- na.omit(Data_index)
colnames(Data) <- c("ID","CHM","Cover","SAR")
Data$Volume <- Data$CHM*5625
Data$Cover <- Data$Cover/5625
#Data$AGBD <- 10.35*Data$CHM*Data$Cover - 5.9236
Data$AGBD <- 9.0665*Data$CHM*Data$Cover

Data <- Data[is.finite(rowSums(Data)),]
Data <- Data[Data$CHM > 0,]
Data <- Data[Data$Cover > 0,]
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]


#Cover
reg2 <- lm(log(Cover) ~ SAR, data = Data)
pred <- exp((coef(reg2)[1] + coef(reg2)[2]*Data$SAR))
coef(reg2)[2]
coef(reg2)[1]
round(summary(reg2)$adj.r.squared,2)
mean(pred - Data$Cover)
sqrt(mean((pred - Data$Cover)^2))
100*mean(pred - Data$Cover)/mean(Data$Cover)
100*sqrt(mean((pred - Data$Cover)^2))/mean(Data$Cover)

#AGBD
Data$AGBD <- Data$AGBD 
reg3 <- lm(AGBD ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
Data$AGBD <- Data$AGBD 
pred <- pred 
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)






# ------------------------------------------------------------------------------------------------ #
#output GEDI volume, cover, and SAR backscatter in one csv file based on SAR date 
#update version 20230907

#20140906
#20160903
#20170902
#20180929
#20190928
#20200926
#20210925
#20220924

sar_date <- "20180929.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_0906\\GEDI_predictions\\Table"
dir_sar <- file.path(dir,paste0("ScanSAR_GEDI_",sar_date))
Data = read.csv(dir_sar)
Data = subset(Data, select=c("RASTERVALU","AGBD"))
colnames(Data) <- c("SAR","AGBD")
Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]


#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)
n<- nrow(Data)





ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 90),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  labs(x="GEDI SAS AGBD (Mg/ha)", 
       y="Mean of ScanSAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  #theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))



sar_date <- "20170702.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_0906\\GEDI_predictions\\Table"
dir_sar <- file.path(dir,paste0("FineBeam_GEDI_",sar_date))
Data = read.csv(dir_sar)
Data <- na.omit(Data)
Data = subset(Data, select=c("RASTERVALU","AGBD"))
colnames(Data) <- c("SAR","AGBD")
Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]

#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)

#Exp(0.452008*"PALSAR2_ScanSAR_HV_mtf_5_db_2018-09-29.tif"+11.09319) 



# ------------------------------------------------------------------------------------------------ #
list_ini <- 0
list_ini2 <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table"
sar_date <- "_20170706_75m.xlsx"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")
for (l in list_dir){
  print(l)

  dir_sar <- file.path(dir,paste0(l,"_sar_1ha.csv"))
  Data = read.csv(dir_sar)
  Data_index = subset(Data, select=c(1,5))
  colnames(Data_index) <- c("ID_1ha","SAR")
  
  dir_chm <- file.path(dir,paste0(l,"_chm_1ha.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))

  dir_cover <- file.path(dir,paste0(l,"_cc_1ha.csv"))
  Data = read.csv(dir_cover)
  Data <- Data[Data$AREA >9900,]
  Data_cover = subset(Data, select=c(2,5))

  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_index, Data_chm, by.x = "ID_1ha", by.y = "FID")
  Data_index <- merge(Data_index, Data_cover, by.x = "ID_1ha", by.y = "FID")
  colnames(Data_index) <- c("ID_1ha","SAR","CHM","Cover")
  Data_index <- na.omit(Data_index)

  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
  Data_index <- Data_index[Data_index$SAR > -25 & Data_index$SAR < -10,]
  Data_index$Cover <- Data_index$Cover/10000
  Data_index$HCC <- Data_index$CHM*Data_index$Cover
  
  #data merge
  list_ini2 <- rbind(list_ini2, Data_index)
}

file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_1ha.csv")
write.csv(list_ini2, file = file_path)



# ------------------------------------------------------------------------------------------------ #
list_ini <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")

for (l in list_dir){
  print(l)
  dir_index <- file.path(dir,paste0(l,"_GEDI_1ha.csv"))
  Data_index = read.csv(dir_index)
  Data_index <- Data_index[Data_index$Join_Count >1,]
  Data_index = Data_index[,c("OID_","sensitivit","GEDI_Cover","GEDI_FHD",
                             "RH1_50","RH1_75","RH1_90","RH1_98","AGBD")]
  
  dir_chm <- file.path(dir,paste0(l,"_GEDI_chm_1ha.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_GEDI_cc_1ha.csv"))
  Data = read.csv(dir_cover)
  Data <- Data[Data$AREA >9900,]
  Data_cover = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_index, Data_chm, by.x = "OID_", by.y = "FID")
  Data_index <- merge(Data_index, Data_cover, by.x = "OID_", by.y = "FID")
  colnames(Data_index) <- c("ID_1ha","sensitivit","GEDI_Cover","GEDI_FHD",
                            "RH1_50","RH1_75","RH1_90","RH1_98",
                            "AGBD","CHM","Cover")
  Data_index <- na.omit(Data_index)
  Data_index$Cover <- Data_index$Cover/10000
  Data_index$HCC <- Data_index$CHM*Data_index$Cover
  
  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
  
  list_ini <- rbind(list_ini, Data_index)
}

file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_gedi_1ha.csv")
write.csv(list_ini, file = file_path)




# ------------------------------------------------------------------------------------------------ #
#update 10/25/2023
list_ini <- 0
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")

for (l in list_dir){
  dir_chm <- file.path(dir,paste0("ALS_1ha_",l,"_CHM.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0("ALS_1ha_",l,"_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_chm <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  colnames(Data_chm) <- c("FID","CHM","Cover")
  
  #data conversion, 1ha agbd cal
  Data_chm$Cover <- Data_chm$Cover/10000
  Data_chm$AGBD <- 18.839*Data_chm$CHM*Data_chm$Cover -1.213
  
  #merge by site
  list_ini <- rbind(list_ini, Data_chm)
}

file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\ALS_1ha_AGBD_HCC.csv")
write.csv(list_ini, file = file_path)




# ------------------------------------------------------------------------------------------------ #
#update version 20231025
#update version 01212024 1ha ALS~AGBD modeling 
#the equation coefficient can be used to map cover/agbd using SAR backscatter
#2014
#2016
#2017
#2018
#2019
#2020
#2021
#2022
#2023
#sar
sar_date <- "ALS_1ha_2014.csv"
#sar_date <- "Finebeam_1ha.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
dir_sar <- file.path(dir,sar_date)
Data = read.csv(dir_sar)
Data = subset(Data, select=c("FID","MEAN"))
colnames(Data) <- c("FID","SAR")

#als
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\ALS_1ha_AGBD_HCC.csv"
Data_als = read.csv(dir)

Data <- merge(Data, Data_als, by.x = "FID", by.y = "FID")

Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$Cover > 0,]
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]

#Cover
reg2 <- lm(log(Cover) ~ SAR, data = Data)
pred <- exp((coef(reg2)[1] + coef(reg2)[2]*Data$SAR))
coef(reg2)[2]
coef(reg2)[1]
round(summary(reg2)$adj.r.squared,2)
mean(pred - Data$Cover)
sqrt(mean((pred - Data$Cover)^2))
100*mean(pred - Data$Cover)/mean(Data$Cover)
100*sqrt(mean((pred - Data$Cover)^2))/mean(Data$Cover)

#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)


bias<-mean(pred - Data$AGBD)
RMSE<-sqrt(mean((pred - Data$AGBD)^2))
rbias<-100*mean(pred - Data$AGBD)/mean(Data$AGBD)
rRMSE<-100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)
Data$preds <- pred

p1 <- ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=61, y=-20,hjust = 0,size = 7,
           label= paste(expression(Logarithmic~model~R^2),": ",
                        round(summary(reg3)$adj.r.squared,3)), parse=TRUE) + 
  annotate("text",x=60,y=-21.5,hjust = 0,size = 7,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(bias,3)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,3)," %")) + 
  annotate("text",x=5,y=-10,hjust = 0,size = 15,family= "A", label= "(a)") + 
  labs(x="Area-based (H*CC) ALS AGBD (Mg/ha)", 
       y="Mean of SAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))



# ------------------------------------------------------------------------------------------------ #
#update 10/25/2023
list_ini <- 0
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
list_dir <- list.dirs(dir, recursive=FALSE,full.names = FALSE)
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")

for (l in list_dir){
  dir_chm <- file.path(dir,paste0("GEDI_1ha_",l,"_CHM.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0("GEDI_1ha_",l,"_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_chm <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  colnames(Data_chm) <- c("FID","CHM","Cover")
  
  #data conversion, 1ha agbd cal
  Data_chm$Cover <- Data_chm$Cover/10000
  Data_chm$AGBD <- 18.839*Data_chm$CHM*Data_chm$Cover -1.213
  
  #merge by site
  list_ini <- rbind(list_ini, Data_chm)
}

file_path <- file.path("E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\GEDI_1ha_AGBD_HCC.csv")
write.csv(list_ini, file = file_path)



# ------------------------------------------------------------------------------------------------ #
#update version 20231025
#GEDI

#2014
#2016
#2017
#2018
#2019
#2020
#2021
#2022
#sar
sar_date <- "GEDI_1ha_2017.csv"
sar_date <- "Finebeam_1ha_gedi.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
dir_sar <- file.path(dir,sar_date)
Data_sar = read.csv(dir_sar)
#Data_sar = subset(Data_sar, select=c("OID_","MEAN"))
Data_sar = subset(Data_sar, select=c("FID","MEAN"))
colnames(Data_sar) <- c("FID","SAR")
#GEDI merged 
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\GEDI_footprint_GEDI_1ha_AGBD_HCC.csv"
Data_gedi = read.csv(dir)
#L4A
Data <- merge(Data_sar, Data_gedi, by.x = "FID", by.y = "FID")
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]
Data <- Data[Data$Cover > 0,]
Data <- Data[Data$AGBD > 0,]
Data$AGBD <- Data$GEDI_AGBD
#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 
coef(reg3)[2]
coef(reg3)[1]
round(summary(reg3)$adj.r.squared,2)
mean(pred - Data$AGBD)
sqrt(mean((pred - Data$AGBD)^2))
100*mean(pred - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)


# ------------------------------------------------------------------------------------------------ #
#2014
#2016
#2017
#2018
#2019
#2020
#2021
#2022

#sar
#sar_date <- "GEDI_1ha_2022.csv"
sar_date <- "Finebeam_1ha_gedi.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
dir_sar <- file.path(dir,sar_date)
Data_sar = read.csv(dir_sar)
Data_sar = subset(Data_sar, select=c("FID","MEAN"))
#Data_sar = subset(Data_sar, select=c("OID_","MEAN"))
colnames(Data_sar) <- c("FID","SAR")

#SAS - RF
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\GEDI_footprint_GEDI_1ha_AGBD_HCC.csv"
Data = read.csv(dir)
Data <- Data[Data$AGBD >= 0,]
Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1
random_sample <- createDataPartition(Data$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data[random_sample, ]
testing_dataset <- Data[-random_sample, ]

#Local updated GEDI SAS RF
repeat_cv <- trainControl(method='cv', number=5)
reg_rf <- train(
  AGBD~sensitivit+
    GEDI_Cover + GEDI_FHD + 
    RH1_50 + RH1_75 + RH1_90 + RH1_98,
  data=training_dataset, 
  method='rf', 
  ntree=100,
  nodesize = 4,
  trControl=repeat_cv)

testing_dataset$preds <- predict(object=reg_rf, 
                                 newdata=testing_dataset[, ])

Data_gedi <- testing_dataset
Data <- merge(Data_sar, Data_gedi, by.x = "FID", by.y = "FID")
hist(Data$preds)
#Data$AGBD <- Data$preds
#AGBD
reg3 <- lm(log(preds) ~ SAR, data = Data)
round(summary(reg3)$adj.r.squared,2)
mean(Data$preds - Data$AGBD)
sqrt(mean((Data$preds - Data$AGBD)^2))
100*mean(Data$preds - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((Data$preds - Data$AGBD)^2))/mean(Data$AGBD)

R2 <- 0.44
RMSE<- 12.48
rRMSE<- 56.3
bias<- 0.19
rbias<- 0.9
  
p2 <- ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=61, y=-20,hjust = 0,size = 7,
           label= paste(expression(Logarithmic~model~R^2),": ",
                        round(R2,3)), parse=TRUE) + 
  annotate("text",x=60,y=-21.5,hjust = 0,size = 7,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(bias,3)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,3)," %")) + 
  annotate("text",x=5,y=-10,hjust = 0,size = 15,family= "A", label= "(b)") + 
  labs(x="GEDI (SAS RF)-estimated AGBD (Mg/ha)", 
       y="Mean of ScanSAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(p1,p2)
out = "E:\\ScanSAR\\ScanSAR\\Result_1004\\figure\\Scantterplot_SAR_AGBD.jpg"
ggsave(out,height=12, width=24, dpi=600)




# ------------------------------------------------------------------------------------------------ #
#2014
#2016
#2017
#2018
#2019
#2020
#2021
#2022

#sar
sar_date <- "GEDI_1ha_2017.csv"
sar_date <- "Finebeam_1ha_gedi.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
dir_sar <- file.path(dir,sar_date)
Data_sar = read.csv(dir_sar)
Data_sar = subset(Data_sar, select=c("FID","MEAN"))
#Data_sar = subset(Data_sar, select=c("OID_","MEAN"))
colnames(Data_sar) <- c("FID","SAR")

#SAS - GLM
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\GEDI_footprint_GEDI_1ha_AGBD_HCC.csv"
Data = read.csv(dir)
Data <- Data[Data$AGBD >= 0,]
Data$RH1_50 <- Data$RH1_50 +100
Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 150, by=1))) %>% 
  arrange(as.numeric(gr))
i <- 1
random_sample <- createDataPartition(Data$AGBD, p = 0.7, list = FALSE)
training_dataset  <- Data[random_sample, ]
testing_dataset <- Data[-random_sample, ]

#Local updated GEDI SAS GLM
repeat_cv <- trainControl(method='cv', number=5)
reg_glm <- train(
  AGBD ~ GEDI_FHD + RH1_50 + RH1_98, 
  data=training_dataset, 
  method='glm', 
  family = Gamma(link = "identity"),
  start=c(0,0.1,100,5),
  trControl=repeat_cv)

testing_dataset$preds <- predict(object=reg_glm, 
                                 newdata=testing_dataset[, ])

Data_gedi <- testing_dataset
Data <- merge(Data_sar, Data_gedi, by.x = "FID", by.y = "FID")
hist(Data$preds)
#Data$AGBD <- Data$preds
#AGBD
reg3 <- lm(log(preds) ~ SAR, data = Data)
round(summary(reg3)$adj.r.squared,2)
mean(Data$preds - Data$AGBD)
sqrt(mean((Data$preds - Data$AGBD)^2))
100*mean(Data$preds - Data$AGBD)/mean(Data$AGBD)
100*sqrt(mean((Data$preds - Data$AGBD)^2))/mean(Data$AGBD)







# ------------------------------------------------------------------------------------------------ #
#make histogram of biomass map 

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_AGBD_als_Clip.tif"
df_agbd <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_agbd) <- c("agbd")
df_agbd <- df_agbd[df_agbd$agbd > 0 & df_agbd$agbd < 100,]
df_agbd<- data.frame(df_agbd)
df_agbd <- df_agbd %>% group_by(gr=cut(df_agbd, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_agbd$gr)
num <- round(num/nrow(df_agbd)*100,2)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

p1<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=30), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 30))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(0, 30, 4),breaks =  seq(0, 30, 4))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="SAR-based AGBD (Field-ALS-SAR) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_AGBD_rf_Clip.tif"
df_agbd <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_agbd) <- c("agbd")
df_agbd <- df_agbd[df_agbd$agbd > 0 & df_agbd$agbd < 100,]
df_agbd<- data.frame(df_agbd)
df_agbd <- df_agbd %>% group_by(gr=cut(df_agbd, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_agbd$gr)
num <- round(num/nrow(df_agbd)*100,2)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

p2<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=30), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 30))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(0, 30, 4),breaks =  seq(0, 30, 4))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="SAR-based AGBD (Field-ALS-GEDI_SAS_RF-SAR) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_AGBD_glm_Clip.tif"
df_agbd <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_agbd) <- c("agbd")
df_agbd <- df_agbd[df_agbd$agbd > 0 & df_agbd$agbd < 100,]
df_agbd<- data.frame(df_agbd)
df_agbd <- df_agbd %>% group_by(gr=cut(df_agbd, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_agbd$gr)
num <- round(num/nrow(df_agbd)*100,2)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

p3<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=30), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 30))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(c)") + 
  scale_y_continuous(minor_breaks = seq(0, 30, 4),breaks =  seq(0, 30, 4))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="SAR-based AGBD (Field-ALS-GEDI_SAS_GLM-SAR) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(p1,p2,p3,nrow=3)
out = "E:\\ScanSAR\\ScanSAR\\Result_1004\\figure\\Bar_agbd_knp.jpg"
ggsave(out,height=20, width=20, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#make histogram of uncertainty -- abs and pct
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100_knp\\SAR_uncertainty_2018_1h.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(-0.01, 50, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,45,5)
l <- cbind(x,num)
df <- data.frame(l)

p1<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=50), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 50))+
  annotate("text",x=45,y=40,hjust = 0,size = 15, label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="Relative AGBD uncertainty (ALS) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))


dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100_knp\\SAR_uncertainty_2018_1h_rf.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(-0.01, 50, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,45,5)
l <- cbind(x,num)
df <- data.frame(l)

p2<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=50), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 50))+
  annotate("text",x=45,y=40,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="Relative AGBD uncertainty (GEDI_RF) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100_knp\\SAR_uncertainty_2018_1h_glm.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(-0.01, 50, by=5))) %>% arrange(as.numeric(gr))

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,45,5)
l <- cbind(x,num)
df <- data.frame(l)

p3<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=50), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 50))+
  annotate("text",x=45,y=40,hjust = 0,size = 15, label= "(c)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,50,5),digits = 1),
                     breaks = round(seq(0,50,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="Relative AGBD uncertainty (GEDI_GLM) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_uncertainty_als.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(-0.01, 6, by=1))) %>% arrange(as.numeric(gr))

mean(df_uncertainty$uncertainty)
sd(df_uncertainty$uncertainty)

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,5,1)
l <- cbind(x,num)
df <- data.frame(l)


p4<- ggplot(df, aes(x=x+0.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=80), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  annotate("text",x=5,y=60,hjust = 0,size = 15, label= "(d)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,6,1),digits = 1),
                     breaks = round(seq(0,6,1),digits = 1))+
  labs(y="Frequency (%)", 
       x="Absolute AGBD uncertainty (ALS) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))


dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_uncertainty_rf.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(-0.01, 6, by=1))) %>% arrange(as.numeric(gr))

mean(df_uncertainty$uncertainty)
sd(df_uncertainty$uncertainty)

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,5,1)
l <- cbind(x,num)
df <- data.frame(l)


p5<- ggplot(df, aes(x=x+0.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=80), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  annotate("text",x=5,y=60,hjust = 0,size = 15, label= "(e)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,6,1),digits = 1),
                     breaks = round(seq(0,6,1),digits = 1))+
  labs(y="Frequency (%)", 
       x="Absolute AGBD uncertainty (GEDI_RF) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))



dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\MC_100\\SAR_uncertainty_glm.tif"
df_uncertainty <- data.frame(rasterToPoints(raster(dir))[,3])
colnames(df_uncertainty) <- c("uncertainty")
df_uncertainty<- data.frame(df_uncertainty)
df_uncertainty <- df_uncertainty %>% group_by(gr=cut(uncertainty, breaks= seq(0, 6, by=1))) %>% arrange(as.numeric(gr))

mean(df_uncertainty$uncertainty)
sd(df_uncertainty$uncertainty)

num <- table(df_uncertainty$gr)
num <- round(num/nrow(df_uncertainty)*100,2)
x <- seq(0,5,1)
l <- cbind(x,num)
df <- data.frame(l)


p6<- ggplot(df, aes(x=x+0.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=80), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 80))+
  annotate("text",x=5,y=60,hjust = 0,size = 15, label= "(f)") + 
  scale_y_continuous(minor_breaks = seq(0,80,10),breaks =  seq(0,80,10))+
  scale_x_continuous(minor_breaks = round(seq(0,6,1),digits = 1),
                     breaks = round(seq(0,6,1),digits = 1))+
  labs(y="Frequency (%)", 
       x="Absolute AGBD uncertainty (GEDI_GLM) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=35))+
  theme(plot.title = element_text(hjust = 0.5))



ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2)
out = "E:\\ScanSAR\\ScanSAR\\Result_1004\\figure\\Bar_uncertainty_knp.jpg"
ggsave(out,height=20, width=40, dpi=600)





#histogram of ALS, GEDI SAS RF and GLM
#als
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_MC\\ALS.csv"
df <- read.csv(dir)
colnames(df) <- c("ID","AGBD")
df<- data.frame(df)
df <- df %>% group_by(gr=cut(AGBD, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))
num <- table(df$gr)
num <- round(num/nrow(df)*100,1)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

p1<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=25), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 25))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(a)") + 
  scale_y_continuous(minor_breaks = seq(0, 100, 5),breaks =  seq(0, 100, 5))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="ALS AGBD (MCH*CC) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))


#gedi -rf
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_MC\\GEDI_RF.csv"
df <- read.csv(dir)
colnames(df) <- c("ID","AGBD")
df<- data.frame(df)
df <- df %>% group_by(gr=cut(AGBD, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))
num <- table(df$gr)
num <- round(num/nrow(df)*100,1)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

p2<- ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=25), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 25))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(0, 100, 5),breaks =  seq(0, 100, 5))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="GEDI AGBD (GLM) (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))


#gedi -glm
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_MC\\GEDI_GLM.csv"
df <- read.csv(dir)
colnames(df) <- c("ID","AGBD")
df<- data.frame(df)
df <- df %>% group_by(gr=cut(AGBD, breaks= seq(0, 100, by=5))) %>% arrange(as.numeric(gr))
num <- table(df$gr)
num <- round(num/nrow(df)*100,1)
x <- seq(0,95,5)
l <- cbind(x,num)
df <- data.frame(l)

ggplot(df, aes(x=x+2.5, y=num,group = x)) + 
  geom_bar(stat='identity')+
  geom_text(aes(label = num,y=25), size = 8)+
  theme_bw()+
  coord_cartesian(ylim = c(0, 25))+
  annotate("text",x=90,y=15,hjust = 0,size = 15, label= "(b)") + 
  scale_y_continuous(minor_breaks = seq(0, 100, 5),breaks =  seq(0, 100, 5))+
  scale_x_continuous(minor_breaks = round(seq(0,100,5),digits = 1),
                     breaks = round(seq(0,100,5),digits = 1))+
  labs(y="Frequency (%)", 
       x="SAR-based AGBD uncertainty (Field-ALS-GEDI_SAS_GLM-SAR) (%)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))



ggarrange(p1,p2,nrow=2)
out = "E:\\ScanSAR\\ScanSAR\\Result_1004\\figure\\Train.jpg"
ggsave(out,height=15, width=20, dpi=600)





# ------------------------------------------------------------------------------------------------ #
#Field-ALS-SAR AGBD vs. ALS AGBD
list_ini <- 0
list_ini2 <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1115\\table"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")
#list_dir <- c("Welverdiendt")

dir_hcc_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_25m.csv"
data_hcc_sar = read.csv(dir_hcc_sar)

for (l in list_dir){
  dir_chm <- file.path(dir,paste0(l,"_val_H.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))

  dir_cover <- file.path(dir,paste0(l,"_val_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))

  dir_AGBD <- file.path(dir,paste0(l,"_val_FAS_AGBD.csv"))
  Data = read.csv(dir_AGBD)
  Data_AGBD = subset(Data, select=c(2,5))

  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  Data_index <- merge(Data_index, Data_AGBD, by.x = "FID", by.y = "FID")
  colnames(Data_index) <- c("FID","CHM","Cover","AGBD")
  Data_index <- na.omit(Data_index)
  
  Data_index$Cover <- Data_index$Cover/10000
  Data_index$ALS_AGBD <- 18.839*Data_index$CHM*Data_index$Cover - 1.213  

  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
}

Data <- Data_index
R2 <- round(cor(Data$AGBD, Data$ALS_AGBD)^2,2)
bias <- mean(Data$AGBD - Data$ALS_AGBD)
RMSE <- sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))
rbias <- 100*mean(Data$AGBD - Data$ALS_AGBD)/mean(Data$ALS_AGBD)
rRMSE <- 100*sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))/mean(Data$ALS_AGBD)
mean(Data_index$AGBD)
p1<- ggplot(Data_index, aes(x=ALS_AGBD, y=AGBD))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=1,y=90,hjust = 0,size = 9,
           label= paste(expression(Linear~model~R^2),": ",
                        round(R2,3)), parse=TRUE) + 
  annotate("text",x=0,y=80,hjust = 0,size = 9,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,2)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,2),"%",
             "\n" , " Bias: ", round(bias,2)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,2)," %")) + 
  annotate("text",x=75,y=95,hjust = 0,size = 15,label= "(a)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(y="Three-stage AGBD (Mg/ha)", 
       x="ALS-based (MCH*CC) AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


Data_index <- Data_index[Data_index$ALS_AGBD < 60,]
breakbin = seq(-40,40,5)
Data_index$bias <- Data_index$AGBD - Data_index$ALS_AGBD
Data_index$group_RH <- cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),dig.lab=1)
Data_index <- na.omit(Data_index)
bi1<- ggplot(Data_index, aes(x=group_RH, y=bias, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.1) +
  geom_boxplot()+
  coord_cartesian(ylim = c(-35, 25))+
  scale_y_continuous(minor_breaks = breakbin,breaks = breakbin)+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="Bias of three-stage AGBD (Mg/ha)")+
  annotate("text",x=10,y=20,hjust = 0,size = 15,label= "(a)") + 
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))

a1 <- tapply(Data_index$bias, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                                  dig.lab=1), mean)
b1 <- tapply(Data_index$ALS_AGBD, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                                      dig.lab=1), mean)
RB1 <- a1/b1*100

a1
b1
RB1

x1 <- seq(0,55,5)
l1 <- cbind(x1,RB1)
df1 <- data.frame(l1)

r1 <- ggplot(df1, aes(x=x1+2.5, y=RB1)) + 
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="%Bias of three-stage AGBD (%)")+
  annotate("text",x=50,y=90,hjust = 0,size = 15,label= "(b)") + 
  coord_cartesian(ylim = c(-50, 100))+
  scale_y_continuous(minor_breaks = seq(-200, 100, 10),breaks = seq(-200, 100, 10))+
  scale_x_continuous(minor_breaks = seq(0,60,5),breaks = seq(0,60,5))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))

# ------------------------------------------------------------------------------------------------ #
#Field-ALS-GEDI_RF-SAR AGBD vs. ALS AGBD
list_ini <- 0
list_ini2 <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1115\\table"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")
#list_dir <- c("Welverdiendt")

dir_hcc_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_25m.csv"
data_hcc_sar = read.csv(dir_hcc_sar)

for (l in list_dir){
  dir_chm <- file.path(dir,paste0(l,"_val_H.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_val_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  
  dir_AGBD <- file.path(dir,paste0(l,"_val_FAS_RF_AGBD.csv"))
  Data = read.csv(dir_AGBD)
  Data_AGBD = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  Data_index <- merge(Data_index, Data_AGBD, by.x = "FID", by.y = "FID")
  colnames(Data_index) <- c("FID","CHM","Cover","AGBD")
  Data_index <- na.omit(Data_index)
  
  Data_index$Cover <- Data_index$Cover/10000
  Data_index$ALS_AGBD <- 18.839*Data_index$CHM*Data_index$Cover - 1.213  
  
  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
  
}
Data <- Data_index
R2 <- round(cor(Data$AGBD, Data$ALS_AGBD)^2,2)
bias <- mean(Data$AGBD - Data$ALS_AGBD)
RMSE <- sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))
rbias <- 100*mean(Data$AGBD - Data$ALS_AGBD)/mean(Data$ALS_AGBD)
rRMSE <- 100*sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))/mean(Data$ALS_AGBD)
mean(Data_index$AGBD)
p2<- ggplot(Data_index, aes(x=ALS_AGBD, y=AGBD))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=1,y=90,hjust = 0,size = 9,
           label= paste(expression(Linear~model~R^2),": ",
                        round(R2,3)), parse=TRUE) + 
  annotate("text",x=0,y=80,hjust = 0,size = 9,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,2)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,2),"%",
             "\n" , " Bias: ", round(bias,2)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,2)," %")) + 
  annotate("text",x=75,y=95,hjust = 0,size = 15,label= "(b)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(y="Four-stage (RF) AGBD (Mg/ha)", 
       x="ALS-based (MCH*CC) AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


Data_index <- Data_index[Data_index$ALS_AGBD < 60,]
breakbin = seq(-40,40,5)
Data_index$bias <- Data_index$AGBD - Data_index$ALS_AGBD
Data_index$group_RH <- cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),dig.lab=1)
Data_index <- na.omit(Data_index)
bi2<- ggplot(Data_index, aes(x=group_RH, y=bias, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.1) +
  geom_boxplot()+
  coord_cartesian(ylim = c(-35, 25))+
  scale_y_continuous(minor_breaks = breakbin,breaks = breakbin)+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="Bias of four-stage (RF) AGBD (Mg/ha)")+
  annotate("text",x=10,y=20,hjust = 0,size = 15,label= "(c)") + 
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


a2 <- tapply(Data_index$bias, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                                  dig.lab=1), mean)
b2 <- tapply(Data_index$ALS_AGBD, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                                      dig.lab=1), mean)
RB2 <- a2/b2*100

a2
b2
RB2

x2 <- seq(0,55,5)
l2 <- cbind(x2,RB2)
df2 <- data.frame(l2)

r2 <- ggplot(df2, aes(x=x2+2.5, y=RB2)) + 
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="%Bias of four-stage (RF) AGBD (%)")+
  annotate("text",x=50,y=90,hjust = 0,size = 15,label= "(d)") + 
  coord_cartesian(ylim = c(-50, 100))+
  scale_y_continuous(minor_breaks = seq(-200, 100, 10),breaks = seq(-200, 100, 10))+
  scale_x_continuous(minor_breaks = seq(0,60,5),breaks = seq(0,60,5))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


# ------------------------------------------------------------------------------------------------ #
#Field-ALS-GEDI_GLM-SAR AGBD vs. ALS AGBD
list_ini <- 0
list_ini2 <- 0

dir <- "E:\\ScanSAR\\ScanSAR\\Result_1115\\table"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")

dir_hcc_sar <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\result\\hcc_sar_25m.csv"
data_hcc_sar = read.csv(dir_hcc_sar)

for (l in list_dir){
  dir_chm <- file.path(dir,paste0(l,"_val_H.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_val_CC.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  
  dir_AGBD <- file.path(dir,paste0(l,"_val_FAS_GLM_AGBD.csv"))
  Data = read.csv(dir_AGBD)
  Data_AGBD = subset(Data, select=c(2,5))
  
  #concat chm, cc, sar to index dataframe
  Data_index <- merge(Data_chm, Data_cover, by.x = "FID", by.y = "FID")
  Data_index <- merge(Data_index, Data_AGBD, by.x = "FID", by.y = "FID")
  colnames(Data_index) <- c("FID","CHM","Cover","AGBD")
  Data_index <- na.omit(Data_index)
  
  Data_index$Cover <- Data_index$Cover/10000
  Data_index$ALS_AGBD <- 18.839*Data_index$CHM*Data_index$Cover - 1.213  
  
  #data clean
  Data_index <- Data_index[is.finite(rowSums(Data_index)),]
  Data_index <- Data_index[Data_index$CHM > 0,]
  Data_index <- Data_index[Data_index$Cover > 0,]
  
}
Data <- Data_index
R2 <- round(cor(Data$AGBD, Data$ALS_AGBD)^2,2)
bias <- mean(Data$AGBD - Data$ALS_AGBD)
RMSE <- sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))
rbias <- 100*mean(Data$AGBD - Data$ALS_AGBD)/mean(Data$ALS_AGBD)
rRMSE <- 100*sqrt(mean((Data$AGBD - Data$ALS_AGBD)^2))/mean(Data$ALS_AGBD)

p3<- ggplot(Data_index, aes(x=ALS_AGBD, y=AGBD))+ 
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(0, 100))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=1,y=90,hjust = 0,size = 9,
           label= paste(expression(Linear~model~R^2),": ",
                        round(R2,3)), parse=TRUE) + 
  annotate("text",x=0,y=80,hjust = 0,size = 9,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,2)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,2),"%",
             "\n" , " Bias: ", round(bias,2)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,2)," %")) + 
  annotate("text",x=75,y=95,hjust = 0,size = 15,label= "(f)") + 
  geom_abline(intercept = 0, slope = 1,color="black", linetype="solid", size=1.5)+
  labs(y="Four-stage (GLM) AGBD (Mg/ha)", 
       x="ALS-based (MCH*CC) AGBD (Mg/ha)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))


Data_index <- Data_index[Data_index$ALS_AGBD < 60,]
breakbin = seq(-40,40,5)
Data_index$bias <- Data_index$AGBD - Data_index$ALS_AGBD
Data_index$group_RH <- cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),dig.lab=1)
table(Data_index$group_RH)/nrow(Data_index)*100
Data_index <- na.omit(Data_index)
bi3<- ggplot(Data_index, aes(x=group_RH, y=bias, group=group_RH)) + 
  stat_boxplot(geom ='errorbar', width = 0.1) +
  geom_boxplot()+
  coord_cartesian(ylim = c(-35, 25))+
  scale_y_continuous(minor_breaks = breakbin,breaks = breakbin)+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="Bias of four-stage (GLM) AGBD (Mg/ha)")+
  annotate("text",x=10,y=20,hjust = 0,size = 15,label= "(e)") + 
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))

a3 <- tapply(Data_index$bias, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                            dig.lab=1), mean)
b3 <- tapply(Data_index$ALS_AGBD, cut(Data_index$ALS_AGBD,breaks = seq(0,60,5),
                                      dig.lab=1), mean)
RB3 <- a3/b3*100

a3
b3
RB3

x3 <- seq(0,55,5)
l3 <- cbind(x3,RB3)
df3 <- data.frame(l3)

r3 <- ggplot(df3, aes(x=x3+2.5, y=RB3)) + 
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x="ALS-based (MCH*CC) AGBD (Mg/ha)",
       y="%Bias of four-stage (GLM) AGBD (%)")+
  annotate("text",x=50,y=90,hjust = 0,size = 15,label= "(f)") + 
  coord_cartesian(ylim = c(-50, 100))+
  scale_y_continuous(minor_breaks = seq(-200, 100, 10),breaks = seq(-200, 100, 10))+
  scale_x_continuous(minor_breaks = seq(0,60,5),breaks = seq(0,60,5))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(text=element_text(family="A"))+
  theme(plot.title = element_text(hjust = 0.5))


round(table(Data_index$group_RH)/nrow(Data_index)*100,2)

ggarrange(p1,p2,p3,ncol=3)
out = "E:\\ScanSAR\\ScanSAR\\Result_1115\\scatterplot_val.jpg"
ggsave(out,height=12, width=36, dpi=600)


ggarrange(bi1,bi2,bi3,ncol=3)
out = "E:\\ScanSAR\\ScanSAR\\Result_1115\\biasplot_val.jpg"
ggsave(out,height=12, width=45, dpi=600)

ggarrange(r1,r2,r3,nrow=3)
out = "E:\\ScanSAR\\ScanSAR\\Result_1115\\biasplot_relative_val.jpg"
ggsave(out,height=30, width=15, dpi=600)



ggarrange(bi1,r1,bi2,r2,bi3,r3,nrow=3,ncol=2)
out = "E:\\ScanSAR\\ScanSAR\\Result_1115\\biasplot_6.jpg"
ggsave(out,height=30, width=30, dpi=600)




#25m 
# ------------------------------------------------------------------------------------------------ #
#update 01/11/2024
list_ini <- 0
dir <- "E:\\ScanSAR\\ScanSAR\\Result_20140111\\ALS_25m"
dir_sar <- "E:\\ScanSAR\\ScanSAR\\Result_20140111\\SAR_25m"
list_dir <- c("Agincourt", "Ireagh", "Justicia", "Welverdiendt", 
              "Limpopo1", "Limpopo2", "Limpopo3")

for (l in list_dir){
  print(l)
  dir_chm <- file.path(dir,paste0(l,"_H_25m.csv"))
  Data = read.csv(dir_chm)
  Data_chm = subset(Data, select=c(2,5))
  
  dir_cover <- file.path(dir,paste0(l,"_CC_25m.csv"))
  Data = read.csv(dir_cover)
  Data_cover = subset(Data, select=c(2,5))
  
  dir_sar2 <- file.path(dir_sar,paste0(l,"_SAR_25m.csv"))
  Data = read.csv(dir_sar2)
  Data_sar = subset(Data, select=c(2,5))
  
  
  #concat chm, cc, sar to index dataframe
  Data_chm <- merge(Data_chm, Data_cover, by.x = "TARGET_FID", by.y = "TARGET_FID")
  Data_chm <- merge(Data_chm, Data_sar, by.x = "TARGET_FID", by.y = "TARGET_FID")
  colnames(Data_chm) <- c("TARGET_FID","CHM","Cover","SAR")
  
  #data conversion, 1ha agbd cal
  Data_chm$Cover <- Data_chm$Cover/625
  Data_chm$AGBD <- 9.0665*Data_chm$CHM*Data_chm$Cover
  
  #merge by site
  list_ini <- rbind(list_ini, Data_chm)
}

Data <- list_ini[list_ini$AGBD > 1,]

Data <- Data %>% group_by(gr=cut(AGBD, breaks= seq(0, 100, by=1))) %>% arrange(as.numeric(gr))
Data <- Data %>% group_by(gr) %>% slice_sample(n=500)

#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR))
bias<-mean(pred - Data$AGBD)
RMSE<-sqrt(mean((pred - Data$AGBD)^2))
rbias<-100*mean(pred - Data$AGBD)/mean(Data$AGBD)
rRMSE<-100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)
Data$preds <- pred

summary(Data$AGBD)

p1 <- ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=61, y=-20,hjust = 0,size = 7,
           label= paste(expression(Logarithmic~model~R^2),": ",
                        round(summary(reg3)$adj.r.squared,3)), parse=TRUE) + 
  annotate("text",x=60,y=-21.5,hjust = 0,size = 7,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(bias,3)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,3)," %")) + 
  labs(x="25m area-based (H*CC) ALS AGBD (Mg/ha)", 
       y="Mean of 25m SAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))



#1ha
#sar
sar_date <- "ALS_1ha_2018.csv"
#sar_date <- "Finebeam_1ha.csv"
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats"
dir_sar <- file.path(dir,sar_date)
Data = read.csv(dir_sar)
Data = subset(Data, select=c("FID","MEAN"))
colnames(Data) <- c("FID","SAR")

#als
dir <- "E:\\ScanSAR\\ScanSAR\\Result_1004\\table_stats\\ALS_1ha_AGBD_HCC.csv"
Data_als = read.csv(dir)

Data <- merge(Data, Data_als, by.x = "FID", by.y = "FID")

Data <- Data[Data$AGBD > 0,]
Data <- Data[Data$Cover > 0,]
Data <- Data[Data$SAR > -25 & Data$SAR < -10,]

#AGBD
reg3 <- lm(log(AGBD) ~ SAR, data = Data)
pred <- exp((coef(reg3)[1] + coef(reg3)[2]*Data$SAR)) 

bias<-mean(pred - Data$AGBD)
RMSE<-sqrt(mean((pred - Data$AGBD)^2))
rbias<-100*mean(pred - Data$AGBD)/mean(Data$AGBD)
rRMSE<-100*sqrt(mean((pred - Data$AGBD)^2))/mean(Data$AGBD)
Data$preds <- pred

p2 <- ggplot(Data, aes(x=AGBD, y=SAR))+ 
  geom_pointdensity()+
  theme_bw()+
  geom_smooth(
    method="lm",
    formula = 'y ~ log(x)', 
    se= F, 
    size = 1.5, 
    linetype = "solid",
    colour = "red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 100),ylim = c(-25, -10))+
  scale_color_viridis(direction = 1)+
  annotate("text", x=61, y=-20,hjust = 0,size = 7,
           label= paste(expression(Logarithmic~model~R^2),": ",
                        round(summary(reg3)$adj.r.squared,3)), parse=TRUE) + 
  annotate("text",x=60,y=-21.5,hjust = 0,size = 7,
           label= paste(
             "\n" , " RMSE: ", round(RMSE,3)," Mg/ha",
             "\n" , " rRMSE: ", round(rRMSE,3),"%",
             "\n" , " Bias: ", round(bias,3)," Mg/ha",
             "\n" , " %Bias: ", round(rbias,3)," %")) + 
  labs(x="1ha area-based (H*CC) ALS AGBD (Mg/ha)", 
       y="Mean of 1ha SAR backscatter (2018 Sep) (dB)")+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(text=element_text(size=30))+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p1,p2)
out = "E:\\ScanSAR\\ScanSAR\\Result_20140111\\Scantterplot_SAR_AGBD.jpg"
ggsave(out,height=12, width=24, dpi=600)






