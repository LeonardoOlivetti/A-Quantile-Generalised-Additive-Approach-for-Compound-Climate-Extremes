#############START############

#setwd("..")

#Source

source("functions.R")
source("functions_EVT.R")

#Packages
library(tidyverse) # package for plotting
library(beepr) # beeeeep
library(tseries)
library(qgam)
library(mgcv)
library(forecast)
library(ggpubr)
library(purrr)
library(timetk)
library(gratia)
library(LSTS)
library(zoo)
library(lubridate)
library(tsibble)
library(xtable)
library(itsadug)
library(extRemes)
library(evgam)
library(quantreg)
library(sf)
library(scales)
library(WRTDStidal)

#Notes

#https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html


#load_data
dat<-readRDS("data/X_var.Rda")%>%
  mutate(date=as.character(date))

out<-readRDS("data/Y_var.Rda")

##CLEANING AND STANDARIDSATION##

tot<-left_join(out,dat)%>%
  mutate(fac_lat=as.character(lat),fac_lon=as.character(lon),fac_lat_jet=as.character(lat_jet))%>%
  mutate(date=as_date(date),jet_proximity=abs(lat-lat_jet),lon=ifelse(lon<5,lon+360,lon))%>%
  mutate(year=year(date))%>%
  as_tsibble(index=date,key=c(lat,lon))%>%
  fill_gaps()%>%
  tk_augment_lags(c(jet_strength,NAO,US_temp,lat_jet,jet_proximity,wind,prec), .lags = c(1:2), .names = "auto")%>%
  tk_augment_leads(c(jet_strength,NAO,US_temp,lat_jet,jet_proximity,wind,prec), .lags = c(-1:-2), .names = "auto")%>%
  mutate(time_cont=as.integer(as.Date(date, "%Y-%m-%d")-as.Date("1959-01-01")))%>%
  filter(month(date)<3 | month(date)>10)%>%
  mutate(across(c(where(is.numeric)), scale))%>%
  mutate(month=month(date))%>%
  mutate(across(starts_with("fac"),as.numeric))%>%
  mutate(across(starts_with("fac"),round,digits=2))%>%
  filter(date>"1959-04-01"& date<"2022-08-31")%>%
  group_by(year,month)%>%
  mutate(season_year=ceiling(cur_group_id()/4))%>%
  ungroup()%>%
  mutate(month=as.factor(month(date)))%>%
  mutate(date=as.character(date))

gc()

#Train-test split


set.seed(23)
year_samp<-tot%>%
  distinct(season_year)%>%
  sample_frac(0.5)%>%
  pull()


tot_train<-tot%>%
  filter(season_year%in%year_samp)%>%
  group_by(fac_lon,fac_lat)%>%
  mutate(qu_95_wind=quantile(wind,0.95),
         qu_99_wind=quantile(wind,0.99),
         qu_95_prec=quantile(prec,0.95),
         qu_99_prec=quantile(prec,0.99)
         )%>%
  ungroup()

rm(tot,dat,out)

gc()



########MODELS############

####### Wind ###########

q_base_95_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_95_wind+s(wind_lag1,k=5),data=tot_train,qu=0.95)

saveRDS(q_base_95_wind,file="q_base_95_wind.Rda")

rm(q_base_95_wind)

gc()

q_cold_95_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_95_wind+s(wind_lag1,k=5),data=tot_train,qu=0.95)

saveRDS(q_cold_95_wind,file="q_cold_95_wind.Rda")

rm(q_cold_95_wind)

gc()

q_jet_95_wind<-qgam(wind~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_95_wind+s(wind_lag1,k=5),data=tot_train,qu=0.95)


saveRDS(q_jet_95_wind,file="q_jet_95_wind.Rda")

rm(q_jet_95_wind)

gc()

q_base_99_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_99_wind+s(wind_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_base_99_wind,file="q_base_99_wind.Rda")

rm(q_base_99_wind)

gc()


q_cold_99_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_99_wind+s(wind_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_cold_99_wind,file="q_cold_99_wind.Rda")

rm(q_cold_99_wind)

gc()

q_jet_99_wind<-qgam(wind~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_99_wind+s(wind_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_jet_99_wind,file="q_jet_99_wind.Rda")

rm(q_jet_99_wind)


gc() 


###############PREC##############


q_base_95_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_95_prec+s(prec_lag1,k=5),data=tot_train,qu=0.95)

saveRDS(q_base_95_prec,file="q_base_95_prec.Rda")

rm(q_base_95_prec)

gc()


q_cold_95_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_95_prec+s(prec_lag1,k=5),data=tot_train,qu=0.95)

saveRDS(q_cold_95_prec,file="q_cold_95_prec.Rda")

rm(q_cold_95_prec)

gc()

q_jet_95_prec<-qgam(prec~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_95_prec+s(prec_lag1,k=5),data=tot_train,qu=0.95)

saveRDS(q_jet_95_prec,file="q_jet_95_prec.Rda")

rm(q_jet_95_prec)

gc()


q_base_99_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_99_prec+s(prec_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_base_99_prec,file="q_base_99_prec.Rda")

rm(q_base_99_prec)

gc()

q_cold_99_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_99_prec+s(prec_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_cold_99_prec,file="q_cold_99_prec.Rda")

rm(q_cold_99_prec)

gc()

q_jet_99_prec<-qgam(prec~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_99_prec+s(prec_lag1,k=5),data=tot_train,qu=0.99)

saveRDS(q_jet_99_prec,file="q_jet_99_prec.Rda")

rm(q_jet_99_prec)

gc()
