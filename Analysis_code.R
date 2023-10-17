#Notes
#lat - latitude of the point in Western Europe
#lon - longitude of the point in Western Europe
#prec - daily total precipitation, anomaly in mm at the surface
#wind - daily mean wind speed anomaly at 10m height, meters per second
#temp - daily mean temperature anomaly at 2m height, kelvin 

#Source

source("functions.R")

#Packages
library(tidyverse)
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
library(ggpattern)
library(metR)

#Notes


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
  filter(season_year%in%year_samp)
  
tot_left<-setdiff(tot,tot_train)

year_samp_test<-tot_left%>%
  distinct(season_year)%>%
  sample_frac(0.5)%>%
  pull()

tot_val<-tot_left%>%
  filter(season_year%in%year_samp_test)

tot_test<-setdiff(tot_left,tot_val)

tot_train<-tot_train%>%
  group_by(fac_lon,fac_lat)%>%
  mutate(qu_95_wind=quantile(wind,0.95),
         qu_99_wind=quantile(wind,0.99),
         qu_95_prec=quantile(prec,0.95),
         qu_99_prec=quantile(prec,0.99),
  )%>%
  ungroup()


clim<-tot_train%>%
  select(fac_lat,fac_lon,qu_95_wind,qu_99_wind,qu_95_prec,qu_99_prec)%>%
  distinct(fac_lat,fac_lon,.keep_all = TRUE)
  
tot_val<-tot_val%>%
  left_join(clim,by=c("fac_lat","fac_lon"))

tot_test<-tot_test%>%
  left_join(clim,by=c("fac_lat","fac_lon"))

#########QGAMS###########

#Wind

q_base_95_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_95_wind,data=tot_train,qu=0.95)

q_cold_95_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_95_wind,data=tot_train,qu=0.95)

q_jet_95_wind<-qgam(wind~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_95_wind,data=tot_train,qu=0.95)

q_base_99_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_99_wind,data=tot_train,qu=0.99)

q_cold_99_wind<-qgam(wind~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_99_wind,data=tot_train,qu=0.99)

q_jet_99_wind<-qgam(wind~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_99_wind,data=tot_train,qu=0.99)

#Prec


q_base_95_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_95_prec,data=tot_train,qu=0.95)


q_cold_95_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_95_prec,data=tot_train,qu=0.95)


q_jet_95_prec<-qgam(prec~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_95_prec,data=tot_train,qu=0.95)


q_base_99_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+month+qu_99_prec,data=tot_train,qu=0.99)


q_cold_99_prec<-qgam(prec~
                       s(lat,lon,k=20)+
                       s(time_cont,k=5)+
                       s(US_temp_lag2,k=5)+month+qu_99_prec,data=tot_train,qu=0.99)


q_jet_99_prec<-qgam(prec~
                      s(lat,lon,k=20)+
                      s(time_cont,k=5)+s(US_temp_lag2,k=5)+
                      s(jet_strength_lag1)+s(lat_jet_lag1)+
                      s(jet_proximity_lag1)+
                      s(NAO_lag1,k=5)+month+qu_99_prec,data=tot_train,qu=0.99)


#PREDICTED VALUES QGAMS

pred_q_base_wind_99<-predict(q_base_wind_99,tot_test)
pred_q_base_wind_95<-predict(q_base_wind_95,tot_test)
pred_q_cold_wind_99<-predict(q_cold_wind_99,tot_test)
pred_q_cold_wind_95<-predict(q_cold_wind_95,tot_test)
pred_q_jet_wind_99<-predict(q_jet_wind_99,tot_test)
pred_q_jet_wind_95<-predict(q_jet_wind_95,tot_test)

pred_q_base_prec_99<-predict(q_base_prec_99,tot_test)
pred_q_base_prec_95<-predict(q_base_prec_95,tot_test)
pred_q_cold_prec_99<-predict(q_cold_prec_99,tot_test)
pred_q_cold_prec_95<-predict(q_cold_prec_95,tot_test)
pred_q_jet_prec_99<-predict(q_jet_prec_99,tot_test)
pred_q_jet_prec_95<-predict(q_jet_prec_95,tot_test)


######################Bias analysis(Fiugre 4-5) and comparison to seasonal climatology (Figure 6-7) ##############

#Wind

p_q_base_wind_99<-Quant_eval(pred_q_base_wind_99,wind,0.99,tot_test)
p_q_base_99<-overpred(p_q_base_wind_99,0.99,direc=-1,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_base_wind_95<-Quant_eval(pred_q_base_wind_95,wind,0.95,tot_test)
p_q_base_95<-overpred(p_q_base_wind_95,0.95,direc=-1,low_lim=0,up_lim=0.1,midpo=0.05)
p_q_cold_wind_99<-Quant_eval(pred_q_cold_wind_99,wind,0.99,tot_test)
p_q_cold_99<-overpred(p_q_cold_wind_99,0.99,direc=-1,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_cold_wind_95<-Quant_eval(pred_q_cold_wind_95,wind,0.95,tot_test)
p_q_cold_95<-overpred(p_q_cold_wind_95,0.95,direc=-1,low_lim=0,up_lim=0.1,midpo=0.05)
p_q_jet_wind_99<-Quant_eval(pred_q_jet_wind_99,wind,0.99,tot_test)
p_q_jet_99<-overpred(p_q_jet_wind_99,0.99,direc=-1,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_jet_wind_95<-Quant_eval(pred_q_jet_wind_95,wind,0.95,tot_test)
p_q_jet_95<-overpred(p_q_jet_wind_95,0.95,direc=-1,low_lim=0,up_lim=0.1,midpo=0.05)

a<-annotate_figure(ggarrange(p_q_base_95,
                             p_q_cold_95,
                             p_q_jet_95,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))


b<-annotate_figure(ggarrange(p_q_base_99,
                             p_q_cold_99,
                             p_q_jet_99,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM bias: Daily 10m Wind Speed",face="bold", size=15)
)

R2_q_base_wind_99<-Pseudo_R2_point(pred_q_base_wind_99,wind,0.99,tot_train,tot_test)
p_q_base_99_R2<-Plot_R2(R2_q_base_wind_99,0.99,up_lim = 0.21,low_lim=-0.21)

R2_q_base_wind_95<-Pseudo_R2_point(pred_q_base_wind_95,wind,0.95,tot_train,tot_test)
p_q_base_95_R2<-Plot_R2(R2_q_base_wind_95,0.95,up_lim = 0.21,low_lim=-0.21)

R2_q_cold_wind_99<-Pseudo_R2_point(pred_q_cold_wind_99,wind,0.99,tot_train,tot_test)
p_q_cold_99_R2<-Plot_R2(R2_q_cold_wind_99,0.99,up_lim = 0.21,low_lim=-0.21)

R2_q_cold_wind_95<-Pseudo_R2_point(pred_q_cold_wind_95,wind,0.95,tot_train,tot_test)
p_q_cold_95_R2<-Plot_R2(R2_q_cold_wind_95,0.95,up_lim = 0.21,low_lim=-0.21)

R2_q_jet_wind_99<-Pseudo_R2_point(pred_q_jet_wind_99,wind,0.99,tot_train,tot_test)
p_q_jet_99_R2<-Plot_R2(R2_q_jet_wind_99,0.99,up_lim = 0.21,low_lim = -0.21)

R2_q_jet_wind_95<-Pseudo_R2_point(pred_q_jet_wind_95,wind,0.95,tot_train,tot_test)
p_q_jet_95_R2<-Plot_R2(R2_q_jet_wind_95,0.95,up_lim = 0.21,low_lim=-0.21)


a<-annotate_figure(ggarrange(p_q_base_95_R2,
                             p_q_cold_95_R2,
                             p_q_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))


b<-annotate_figure(ggarrange(p_q_base_99_R2,
                             p_q_cold_99_R2,
                             p_q_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM performance: Daily 10m Wind Speed",face="bold", size=15)
)

#overview

Pseudo_R2_jet_99_QGAM_wind<-Pseudo_R2_tot(pred_q_jet_wind_99,wind,0.99,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_wind<-Pseudo_R2_tot(pred_q_jet_wind_95,wind,0.95,tot_train,tot_test)

Pseudo_R2_jet_99_QGAM_wind
Pseudo_R2_jet_95_QGAM_wind


#Prec

p_q_base_prec_99<-Quant_eval(pred_q_base_prec_99,prec,0.99,tot_test)
p_q_base_99<-overpred(p_q_base_prec_99,0.99,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_base_prec_95<-Quant_eval(pred_q_base_prec_95,prec,0.95,tot_test)
p_q_base_95<-overpred(p_q_base_prec_95,0.95,low_lim=0,up_lim=0.1,midpo=0.05)
p_q_cold_prec_99<-Quant_eval(pred_q_cold_prec_99,prec,0.99,tot_test)
p_q_cold_99<-overpred(p_q_cold_prec_99,0.99,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_cold_prec_95<-Quant_eval(pred_q_cold_prec_95,prec,0.95,tot_test)
p_q_cold_95<-overpred(p_q_cold_prec_95,0.95,low_lim=0,up_lim=0.1,midpo=0.05)
p_q_jet_prec_99<-Quant_eval(pred_q_jet_prec_99,prec,0.99,tot_test)
p_q_jet_99<-overpred(p_q_jet_prec_99,0.99,low_lim=0,up_lim=0.02,midpo=0.01)
p_q_jet_prec_95<-Quant_eval(pred_q_jet_prec_95,prec,0.95,tot_test)
p_q_jet_95<-overpred(p_q_jet_prec_95,0.95,low_lim=0,up_lim=0.1,midpo=0.05)

a<-annotate_figure(ggarrange(p_q_base_95,
                             p_q_cold_95,
                             p_q_jet_95,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_q_base_99,
                             p_q_cold_99,
                             p_q_jet_99,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM bias: Daily Precipitation",face="bold", size=15)
)

R2_q_base_prec_99<-Pseudo_R2_point(pred_q_base_prec_99,prec,0.99,tot_train,tot_test)
p_q_base_99_R2<-Plot_R2(R2_q_base_prec_99,0.99,up_lim = 0.21,low_lim=-0.21,opt="G")

R2_q_base_prec_95<-Pseudo_R2_point(pred_q_base_prec_95,prec,0.95,tot_train,tot_test)
p_q_base_95_R2<-Plot_R2(R2_q_base_prec_95,0.95,up_lim = 0.21,low_lim=-0.21,opt="G")

R2_q_cold_prec_99<-Pseudo_R2_point(pred_q_cold_prec_99,prec,0.99,tot_train,tot_test)
p_q_cold_99_R2<-Plot_R2(R2_q_cold_prec_99,0.99,up_lim = 0.21,low_lim=-0.21,opt="G")

R2_q_cold_prec_95<-Pseudo_R2_point(pred_q_cold_prec_95,prec,0.95,tot_train,tot_test)
p_q_cold_95_R2<-Plot_R2(R2_q_cold_prec_95,0.95,up_lim = 0.21,low_lim=-0.21,opt="G")

R2_q_jet_prec_99<-Pseudo_R2_point(pred_q_jet_prec_99,prec,0.99,tot_train,tot_test)
p_q_jet_99_R2<-Plot_R2(R2_q_jet_prec_99,0.99,up_lim = 0.21,low_lim=-0.21,opt="G")

R2_q_jet_prec_95<-Pseudo_R2_point(pred_q_jet_prec_95,prec,0.95,tot_train,tot_test)
p_q_jet_95_R2<-Plot_R2(R2_q_jet_prec_95,0.95,up_lim = 0.21,low_lim=-0.21,opt="G")



a<-annotate_figure(ggarrange(p_q_base_95_R2,
                             p_q_cold_95_R2,
                             p_q_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_q_base_99_R2,
                             p_q_cold_99_R2,
                             p_q_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM performance: Daily Precipitation",face="bold", size=15)
)

Pseudo_R2_jet_99_QGAM_prec<-Pseudo_R2_tot(pred_q_jet_prec_99,prec,0.99,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_prec<-Pseudo_R2_tot(pred_q_jet_prec_95,prec,0.95,tot_train,tot_test)

Pseudo_R2_jet_99_QGAM_prec
Pseudo_R2_jet_95_QGAM_prec


############Figure 3 - partial effects##############

a<-annotate_figure(ggarrange(
  annotate_figure(
    ggarrange(
      plot_smooth(q_cold_wind_95)+labs(y="Partial effect",x="")+xlim(-2.5,2.5)+ylim(-1,1),
      plot_smooth(q_cold_prec_95,col="dark blue")+labs(y="",x="")+xlim(-2.5,2.5)+ylim(-1,1),
      labels=c("a","b")),text_grob("Cold spell model",size=12,face="bold")),
  
  annotate_figure(
    ggarrange(
      plot_smooth(q_jet_wind_95)+labs(y="",x="")+xlim(-2.5,2.5)+ylim(-1,1),
      plot_smooth(q_jet_prec_95,col="dark blue")+labs(y="",x="")+xlim(-2.5,2.5)+ylim(-1,1),
      labels = c("e","f")
    ),text_grob("Cold spell and jet stream model",size=12,face="bold")
  ),nrow=1),text_grob("95th percentile",size=12,face="bold")
)


b<-annotate_figure(ggarrange(
  ggarrange(
    plot_smooth(q_cold_wind_99)+xlim(-2.5,2.5)+ylim(-1,1)+labs(y="Partial effect",x="lag -2 2m temp. an. (K), std. variable"),
    plot_smooth(q_cold_prec_99,col="dark blue")+xlim(-2.5,2.5)+ylim(-1,1)+labs(y="",x="lag -2 2m temp. an. (K), std. variable"),
    labels = c("c","d")
  ),
  ggarrange(
    plot_smooth(q_jet_wind_99)+xlim(-2.5,2.5)+ylim(-1,1)+labs(y="",x="lag -2 2m temp. an. (K), std. variable"),
    plot_smooth(q_jet_prec_99,col="dark blue")+xlim(-2.5,2.5)+ylim(-1,1)+labs(y="",x="lag -2 2m temp. an. (K), std. variable"),
    labels = c("g","h")
  ),nrow = 1),text_grob("99th percentile",size=12,face="bold")
)

annotate_figure(
  ggarrange(a,b,nrow=2),text_grob("Partial effect of temperature at 2m height in North America on near-surface weather in Western Europe",size=15,face="bold"))


############QREG########

#wind

l_base_wind_95<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(qu_95_wind),tau=0.95,tot_train)
l_cold_wind_95<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(US_temp_lag2)+as.numeric(qu_95_wind),tau=0.95,tot_train)
l_jet_wind_95<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(jet_proximity_lag1)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(NAO_lag1)+as.numeric(qu_95_wind),tau=0.95,tot_train)

l_base_wind_99<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(qu_99_wind),tau=0.99,tot_train)
l_cold_wind_99<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(US_temp_lag2)+as.numeric(qu_99_wind),tau=0.99,tot_train)
l_jet_wind_99<-rq(wind~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(jet_proximity_lag1)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(NAO_lag1)+as.numeric(qu_99_wind),tau=0.99,tot_train)


pred_l_base_wind_99<-predict(l_base_wind_99,tot_test)
pred_l_base_wind_95<-predict(l_base_wind_95,tot_test)
pred_l_cold_wind_99<-predict(l_cold_wind_99,tot_test)
pred_l_cold_wind_95<-predict(l_cold_wind_95,tot_test)
pred_l_jet_wind_99<-predict(l_jet_wind_99,tot_test)
pred_l_jet_wind_95<-predict(l_jet_wind_95,tot_test)


p_l_base_wind_99<-Quant_eval(pred_l_base_wind_99,wind,0.99,tot_test)
p_l_base_99<-overpred(p_l_base_wind_99,0.99,direc=-1)
p_l_base_wind_95<-Quant_eval(pred_l_base_wind_95,wind,0.95,tot_test)
p_l_base_95<-overpred(p_l_base_wind_95,0.95,direc=-1)
p_l_cold_wind_99<-Quant_eval(pred_l_cold_wind_99,wind,0.99,tot_test)
p_l_cold_99<-overpred(p_l_cold_wind_99,0.99,direc=-1)
p_l_cold_wind_95<-Quant_eval(pred_l_cold_wind_95,wind,0.95,tot_test)
p_l_cold_95<-overpred(p_l_cold_wind_95,0.95,direc=-1)
p_l_jet_wind_99<-Quant_eval(pred_l_jet_wind_99,wind,0.99,tot_test)
p_l_jet_99<-overpred(p_l_jet_wind_99,0.99,direc=-1)
p_l_jet_wind_95<-Quant_eval(pred_l_jet_wind_95,wind,0.95,tot_test)
p_l_jet_95<-overpred(p_l_jet_wind_95,0.95,direc=-1)

a<-annotate_figure(ggarrange(p_l_base_95,
                             p_l_cold_95,
                             p_l_jet_95,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_l_base_99,
                             p_l_cold_99,
                             p_l_jet_99,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM bias: Daily 10m Wind Speed",face="bold", size=15)
)

R2_l_base_wind_99<-Pseudo_R2_point(pred_l_base_wind_99,wind,0.99,tot_train,tot_test)
p_l_base_99_R2<-Plot_R2(R2_l_base_wind_99,0.99,up_lim = 0.31)

R2_l_base_wind_95<-Pseudo_R2_point(pred_l_base_wind_95,wind,0.95,tot_train,tot_test)
p_l_base_95_R2<-Plot_R2(R2_l_base_wind_95,0.95,up_lim = 0.31)

R2_l_cold_wind_99<-Pseudo_R2_point(pred_l_cold_wind_99,wind,0.99,tot_train,tot_test)
p_l_cold_99_R2<-Plot_R2(R2_l_cold_wind_99,0.99,up_lim = 0.31)

R2_l_cold_wind_95<-Pseudo_R2_point(pred_l_cold_wind_95,wind,0.95,tot_train,tot_test)
p_l_cold_95_R2<-Plot_R2(R2_l_cold_wind_95,0.95,up_lim = 0.31)

R2_l_jet_wind_99<-Pseudo_R2_point(pred_l_jet_wind_99,wind,0.99,tot_train,tot_test)
p_l_jet_99_R2<-Plot_R2(R2_l_jet_wind_99,0.99,up_lim = 0.31)

R2_l_jet_wind_95<-Pseudo_R2_point(pred_l_jet_wind_95,wind,0.95,tot_train,tot_test)
p_l_jet_95_R2<-Plot_R2(R2_l_jet_wind_95,0.95,up_lim = 0.31)



a<-annotate_figure(ggarrange(p_l_base_95_R2,
                             p_l_cold_95_R2,
                             p_l_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_l_base_99_R2,
                             p_l_cold_99_R2,
                             p_l_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM performance: Daily 10m Wind Speed",face="bold", size=15)
)


summary(l_base_wind_95)
summary(l_cold_wind_95)
summary(l_jet_wind_95)

summary(l_base_wind_99)
summary(l_cold_wind_99)
summary(l_jet_wind_99)

Pseudo_R2_tot(pred_l_jet_wind_95,wind,0.95,tot_train,tot_test)
Pseudo_R2_tot(pred_l_jet_wind_99,wind,0.99,tot_train,tot_test)



#prec

l_base_prec_95<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(qu_95_prec),tau=0.95,tot_train)
l_cold_prec_95<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(US_temp_lag2)+as.numeric(qu_95_prec),tau=0.95,tot_train)
l_jet_prec_95<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(jet_proximity_lag1)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(NAO_lag1)+as.numeric(qu_95_prec),tau=0.95,tot_train)

l_base_prec_99<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(qu_99_prec),tau=0.99,tot_train)
l_cold_prec_99<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(US_temp_lag2)+as.numeric(qu_99_prec),tau=0.99,tot_train)
l_jet_prec_99<-rq(prec~as.numeric(lat)*as.numeric(lon)+as.numeric(year)+month+as.numeric(jet_proximity_lag1)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(NAO_lag1)+as.numeric(qu_99_prec),tau=0.99,tot_train)

pred_l_base_prec_99<-predict(l_base_prec_99,tot_test)
pred_l_base_prec_95<-predict(l_base_prec_95,tot_test)
pred_l_cold_prec_99<-predict(l_cold_prec_99,tot_test)
pred_l_cold_prec_95<-predict(l_cold_prec_95,tot_test)
pred_l_jet_prec_99<-predict(l_jet_prec_99,tot_test)
pred_l_jet_prec_95<-predict(l_jet_prec_95,tot_test)


p_l_base_prec_99<-Quant_eval(pred_l_base_prec_99,prec,0.99,tot_test)
p_l_base_99<-overpred(p_l_base_prec_99,0.99,direc=-1)
p_l_base_prec_95<-Quant_eval(pred_l_base_prec_95,prec,0.95,tot_test)
p_l_base_95<-overpred(p_l_base_prec_95,0.95,direc=-1)
p_l_cold_prec_99<-Quant_eval(pred_l_cold_prec_99,prec,0.99,tot_test)
p_l_cold_99<-overpred(p_l_cold_prec_99,0.99,direc=-1)
p_l_cold_prec_95<-Quant_eval(pred_l_cold_prec_95,prec,0.95,tot_test)
p_l_cold_95<-overpred(p_l_cold_prec_95,0.95,direc=-1)
p_l_jet_prec_99<-Quant_eval(pred_l_jet_prec_99,prec,0.99,tot_test)
p_l_jet_99<-overpred(p_l_jet_prec_99,0.99,direc=-1)
p_l_jet_prec_95<-Quant_eval(pred_l_jet_prec_95,prec,0.95,tot_test)
p_l_jet_95<-overpred(p_l_jet_prec_95,0.95,direc=-1)

a<-annotate_figure(ggarrange(p_l_base_95,
                             p_l_cold_95,
                             p_l_jet_95,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_l_base_99,
                             p_l_cold_99,
                             p_l_jet_99,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM bias: Precipitation rate",face="bold", size=15)
)

R2_l_base_prec_99<-Pseudo_R2_point(pred_l_base_prec_99,prec,0.99,tot_train,tot_test)
p_l_base_99_R2<-Plot_R2(R2_l_base_prec_99,0.99,up_lim = 0.31, opt="G")

R2_l_base_prec_95<-Pseudo_R2_point(pred_l_base_prec_95,prec,0.95,tot_train,tot_test)
p_l_base_95_R2<-Plot_R2(R2_l_base_prec_95,0.95,up_lim = 0.31, opt="G")

R2_l_cold_prec_99<-Pseudo_R2_point(pred_l_cold_prec_99,prec,0.99,tot_train,tot_test)
p_l_cold_99_R2<-Plot_R2(R2_l_cold_prec_99,0.99,up_lim = 0.31, opt="G")

R2_l_cold_prec_95<-Pseudo_R2_point(pred_l_cold_prec_95,prec,0.95,tot_train,tot_test)
p_l_cold_95_R2<-Plot_R2(R2_l_cold_prec_95,0.95,up_lim = 0.31, opt="G")

R2_l_jet_prec_99<-Pseudo_R2_point(pred_l_jet_prec_99,prec,0.99,tot_train,tot_test)
p_l_jet_99_R2<-Plot_R2(R2_l_jet_prec_99,0.99,up_lim = 0.31, opt="G")

R2_l_jet_prec_95<-Pseudo_R2_point(pred_l_jet_prec_95,prec,0.95,tot_train,tot_test)
p_l_jet_95_R2<-Plot_R2(R2_l_jet_prec_95,0.95,up_lim = 0.31, opt="G")



a<-annotate_figure(ggarrange(p_l_base_95_R2,
                             p_l_cold_95_R2,
                             p_l_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_l_base_99_R2,
                             p_l_cold_99_R2,
                             p_l_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM performance: Daily precipitation",face="bold", size=15)
)



summary(l_base_prec_95)
summary(l_cold_prec_95)
summary(l_jet_prec_95)

summary(l_base_prec_99)
summary(l_cold_prec_99)
summary(l_jet_prec_99)


Pseudo_R2_tot(pred_l_base_prec_95,prec,0.95,tot_train,tot_test)
Pseudo_R2_tot(pred_l_cold_prec_95,prec,0.95,tot_train,tot_test)
Pseudo_R2_tot(pred_l_jet_prec_95,prec,0.95,tot_train,tot_test)


#############Multiplots - Plot 8-11###########


#WIND

R2_l_base_multi_99<-Pseudo_R2_point_multi(pred_q_base_wind_99,pred_l_base_wind_99,wind,0.99,tot_train,tot_test)
p_l_base_99_R2<-Plot_R2(R2_l_base_multi_99,0.99,up_lim = 0.11, low_lim = -0.11)

R2_l_cold_multi_99<-Pseudo_R2_point_multi(pred_q_cold_wind_99,pred_l_cold_wind_99,wind,0.99,tot_train,tot_test)
p_l_cold_99_R2<-Plot_R2(R2_l_cold_multi_99,0.99,up_lim = 0.11, low_lim = -0.11)

R2_l_jet_multi_99<-Pseudo_R2_point_multi(pred_q_jet_wind_99,pred_l_jet_wind_99,wind,0.99,tot_train,tot_test)
p_l_jet_99_R2<-Plot_R2(R2_l_jet_multi_99,0.99,up_lim = 0.11, low_lim = -0.11)


Pseudo_R2_tot_multi(pred_q_base_wind_99,pred_l_base_wind_99,wind,0.99,tot_train,tot_test)
Pseudo_R2_tot_multi(pred_q_cold_wind_99,pred_l_cold_wind_99,wind,0.99,tot_train,tot_test)
Pseudo_R2_jet_99_QGAM_QREG_wind<-Pseudo_R2_tot_multi(pred_q_jet_wind_99,pred_l_jet_wind_99,wind,0.99,tot_train,tot_test)
Pseudo_R2_jet_99_QGAM_QREG_wind


R2_l_base_multi_95<-Pseudo_R2_point_multi(pred_q_base_wind_95,pred_l_base_wind_95,wind,0.95,tot_train,tot_test)
p_l_base_95_R2<-Plot_R2(R2_l_base_multi_95,0.95,up_lim = 0.11, low_lim = -0.11)


R2_l_cold_multi_95<-Pseudo_R2_point_multi(pred_q_cold_wind_95,pred_l_cold_wind_95,wind,0.95,tot_train,tot_test)
p_l_cold_95_R2<-Plot_R2(R2_l_cold_multi_95,0.95,up_lim = 0.11, low_lim = -0.11)


R2_l_jet_multi_95<-Pseudo_R2_point_multi(pred_q_jet_wind_95,pred_l_jet_wind_95,wind,0.95,tot_train,tot_test)
p_l_jet_95_R2<-Plot_R2(R2_l_jet_multi_95,0.95,up_lim = 0.11, low_lim = -0.11)


a<-annotate_figure(ggarrange(p_l_base_95_R2,
                             p_l_cold_95_R2,
                             p_l_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))
b<-annotate_figure(ggarrange(p_l_base_99_R2,
                             p_l_cold_99_R2,
                             p_l_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM vs QREG performance: Daily 10m Wind Speed",face="bold", size=15)
)

Pseudo_R2_tot_multi(pred_q_base_wind_95,pred_l_base_wind_95,wind,0.95,tot_train,tot_test)
Pseudo_R2_tot_multi(pred_q_cold_wind_95,pred_l_cold_wind_95,wind,0.95,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_QREG_wind<-Pseudo_R2_tot_multi(pred_q_jet_wind_95,pred_l_jet_wind_95,wind,0.95,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_QREG_wind


#PRECIPITATION

R2_l_base_multi_99<-Pseudo_R2_point_multi(pred_q_base_prec_99,pred_l_base_prec_99,prec,0.99,tot_train,tot_test)
p_l_base_99_R2<-Plot_R2(R2_l_base_multi_99,0.99,up_lim = 0.11, low_lim = -0.11, opt="G")


R2_l_cold_multi_99<-Pseudo_R2_point_multi(pred_q_cold_prec_99,pred_l_cold_prec_99,prec,0.99,tot_train,tot_test)
p_l_cold_99_R2<-Plot_R2(R2_l_cold_multi_99,0.99,up_lim = 0.11, low_lim = -0.11,opt="G")


R2_l_jet_multi_99<-Pseudo_R2_point_multi(pred_q_jet_prec_99,pred_l_jet_prec_99,prec,0.99,tot_train,tot_test)
p_l_jet_99_R2<-Plot_R2(R2_l_jet_multi_99,0.99,up_lim = 0.11, low_lim = -0.11,opt="G")

Pseudo_R2_tot_multi(pred_q_base_prec_99,pred_l_base_prec_99,prec,0.99,tot_train,tot_test)
Pseudo_R2_tot_multi(pred_q_cold_prec_99,pred_l_cold_prec_99,prec,0.99,tot_train,tot_test)
Pseudo_R2_jet_99_QGAM_QREG_prec<-Pseudo_R2_tot_multi(pred_q_jet_prec_99,pred_l_jet_prec_99,prec,0.99,tot_train,tot_test)
Pseudo_R2_jet_99_QGAM_QREG_prec


R2_l_base_multi_95<-Pseudo_R2_point_multi(pred_q_base_prec_95,pred_l_base_prec_95,prec,0.95,tot_train,tot_test)
p_l_base_95_R2<-Plot_R2(R2_l_base_multi_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")


R2_l_cold_multi_95<-Pseudo_R2_point_multi(pred_q_cold_prec_95,pred_l_cold_prec_95,prec,0.95,tot_train,tot_test)
p_l_cold_95_R2<-Plot_R2(R2_l_cold_multi_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")


R2_l_jet_multi_95<-Pseudo_R2_point_multi(pred_q_jet_prec_95,pred_l_jet_prec_95,prec,0.95,tot_train,tot_test)
p_l_jet_95_R2<-Plot_R2(R2_l_jet_multi_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")

Pseudo_R2_tot_multi(pred_q_base_prec_95,pred_l_base_prec_95,prec,0.95,tot_train,tot_test)
Pseudo_R2_tot_multi(pred_q_cold_prec_95,pred_l_cold_prec_95,prec,0.95,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_QREG_prec<-Pseudo_R2_tot_multi(pred_q_jet_prec_95,pred_l_jet_prec_95,prec,0.95,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_QREG_prec

a<-annotate_figure(ggarrange(p_l_base_95_R2,
                             p_l_cold_95_R2,
                             p_l_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_l_base_99_R2,
                             p_l_cold_99_R2,
                             p_l_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM vs QREG performance: Daily Precipitation",face="bold", size=15)
)


################ POT ################

#Wind

da<-POT_declust_out(tot_train,wind,0.8,over=TRUE)
da$dat$excesses<-da$dat$key_var-da$threshold

POT_base_wind_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(qu_99_wind),~1),family="gpd",da$dat)
R2_POT_base_wind_99<-Pseudo_R2_point_EVT(POT_base_wind_99,wind,0.99,train=tot_train,test=tot_test)
p_POT_base_99_R2<-Plot_R2(R2_POT_base_wind_99,0.99,up_lim = 0.21, low_lim = -0.21)
Pseudo_R2_tot_multi_EVT(pred_q_base_wind_99,POT_base_wind_99,wind,0.99,tot_train,tot_test)

R2_multi_POT_base_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_base_wind_99,POT_base_wind_99,wind,0.99,tot_train,tot_test)
p_q_base_99_R2<-Plot_R2(R2_multi_POT_base_R2_99,0.99,up_lim = 0.11, low_lim = -0.11)

POT_cold_wind_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(qu_99_wind),~1),family="gpd",da$dat)
R2_POT_cold_wind_99<-Pseudo_R2_point_EVT(POT_cold_wind_99,wind,0.99,tot_train,tot_test)
p_POT_cold_99_R2<-Plot_R2(R2_POT_cold_wind_99,0.99,up_lim = 0.62)
Pseudo_R2_tot_multi_EVT(pred_q_cold_wind_99,POT_cold_wind_99,wind,0.99,tot_train,tot_test)

R2_multi_POT_cold_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_cold_wind_99,POT_cold_wind_99,wind,0.99,tot_train,tot_test)
p_q_cold_99_R2<-Plot_R2(R2_multi_POT_cold_R2_99,0.99,up_lim = 0.11, low_lim = -0.11)

da$dat<-da$dat%>%
  drop_na(NAO_lag1)

POT_jet_wind_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(lat_jet_lag1)+as.numeric(NAO_lag1)+as.numeric(jet_proximity_lag1)+as.numeric(qu_99_wind),~1),family="gpd",da$dat)
R2_POT_jet_wind_99<-Pseudo_R2_point_EVT(POT_jet_wind_99,wind,0.99,tot_train,tot_test)
p_POT_jet_99_R2<-Plot_R2(R2_POT_jet_wind_99,0.99,up_lim = 0.62)
Pseudo_R2_jet_99_QGAM_POT_wind<-Pseudo_R2_tot_multi_EVT(pred_q_jet_wind_99,POT_jet_wind_99,wind,0.99,tot_train,tot_test)

R2_multi_POT_jet_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_jet_wind_99,POT_jet_wind_99,wind,0.99,tot_train,tot_test)
p_q_jet_99_R2<-Plot_R2(R2_multi_POT_jet_R2_99,0.99,up_lim = 0.11, low_lim = -0.11)

Pseudo_R2_jet_99_QGAM_POT_wind

da<-POT_declust_out(tot_train,wind,0.8,over=TRUE)
da$dat$excesses<-da$dat$key_var-da$threshold

POT_base_wind_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(qu_95_wind),~1),family="gpd",da$dat)
R2_POT_base_wind_95<-Pseudo_R2_point_EVT(POT_base_wind_95,wind,0.95,tot_train,tot_test)
p_POT_base_95_R2<-Plot_R2(R2_POT_base_wind_95,0.95,up_lim = 0.21, low_lim = -0.21)
Pseudo_R2_tot_multi_EVT(pred_q_base_wind_95,POT_base_wind_95,wind,0.95,tot_train,tot_test)

R2_multi_POT_base_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_base_wind_95,POT_base_wind_95,wind,0.95,tot_train,tot_test)
p_q_base_95_R2<-Plot_R2(R2_multi_POT_base_R2_95,0.95,up_lim = 0.11, low_lim = -0.11)

POT_cold_wind_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(qu_95_wind),~1),family="gpd",da$dat)
R2_POT_cold_wind_95<-Pseudo_R2_point_EVT(POT_cold_wind_95,wind,0.95,tot_train,tot_test)
p_POT_cold_95_R2<-Plot_R2(R2_POT_cold_wind_95,0.95,up_lim = 0.62)
Pseudo_R2_tot_multi_EVT(pred_q_cold_wind_95,POT_cold_wind_95,wind,0.95,tot_train,tot_test)

R2_multi_POT_cold_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_cold_wind_95,POT_cold_wind_95,wind,0.95,tot_train,tot_test)
p_q_cold_95_R2<-Plot_R2(R2_multi_POT_cold_R2_95,0.95,up_lim = 0.11, low_lim = -0.11)

da$dat<-da$dat%>%
  drop_na(NAO_lag1)

POT_jet_wind_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(lat_jet_lag1)+as.numeric(NAO_lag1)+as.numeric(jet_proximity_lag1)+as.numeric(qu_95_wind),~1),family="gpd",da$dat)
R2_POT_jet_wind_95<-Pseudo_R2_point_EVT(POT_jet_wind_95,wind,0.95,tot_train,tot_test)
p_POT_jet_95_R2<-Plot_R2(R2_POT_jet_wind_95,0.95,up_lim = 0.62)
Pseudo_R2_jet_95_QGAM_POT_wind<-Pseudo_R2_tot_multi_EVT(pred_q_jet_wind_95,POT_jet_wind_95,wind,0.95,tot_train,tot_test)
Pseudo_R2_jet_95_QGAM_POT_wind

R2_multi_POT_jet_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_jet_wind_95,POT_jet_wind_95,wind,0.95,tot_train,tot_test)
p_q_jet_95_R2<-Plot_R2(R2_multi_POT_jet_R2_95,0.95,up_lim = 0.11, low_lim = -0.11)

a<-annotate_figure(ggarrange(p_q_base_95_R2,
                             p_q_cold_95_R2,
                             p_q_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_q_base_99_R2,
                             p_q_cold_99_R2,
                             p_q_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM vs POT performance: Daily 10m Wind Speed",face="bold", size=15)
)

#Prec

da<-POT_declust_out(tot_train,prec,0.8,over=TRUE,declust=TRUE)
da$dat$excesses<-da$dat$key_var-da$threshold

POT_base_prec_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(qu_99_prec),~1),family="gpd",da$dat)
R2_POT_base_prec_99<-Pseudo_R2_point_EVT(POT_base_prec_99,prec,0.99,tot_train,tot_test)
p_POT_base_99_R2<-Plot_R2(R2_POT_base_prec_99,0.99,up_lim = 0.21, low_lim = -0.21)
Pseudo_R2_tot_multi_EVT(pred_q_base_prec_99,POT_base_prec_99,prec,0.99,tot_train,tot_test)

R2_multi_POT_base_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_base_prec_99,POT_base_prec_99,prec,0.99,tot_train,tot_test)
p_q_base_99_R2<-Plot_R2(R2_multi_POT_base_R2_99,0.99,up_lim = 0.11, low_lim = -0.11,opt="G")

POT_cold_prec_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(qu_99_prec),~1),family="gpd",da$dat)
R2_POT_cold_prec_99<-Pseudo_R2_point_EVT(POT_cold_prec_99,prec,0.99,tot_train,tot_test)
p_POT_cold_99_R2<-Plot_R2(R2_POT_cold_prec_99,0.99,up_lim = 0.62)
Pseudo_R2_tot_multi_EVT(pred_q_cold_prec_99,POT_cold_prec_99,prec,0.99,tot_train,tot_test)


R2_multi_POT_cold_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_cold_prec_99,POT_cold_prec_99,prec,0.99,tot_train,tot_test)
p_q_cold_99_R2<-Plot_R2(R2_multi_POT_cold_R2_99,0.99,up_lim = 0.11, low_lim = -0.11,opt="G")
Pseudo_R2_tot_multi_EVT(pred_q_cold_prec_99,POT_cold_prec_99,prec,0.99,tot_train,tot_test)

da$dat<-da$dat%>%
  drop_na(NAO_lag1)

POT_jet_prec_99<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(lat_jet_lag1)+as.numeric(NAO_lag1)+as.numeric(jet_proximity_lag1)+as.numeric(qu_99_prec),~1),family="gpd",da$dat)
R2_POT_jet_prec_99<-Pseudo_R2_point_EVT(POT_jet_prec_99,prec,0.99,tot_train,tot_test)
p_POT_jet_99_R2<-Plot_R2(R2_POT_jet_prec_99,0.99,up_lim = 0.62)
Pseudo_R2_jet_99_QGAM_POT_prec<-Pseudo_R2_tot_multi_EVT(pred_q_jet_prec_99,POT_jet_prec_99,prec,0.99,tot_train,tot_test)

R2_multi_POT_jet_R2_99<-Pseudo_R2_point_multi_EVT(pred_q_jet_prec_99,POT_jet_prec_99,prec,0.99,tot_train,tot_test)
p_q_jet_99_R2<-Plot_R2(R2_multi_POT_jet_R2_99,0.99,up_lim = 0.11, low_lim = -0.11,opt="G")
Pseudo_R2_jet_99_QGAM_POT_prec

da<-POT_declust_out(tot_train,prec,0.8,over=TRUE,declust=TRUE)
da$dat$excesses<-da$dat$key_var-da$threshold

POT_base_prec_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(qu_95_prec),~1),family="gpd",da$dat)
R2_POT_base_prec_95<-Pseudo_R2_point_EVT(POT_base_prec_95,prec,0.95,tot_train,tot_test)
p_POT_base_95_R2<-Plot_R2(R2_POT_base_prec_95,0.95,up_lim = 0.21, low_lim = -0.21)
Pseudo_R2_tot_multi_EVT(pred_q_base_prec_95,POT_base_prec_95,prec,0.95,tot_train,tot_test)

R2_multi_POT_base_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_base_prec_95,POT_base_prec_95,prec,0.95,tot_train,tot_test)
p_q_base_95_R2<-Plot_R2(R2_multi_POT_base_R2_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")


POT_cold_prec_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(qu_95_prec),~1),family="gpd",da$dat)
R2_POT_cold_prec_95<-Pseudo_R2_point_EVT(POT_cold_prec_95,prec,0.95,tot_train,tot_test)
p_POT_cold_95_R2<-Plot_R2(R2_POT_cold_prec_95,0.95,up_lim = 0.62)
Pseudo_R2_tot_multi_EVT(pred_q_cold_prec_95,POT_cold_prec_95,prec,0.95,tot_train,tot_test)

R2_multi_POT_cold_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_cold_prec_95,POT_cold_prec_95,prec,0.95,tot_train,tot_test)
p_q_cold_95_R2<-Plot_R2(R2_multi_POT_cold_R2_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")

da$dat<-da$dat%>%
  drop_na(NAO_lag1)

POT_jet_prec_95<-evgam(list(excesses~as.numeric(lat)*as.numeric(lon)+month+as.numeric(year)+as.numeric(US_temp_lag2)+as.numeric(jet_strength_lag1)+as.numeric(lat_jet_lag1)+as.numeric(NAO_lag1)+as.numeric(jet_proximity_lag1)+as.numeric(qu_95_prec),~1),family="gpd",da$dat)
R2_POT_jet_prec_95<-Pseudo_R2_point_EVT(POT_jet_prec_95,prec,0.95,tot_train,tot_test)
p_POT_jet_95_R2<-Plot_R2(R2_POT_jet_prec_95,0.95,up_lim = 0.62)
Pseudo_R2_jet_95_QGAM_POT_prec<-Pseudo_R2_tot_multi_EVT(pred_q_jet_prec_95,POT_jet_prec_95,prec,0.95,tot_train,tot_test)

R2_multi_POT_jet_R2_95<-Pseudo_R2_point_multi_EVT(pred_q_jet_prec_95,POT_jet_prec_95,prec,0.95,tot_train,tot_test)
p_q_jet_95_R2<-Plot_R2(R2_multi_POT_jet_R2_95,0.95,up_lim = 0.11, low_lim = -0.11,opt="G")
Pseudo_R2_jet_95_QGAM_POT_prec

a<-annotate_figure(ggarrange(p_q_base_95_R2,
                             p_q_cold_95_R2,
                             p_q_jet_95_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("a","b","c")),
                   text_grob("95th percentile",face="bold", size=12))

b<-annotate_figure(ggarrange(p_q_base_99_R2,
                             p_q_cold_99_R2,
                             p_q_jet_99_R2,
                             common.legend = TRUE,legend="right",nrow=1,labels=c("d","e","f")),
                   text_grob("99th percentile",face="bold", size=12))

annotate_figure(
  ggarrange(a,b,nrow=2),
  text_grob("QGAM vs POT performance: Daily Precipitation",face="bold", size=15)
)


##########Summary table###########

R2_wind_tab<-matrix(c(Pseudo_R2_jet_95_QGAM_wind,Pseudo_R2_jet_99_QGAM_wind,Pseudo_R2_jet_95_QGAM_POT_wind,Pseudo_R2_jet_99_QGAM_POT_wind,Pseudo_R2_jet_95_QGAM_QREG_wind,Pseudo_R2_jet_99_QGAM_QREG_wind),nrow=3,ncol=2,byrow=TRUE)%>%
  as.data.frame()

colnames(R2_wind_tab)<-c("95th percentile","99th percentile")
rownames(R2_wind_tab)<-c("Quantile of the seasonal climatology","POT","QREG")

R2_wind_tab%>%
  xtable(digits=4)

Pseudo_R2_jet_95_QGAM_wind
Pseudo_R2_jet_99_QGAM_wind
Pseudo_R2_jet_95_QGAM_QREG_wind
Pseudo_R2_jet_99_QGAM_QREG_wind

R2_prec_tab<-matrix(c(Pseudo_R2_jet_95_QGAM_prec,Pseudo_R2_jet_99_QGAM_prec,Pseudo_R2_jet_95_QGAM_POT_prec,Pseudo_R2_jet_99_QGAM_POT_prec,Pseudo_R2_jet_95_QGAM_QREG_prec,Pseudo_R2_jet_99_QGAM_QREG_prec),nrow=3,ncol=2,byrow=TRUE)%>%
  as.data.frame()

colnames(R2_prec_tab)<-c("95th percentile","99th percentile")
rownames(R2_prec_tab)<-c("Quantile of the seasonal climatology","POT","QREG")

R2_prec_tab%>%
  xtable(digits=4)

Pseudo_R2_jet_95_QGAM_prec
Pseudo_R2_jet_99_QGAM_prec
Pseudo_R2_jet_95_QGAM_QREG_prec
Pseudo_R2_jet_99_QGAM_QREG_prec



