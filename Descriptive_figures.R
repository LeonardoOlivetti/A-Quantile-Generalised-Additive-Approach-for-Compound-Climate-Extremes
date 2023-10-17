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
library(LSTS)
library(zoo)
library(lubridate)
library(tsibble)
library(itsadug)
library(extRemes)
library(evgam)
library(scales)


#load_data
dat<-readRDS("data/Y_var.Rda")%>%
  mutate(date=as.character(date))

out<-readRDS("data/X_var.Rda")

#Cleaning

tot<-left_join(out,dat)%>%
  mutate(fac_lat=as.character(lat),fac_lon=as.character(lon),fac_lat_jet=as.character(lat_jet))%>%
  mutate(date=as_date(date),jet_proximity=abs(lat-lat_jet),lon=ifelse(lon<5,lon+360,lon))%>%
  mutate(year=year(date))%>%
  as_tsibble(index=date,key=c(lat,lon))%>%
  fill_gaps()%>%
  tk_augment_lags(c(jet_strength,NAO,US_temp,lat_jet,jet_proximity,wind,prec), .lags = c(1:15), .names = "auto")%>%
  tk_augment_leads(c(jet_strength,NAO,US_temp,lat_jet,jet_proximity,wind,prec), .lags = c(-1:-15), .names = "auto")%>%
  mutate(time_cont=as.integer(as.Date(date, "%Y-%m-%d")-as.Date("1959-01-01")))%>%
  filter(month(date)<3 | month(date)>10)%>%
  mutate(month=month(date))%>%
  mutate(across(starts_with("fac"),as.numeric))%>%
  mutate(across(starts_with("fac"),round,digits=2))%>%
  filter(date>"1959-3-31"& date<"2022-09-01")%>%
  group_by(year,month)%>%
  mutate(season_year=ceiling(cur_group_id()/4))%>%
  ungroup()%>%
  mutate(month=as.factor(month(date)))%>%
  mutate(date=as.character(date))

gc()

###########Figure 2 - right side ###############

wind_comp_mean<-POT_cold_spell_monte_mean(tot,US_temp,wind,0.05,B=20000)
prec_comp_mean<-POT_cold_spell_monte_mean(tot,US_temp,prec,0.05,B=20000)
wind_comp_95<-POT_cold_spell_monte_qu(tot,US_temp,wind,qu_out=0.95,B=20000)
prec_comp_95<-POT_cold_spell_monte_qu(tot,US_temp,prec,0.05,qu_out=0.95,B=20000)
wind_comp_99<-POT_cold_spell_monte_qu(tot,US_temp,wind,0.05,qu_out=0.99,B=20000)
prec_comp_99<-POT_cold_spell_monte_qu(tot,US_temp,prec,0.05,qu_out=0.99,B=20000)

c<-wind_comp_mean$pl+labs(x="lag(days)",y="daily 10m wind an.(m/s)")+theme(axis.title.y = element_text (size=11))
d<-prec_comp_mean$pl+labs(x="lag(days)",y="daily prec. an.(mm/day)")+theme(axis.title.y = element_text( size = 11))
cd<-annotate_figure(
  ggarrange(c,d,hjust=-0.4,labels=c("c","d")),text_grob("Mean",size=16,face="bold")
)

e<-wind_comp_95$pl+labs(x="lag(days)",y="daily 10m wind an.(m/s)")+theme(axis.title.y = element_text (size=11))
f<-prec_comp_95$pl+labs(x="lag(days)",y="daily prec. an.(mm/day)")+theme(axis.title.y = element_text( size = 11))
ef<-annotate_figure(ggarrange(e,f,hjust=-0.4,labels=c("e","f")),
                    text_grob("95th quantile",size=16,face="bold"))

g<-wind_comp_99$pl+labs(x="lag(days)",y="daily 10m wind an.(m/s)")+theme(axis.title.y = element_text (size=11))
h<-prec_comp_99$pl+labs(x="lag(days)",y="daily prec. an.(mm/day)")+theme(axis.title.y = element_text( size = 11))
gh<-annotate_figure(ggarrange(g,h,hjust=-0.4,labels=c("g","h")),
                    text_grob("99th quantile",size=16,face="bold"))

cdefgh<-annotate_figure(
  ggarrange(cd,ef,gh,nrow=3))

#########Figure 2 - Non-overlapping densities######

library(cowplot)
library(ggridges)


a<-wind_comp_mean$dat$wind_lead2
b<-wind_comp_mean$dat$wind_lead1
c<-wind_comp_mean$dat$wind
d<-wind_comp_mean$dat$wind_lag2
e<-tot$wind

a<-c(a,rep(NA, length(e)-length(a)))
b<-c(b,rep(NA, length(e)-length(b)))
c<-c(c,rep(NA, length(e)-length(c)))
d<-c(d,rep(NA, length(e)-length(d)))

da<-data.frame(a,b,c,d,e)%>%
  rename("climatology"=e,"two days before cold spell"=d, "same day as cold spell"=c, "one day after cold spell"=b, "two days after cold spell"=a)%>%
  pivot_longer(cols=1:5)%>%
  mutate(name=factor(name,levels=c("two days after cold spell","one day after cold spell","same day as cold spell","two days before cold spell","climatology")))


wind_dens<-da%>%
  ggplot(aes(x = value, y = name)) +
  stat_density_ridges(aes(fill = factor(after_stat(quantile))),
                      bandwidth = 0.5,
                      color="black",
                      geom = "density_ridges_gradient",
                      calc_ecdf = TRUE,
                      quantiles = c(0.95, 0.99),
                      scale=1,
                      quantile_fun = quantile,
                      quantile_lines = TRUE,
  ) +
  stat_density_ridges(
    bandwidth = 0.5,
    vline_color="black",
    vline_linetype="dashed",
    vline_size=0.8,
    geom = "density_ridges_gradient",
    quantile_lines = TRUE,
    scale=1,
    fill = NA,
    quantile_fun = mean,
  )+
  scale_fill_manual(
    name = "Area under the curve", values = c("#F6F3E7","yellow", "#FF0000A0"),
    labels = c("(0, 0.95]","(0.95, 0.99]", "(0.99, 1]")
  )+
  theme_minimal_vgrid()+
  scale_x_continuous(breaks = c(seq(-4,6,by=2)))+
  scale_y_discrete(expand = expand_scale(add = c(0.2, 1)))+
  labs(y="",x="daily precipitation anomalies (mm)")+
  theme(axis.title.x = element_text (size=11))+
  geom_vline(aes(xintercept = 50,linetype="dashed"),size=0.6)+
  scale_linetype_manual(name="", values=c(dashed="dashed"),labels=c("mean"))+
  scale_size_manual(name="")+
  guides(fill=guide_legend(order=1))+
  coord_cartesian(xlim=c(-4,7))


a<-prec_comp_mean$dat$prec_lead2
b<-prec_comp_mean$dat$prec_lead1
c<-prec_comp_mean$dat$prec
d<-prec_comp_mean$dat$prec_lag2
e<-tot$prec

a<-c(a,rep(NA, length(e)-length(a)))
b<-c(b,rep(NA, length(e)-length(b)))
c<-c(c,rep(NA, length(e)-length(c)))
d<-c(d,rep(NA, length(e)-length(d)))


da<-data.frame(a,b,c,d,e)%>%
  rename("climatology"=e,"two days before cold spell"=d, "same day as cold spell"=c, "one day after cold spell"=b, "two days after cold spell"=a)%>%
  pivot_longer(cols=1:5)%>%
  mutate(name=factor(name,levels=c("two days after cold spell","one day after cold spell","same day as cold spell","two days before cold spell","climatology")))

library(ggallin)

prec_dens<-da%>%
  ggplot(aes(x = value, y = name)) +
  stat_density_ridges(aes(fill = factor(after_stat(quantile))),
                      bandwidth=6,
    color="black",
    #size=0.7,
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.95, 0.99),
    scale=1,
    quantile_fun = quantile,
    quantile_lines = TRUE,
    rel_min_height=0.001
  ) +
  stat_density_ridges(
    bandwidth=6,
    vline_color="black",
    calc_ecdf = TRUE,
    vline_linetype="dashed",
    vline_size=0.8,
    geom = "density_ridges_gradient",
    quantile_lines = TRUE,
    scale=1,
    fill = NA,
    quantile_fun = mean,
    rel_min_height=0.001
  )+
  scale_fill_manual(
    name = "Area under the curve", values = c("#F6F3E7","yellow", "#FF0000A0"),
    labels = c("(0, 0.95]","(0.95, 0.99]", "(0.99, 1]")
  )+
  theme_minimal_vgrid()+
  scale_x_continuous(breaks = c(seq(-5,30,by=10)))+
  scale_y_discrete(expand = expand_scale(add = c(0.0, 0.0)))+
  labs(y="",x="daily precipitation anomalies (mm)")+
  theme(axis.title.x = element_text (size=11))+
  geom_vline(aes(xintercept = 150,linetype="dashed"),size=0.6)+
  scale_linetype_manual(name="", values=c(dashed="dashed"),labels=c("mean"))+
  scale_size_manual(name="")+
  guides(fill=guide_legend(order=1))+
  coord_cartesian(xlim=c(-5,30))


ab<-annotate_figure(ggarrange(wind_dens,prec_dens,nrow=2,labels=c("a","b"),common.legend = TRUE, legend = "bottom"))

#Figure 2 - panel 
annotate_figure(
  ggarrange(ab,cdefgh, ncol = 2),
  text_grob("Association between cold winter days in North America and near-surface weather in Western Europe", size=20,face="bold"
            )
  )

######## Figure 1 ##############

#Requires US_temps (2m temperature anomalies over North America) as input  

cold_days<-POT_cold_spell_monte_mean(tot,US_temp,wind,0.05,B=1)$dat%>%
  distinct(date)

cold_days<-cold_days%>%mutate(date=as.character(date))%>%pull(date)

Temps<-readRDS("data/daily_2m_temp.RDS")%>%
  filter(time>"1959-4-01"& time<"2022-09-01")%>%
  mutate(longitude=if_else(longitude<0, longitude+360,longitude))%>%
  rename(lon=longitude,lat=latitude,date=time,air=t2m)


Land_Sea_05<-hyper_tibble("data/land_sea_05.nc")%>%
  select(-time)%>%
  mutate(lsm=if_else(lsm>0.5,1,0))%>%
  rename(land=lsm,lon=longitude,lat=latitude)

US_temps<-Temps%>%
  filter(date>"1959-4-01"& date<"2022-09-01")%>%
  left_join(Land_Sea_05)%>%
  filter(lon>=260&lon<=290&lat>=30&lat<=45)%>%
  filter(land==1)%>%
  group_by(lat,lon)%>%
  mutate(roll_air=rollmean(air,7,fill=NA))%>%
  group_by(lat,lon,month(date),day(date))%>%
  mutate(clim_air=mean(roll_air,na.rm=TRUE))%>%
  mutate(air_an=air-clim_air)%>%
  ungroup()%>%
  mutate(date=as.character(date))%>%
  filter(date%in%cold_days)%>%
  group_by(lat,lon)%>%
  summarise(air_an=mean(air_an,na.rm=TRUE))%>%
  ungroup()

gc()

fig_1(US_temps,air_an,lon,lat,-105,-65,25,60)+
  labs(title= "Cold spell region, t2m anomalies",fill="t2m anomaly (K)",
       x="longitude",y="latitude")+
  theme(plot.title = element_text(hjust = 0.5))

