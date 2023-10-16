
#load_three-dimensional_nc
load_nc_three<-function(file){
  main<-hyper_tibble(file)
  main<-main%>%
    filter(lat>=10 & lat <80 & (lon<50 | lon>240))
  main$time<-as.POSIXct("1800-01-01 00:00")+as.difftime(main$time%>%unlist(),units="hours")
  main<-main%>%
    mutate(date=format(time,format="%Y-%m-%d"))%>%
    mutate(month=format(as.Date(date), format="%m"))%>%
    mutate(month=as.numeric(month))%>%
    filter(month<4 | month>9)%>%
    select(-month,-time)
  return(main)
}


#load_four-dimensional_nc
load_nc_four<-function(file){
  main<-hyper_tibble(file)
  main<-main%>%
    filter(level==925 | level ==500 |level==250)%>%
    filter(lat>20 & lat <80 & (lon<50 | lon>240))
  main$time<-as.POSIXct("1800-01-01 00:00")+as.difftime(main$time%>%unlist(),units="hours")
  main<-main%>%
    mutate(date=format(time,format="%Y-%m-%d"))%>%
    mutate(month=format(as.Date(date), format="%m"))%>%
    mutate(month=as.numeric(month))%>%
    filter(month<4 | month>9)%>%
    select(-month,-time)
  
  gc()
  return(main)
}

#load_four-dimensional_nc
load_nc_four_850<-function(file){
  main<-hyper_tibble(file)
  main<-main%>%
    filter(level==850)%>%
    filter(lat==30 | lat ==40 |lat==70 | lat==80)
  main$time<-as.POSIXct("1800-01-01 00:00")+as.difftime(main$time%>%unlist(),units="hours")
  main<-main%>%
    mutate(date=format(time,format="%Y-%m-%d"))%>%
    mutate(month=format(as.Date(date), format="%m"))%>%
    mutate(month=as.numeric(month))%>%
    filter(month<4 | month>9)%>%
    select(-month,-time)
  
  gc()
  return(main)
}


#Geom_line_month

Line_month<-function(var1,var2,var3,year_number,month_number){
  a<-tibble(var1,var2,var3)
  colnames(a)<-c("date","fitted","actual")
  
  a%>%
    filter(year(date)==year_number&month(date)==month_number)%>%
    pivot_longer(2:3)%>%
    ggplot()+
    geom_line(aes(x=as.Date(date),y=value, color=name),lwd=1)+
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    #facet_wrap(vars(year(date),month(date)))+
    theme_minimal()+
    labs(x="Date",y="Value",title="Predicted vs Actual",color="")
}

#thanks to https://www.statology.org/remove-outliers-from-multiple-columns-in-r/ for the clean function

outliers <- function(x) {
  
  Q1 <- quantile(x, probs=.25,na.rm=TRUE)
  Q3 <- quantile(x, probs=.75,na.rm=TRUE)
  iqr = Q3-Q1
  
  upper_limit = Q3 + (iqr*3)
  lower_limit = Q1 - (iqr*3)
  
  x > upper_limit | x < lower_limit
}

outliers_low <- function(x) {
  
  Q1 <- quantile(x, probs=.25,na.rm=TRUE)
  Q3 <- quantile(x, probs=.75,na.rm=TRUE)
  iqr = Q3-Q1
  
  upper_limit = Q3 + (iqr*1.5)
  lower_limit = Q1 - (iqr*1.5)
  
  x > upper_limit | x < lower_limit
}

#RMSE plot

RMSE_plot<-function(model,variable){
  srd_place<-test%>%
    mutate(pred=predict(model,test)%>%as.numeric())%>%
    mutate(srd_error=({{variable}}-pred)^2)%>%
    group_by(lon_fac,lat_fac)%>%
    select(srd_error)%>%
    summarise(srd_error=sqrt(mean(srd_error)))%>%
    rename(lon=lon_fac,lat=lat_fac,RMSE=srd_error)%>%
    mutate(lon=lon-360)
  
  mid=mean(srd_place$RMSE)
  
  world<-map_data("world")
  
  pl<-ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=srd_place,aes(lon,lat,color=RMSE),size=4)+
    coord_sf(xlim = c(-10, 2), ylim = c(35, 60), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_color_gradient2( midpoint=mid,low ="white",mid="orange",
                           high = "red", space = "Lab" )+
    #guides(color="none")+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
  
  return(pl)
}


#RMSE line plot

RMSE_line<-function(model,variable){
  srd_time<-test%>%
    mutate(pred=predict(model,test)%>%as.numeric())%>%
    mutate(srd_error=({{variable}}-pred)^2)%>%
    group_by(year(date))%>%
    select(srd_error)%>%
    summarise(srd_error=sqrt(mean(srd_error)))%>%
    rename(Year="year(date)",RMSE=srd_error)%>%
    ggplot(aes(x=Year,y=RMSE))+
    geom_line(color="black",lwd=1.5)+
    theme_bw()+
    ylim(0,2)
  
  return(srd_time)
}

#RMSE line facet plot

RMSE_line_facet<-function(model_1,variable_1,model_2,variable_2){
  srd_time<-test%>%
    mutate(pred_1=predict(model_1,test)%>%as.numeric(),
           pred_2=predict(model_2,test)%>%as.numeric())%>%
    mutate(srd_error_1=({{variable_1}}-pred_1)^2,
           srd_error_2=({{variable_2}}-pred_2)^2,)%>%
    group_by(year(date),lat_fac,lon_fac)%>%
    select(srd_error_1,srd_error_2)%>%
    summarise(srd_error_1=sqrt(mean(srd_error_1)),
              srd_error_2=sqrt(mean(srd_error_2)))%>%
    pivot_longer(srd_error_1:srd_error_2,names_to = "var",values_to = "RMSE")%>%
    rename(Year="year(date)")%>%
    ggplot(aes(x=Year,y=RMSE,color=var))+
    geom_line(lwd=1.5)+
    ylim(0,1.5)+
    facet_wrap(vars(lat_fac,lon_fac))+
    theme_bw()
}

#Jet Prediction

Jet_wind_pred_plot<-function(model,jet_pos){
  
  test_2<-test%>%
    mutate(lat_jet_lag1=jet_pos)
  
  srd_place<-test_2%>%
    mutate(pred=predict(model,test_2)%>%as.numeric())%>%
    group_by(lon_fac,lat_fac)%>%
    select(pred)%>%
    summarise(pred=mean(pred))%>%
    rename(lon=lon_fac,lat=lat_fac)%>%
    mutate(lon=lon-360)
  
  mid=mean(srd_place$pred)
  
  world<-map_data("world")
  
  pl<-ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=srd_place,aes(lon,lat,color=pred),size=4)+
    coord_sf(xlim = c(-10, 2), ylim = c(35, 60), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_color_gradient2( limits=c(0,2.5),midpoint=mid,low = "white",mid="orange",
                           high = "red", space = "Lab" )+
    #guides(color="none")+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
  
  return(pl)
}

Jet_slp_pred_plot<-function(model,jet_pos){
  
  test_2<-test%>%
    mutate(lat_jet_lag1=jet_pos)
  
  srd_place<-test_2%>%
    mutate(pred=predict(model,test_2)%>%as.numeric())%>%
    group_by(lon_fac,lat_fac)%>%
    select(pred)%>%
    summarise(pred=mean(pred))%>%
    rename(lon=lon_fac,lat=lat_fac)%>%
    mutate(lon=lon-360)
  
  #mid=mean(srd_place$pred)
  
  #high=max(srd_place$pred)
  
  #low=min(srd_place$pred)
  
  world<-map_data("world")
  
  pl<-ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=srd_place,aes(lon,lat,color=pred),size=4)+
    coord_sf(xlim = c(-10, 2), ylim = c(35, 60), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_color_gradient2( limits=c(-1.1,0.5),midpoint=0,low = "red",mid="orange",
                           high = "white", space = "Lab" )+
    #guides(color="none")+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
  
  return(pl)
}

#Qgam geo control plot

#RMSE plot

QGAM_geo_plot<-function(model,variable,qu,test){
  srd_place<-test%>%
    mutate(pred=predict(model,test)%>%as.numeric())%>%
    mutate(srd_error=ifelse({{variable}}-pred>0,0,1))%>%
    group_by(fac_lon,fac_lat)%>%
    select(srd_error)%>%
    summarise(srd_error=mean(srd_error))%>%
    rename(lon=fac_lon,lat=fac_lat,perc=srd_error)%>%
    mutate(lon=lon-360)
  
  mid=mean(srd_place$perc)
  
  world<-map_data("world")
  
  pl<-ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=srd_place,aes(lon,lat,color=perc),size=4)+
    coord_sf(xlim = c(-10, 6), ylim = c(35, 55), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_color_gradient2( midpoint=mid,low ="white",mid="orange",
                           high = "red", space = "Lab" )+
    labs(color="Percentage")+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
  
  return(pl)
}

#Qgam line plot

QGAM_line<-function(model,variable,qu,test){
  srd_time<-test%>%
    mutate(pred=predict(model,test)%>%as.numeric())%>%
    mutate(srd_error=ifelse({{variable}}-pred>0,0,1))%>%
    group_by(year(date))%>%
    select(srd_error)%>%
    summarise(srd_error=mean(srd_error))%>%
    rename(perc=srd_error,Year="year(date)",RMSE=srd_error)%>%
    ggplot(aes(x=Year,y=RMSE))+
    geom_line(color="black",lwd=1.5)+
    geom_hline(yintercept=qu, linetype="dashed")+
    labs(y="Proportion")+
    theme_bw()+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1))
  
  return(srd_time)
}


#QGAM line facet plot

QGAM_line_facet<-function(model_1,variable_1,model_2,variable_2,test){
  srd_time<-test%>%
    mutate(pred_1=predict(model_1,test)%>%as.numeric(),
           pred_2=predict(model_2,test)%>%as.numeric())%>%
    mutate(srd_error_1=ifelse({{variable_1}}-pred_1>0,0,1),
           srd_error_2=ifelse({{variable_2}}-pred_2>0,0,1))%>%
    group_by(year(date),lat_fac,lon_fac)%>%
    select(srd_error_1,srd_error_2)%>%
    summarise(srd_error_1=mean(srd_error_1),
              srd_error_2=mean(srd_error_2))%>%
    pivot_longer(srd_error_1:srd_error_2,names_to = "var",values_to = "RMSE")%>%
    rename(Year="year(date)")%>%
    ggplot(aes(x=Year,y=RMSE,color=var))+
    geom_line(lwd=1.5)+
    ylim(0,1)+
    facet_wrap(vars(lat_fac,lon_fac))+
    theme_bw()
}


#Compare gams and qgams' smooths

draw_leo<-function(dat.gam,dat.qgam,variab){
  
  dat_gam<-smooth_estimates(dat.gam)%>%
    filter(smooth==s({{variab}}))
  
  dat_qgam<-smooth_estimates(dat.qgam)%>%
    filter(smooth==s({{variab}}))
  
  ggplot()+
    geom_ribbon(aes(x={{variab}},y=est,ymin=est-se,ymax=est+se),data=dat_gam,fill="steelblue",alpha=0.3)+
    geom_line(aes(x={{variab}},y=est),data=dat_gam,col="black",linetype="dashed")+
    geom_ribbon(aes(x={{variab}},y=est,ymin=est-se,ymax=est+se),data=dat_qgam,fill="steelblue",alpha=0.7)+
    geom_line(aes(x={{variab}},y=est),data=dat_qgam, col="dark red",lwd=1)+
    labs(y="Effect")+
    theme_bw()
  
}


####Quantile Models

Quant_eval<-function(predicted,variable,qu,test){
  srd_place<-test%>%
    mutate(pred=predicted%>%as.numeric())%>%
    mutate(srd_error=ifelse({{variable}}-pred>0,0,1))%>%
    group_by(fac_lon,fac_lat)%>%
    select(srd_error)%>%
    summarise(srd_error=mean(srd_error,na.rm=TRUE))%>%
    rename(lon=fac_lon,lat=fac_lat,perc=srd_error)
}


overpred<-function(dat,qu,direc=1,low_lim,up_lim,midpo){
  
  dat<-dat%>%
    mutate(lon=ifelse(as.numeric(lon)>300,as.numeric(lon)-360,lon))%>%
    mutate(perc=abs(perc-qu))
  
  
  world<-map_data("world")
  
  ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=dat,aes(lon,lat,fill=perc),size=4,shape=21, stroke=NA)+
    coord_sf(xlim = c(-10, 1), ylim = c(35, 52.5), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank(),legend.title = element_text(size = 12))+
    #scale_color_viridis_c(limits=c(0,(1-qu)),oob=squish,option=opt,direction=direc,labels = scales::number_format(accuracy = 0.0001,
    #                                                                                                              decimal.mark = '.'))+
    scale_fill_gradient2(limits=c(low_lim,up_lim),oob=squish,
    low = muted("green"),
  mid = "#F6F3E7",
  high = muted("red"),
  midpoint=midpo)+
    labs(fill=expression(widehat(Bias)))+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))+
    guides(color = guide_colorbar(reverse = TRUE))
}

#expression(abs(hat(tau)-tau))

overpred_2<-function(dat,qu,opt="B",direc=1){
  
  dat<-dat%>%
    mutate(lon=ifelse(as.numeric(lon)>300,as.numeric(lon)-360,lon))%>%
    mutate(perc=abs(perc-qu))
  
  
  world<-map_data("world")
  
  ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "light gray", fill = "light gray", size = 0.5, alpha=0.5
    )+
    geom_tile(data=dat,aes(lon,lat,fill=perc))+
    stat_contour(aes(fill=..level..), geom="polygon") + 
    geom_contour(color="white", alpha=0.5)+
    #geom_path(data = world, aes(x = long, y = lat, group = group))+
    coord_sf(xlim = c(-10, 1), ylim = c(35, 52.5), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank(),legend.title = element_text(size = 12))+
    scale_fill_viridis_c(limits=c(0,(1-qu)),oob=squish,option=opt,direction=direc,labels = scales::number_format(accuracy = 0.0001,
                                                                                                                  decimal.mark = '.'))+
    labs(fill=expression(widehat(Bias)))+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))+
    guides(fill = guide_colorbar(reverse = TRUE))
}


rho <- function(u,tau){pmax(tau * u, (tau - 1) * u)}


Pseudo_R2_tot<-function(predicted,variable,qu,train,test){
  
  
  base<-train%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(qu_var=quantile({{variable}},qu,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    dplyr::select(-{{variable}})%>%
    ungroup()
  
  test<-test%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    left_join(base)%>%
    mutate(a={{variable}}-qu_var,b={{variable}}-predicted%>%as.numeric())%>%
    #group_by(fac_lon,fac_lat)%>%
    #summarise(a=mean(a_error),b=mean(b_error))%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    rename(lon=fac_lon,lat=fac_lat)%>%
    dplyr::select(R2)%>%
    pull()%>%
    .[1]
  
  return(test)
  
  
}


Pseudo_R2_point<-function(predicted,variable,qu,train,test){
  
  
  base<-train%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(qu_var=quantile({{variable}},qu,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    dplyr::select(-{{variable}})%>%
    ungroup()
  
  test<-test%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    left_join(base)%>%
    mutate(a={{variable}}-qu_var,b={{variable}}-predicted%>%as.vector())%>%
    group_by(fac_lon,fac_lat)%>%
    #summarise(a=mean(a_error),b=mean(b_error))%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    rename(lon=fac_lon,lat=fac_lat)
  
  return(test)
  
  
}


Pseudo_R2_point_EVT<-function(model_evt,variable,qu,train,test){
  
  
  base<-train%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(qu_var=quantile({{variable}},qu,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    dplyr::select(-{{variable}})%>%
    ungroup()
  
  test<-test%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    left_join(base)%>%
    mutate(a={{variable}}-qu_var,b={{variable}}-predict(model_evt,test,prob=qu))%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    rename(lon=fac_lon,lat=fac_lat)
  
  return(test)
  
  
}

Pseudo_R2_tot_multi<-function(predicted_1,predicted_2,variable,qu,train,test){
  
  a<-rho(test%>%dplyr::select({{variable}})%>%pull()-predicted_2%>%as.vector(),qu)
  b<-rho(test%>%dplyr::select({{variable}})%>%pull()-predicted_1%>%as.vector(),qu)
  
  R2<-1-sum(b,na.rm=TRUE)/sum(a,na.rm=TRUE)
  
  return(R2)
  
  
}

Pseudo_R2_tot_multi_EVT<-function(predicted_1,model_evt,variable,qu,train,test){
  
  
  "test<-test%>%
    mutate(logscale=predict(model_evt,test)$logscale,shape=predict(model_evt,test)$shape)%>%
    mutate(scale=exp(logscale))%>%
    drop_na(scale,shape)"
  
  R2<-test%>%
    mutate(a={{variable}}-predict(model_evt,test,prob=qu),b={{variable}}-predicted_1%>%as.vector())%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    pull(R2)
  
  '
  test<-test%>%
    mutate(logscale=predict(model_evt,test)$logscale,shape=predict(model_evt,test)$shape)%>%
    mutate(scale=exp(logscale))
  
  
  b<-rho(predict(model_1,test)%>%as.vector()-test%>%dplyr::select({{variable}})%>%pull(),qu)
  a<-rho(qgpd(qu_evt,sigmau = test$scale,xi=test$shape)+da$threshold%>%as.vector()-test%>%dplyr::select({{variable}})%>%pull())
  
  R2<-1-sum(b)/sum(a)'
  
  R2<-R2[1]
  
  return(R2)
  
  
}

Pseudo_R2_point_multi_EVT<-function(predicted_1,model_evt,variable,qu,train,test){
  
  
  
  
  #test<-test%>%
    #mutate(logscale=predict(model_evt,test)$logscale,shape=predict(model_evt,test)$shape)%>%
    #mutate(scale=exp(logscale))%>%
    #drop_na(scale,shape)
  
  test<-test%>%
    mutate(a={{variable}}-predict(model_evt,test,type="quantile",prob=qu),b={{variable}}-predicted_1%>%as.vector())%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    dplyr::select(fac_lat,fac_lon,b,a,b_rho,a_rho,R2)%>%
    rename(lon=fac_lon,lat=fac_lat)
  
 # R2<-goodfit((test$b),(test$a),qu)
  
  return(test)
  
  
}

Pseudo_R2_point_multi<-function(predicted_1,predicted_2,variable,qu,train,test){
  
  
  base<-train%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(qu_var=quantile({{variable}},qu,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    dplyr::select(-{{variable}})%>%
    ungroup()
  
  test<-test%>%
    dplyr::select(fac_lon,fac_lat,{{variable}})%>%
    left_join(base)%>%
    mutate(a={{variable}}-predicted_2%>%as.vector(),b={{variable}}-predicted_1%>%as.vector())%>%
    group_by(fac_lon,fac_lat)%>%
    #summarise(a=mean(a_error),b=mean(b_error))%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    rename(lon=fac_lon,lat=fac_lat)
  
  #R2<-goodfit((test$b),(test$a),qu)
  
  return(test)
  
  
}





Plot_R2<-function(dat,qu,opt="B",direc=1,low_lim=0,up_lim=0.75){
  
  dat<-dat%>%
    mutate(lon=ifelse(as.numeric(lon)>300,as.numeric(lon)-360,lon))%>%
    mutate(R2=as.numeric(R2))%>%
    mutate(pos=ifelse(R2>0,1,0)%>%as.character())
  #mutate(R2=if_else(is.na(R2),0,R2))
  
  world<-map_data("world")
  
  ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=dat,aes(lon,lat,fill=R2),size=4,shape=21,stroke=NA)+
    #geom_point(data=dat%>%filter(pos==0),aes(lon,lat),size=3.5,shape=4)+
    coord_sf(xlim = c(-10, 1), ylim = c(35, 52.5), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank(),legend.title = element_text(size = 12))+
    #scale_fill_viridis_c(limits=c(low_lim,up_lim),oob=squish,option=opt,direction=direc,labels = scales::number_format(accuracy = 0.0001,
    #                                                                                                                    decimal.mark = "."))+
    scale_fill_gradient2(limits=c(low_lim,up_lim),oob=squish,
                         low = muted("red"),
                         mid = "#F6F3E7",
                         high = muted("green"))+
    labs(fill=expression(paste("Pseudo ", R^{2})))+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
    #scale_shape_manual(values = c(4))+
    #guides(shape="none")
  
}


plot_smooth<-function(mod,col="dark red"){
  a<-summary(mod)
  
  label<-paste("P-value:",as.character(format(round(a$s.table[,"p-value"], 3), nsmall = 3)%>%.["s(US_temp_lag2)"]%>%unlist()))
  
  mod%>%
    smooth_estimates(smooth="s(US_temp_lag2)")%>%
    add_confint()%>%
    ggplot(aes(y = est, x = US_temp_lag2)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                alpha = 0.2, fill = "steelblue")+ 
    geom_line(colour = col, size = 1.5) +
    #labs(y = "Partial effect",
    #title = expression("Partial effect of" ~ f(x[2])),
    #x = expression(x[2]))+
    annotate("text",x=1,y=0.9, label=label,size=4)+
    theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))
}


sanity_gg<-function(dat,var,longitude,latitude,xlim_west,xlim_east,y_lim_south,y_lim_north){
  
  dat<-dat%>%
    mutate(lon=ifelse(as.numeric(lon)>180,as.numeric(lon)-360,lon))
  
  world<-map_data("world")
  
  ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=dat,aes({{longitude}},{{latitude}},fill={{var}}),size=2.5,shape=22,stroke=NA)+
    coord_sf(xlim = c(xlim_west, xlim_east), ylim = c(y_lim_south, y_lim_north), expand = FALSE)+
    theme_bw()+
    theme(legend.title = element_text(size = 12),plot.title = element_text(size=15,face = "bold"))+
    guides(color="none")+
    scale_fill_gradient(high="light blue",low="#132B43")
}

