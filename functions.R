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
  
  R2<-test%>%
    mutate(a={{variable}}-predict(model_evt,test,prob=qu),b={{variable}}-predicted_1%>%as.vector())%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    pull(R2)
  
  R2<-R2[1]
  
  return(R2)
  
  
}

Pseudo_R2_point_multi_EVT<-function(predicted_1,model_evt,variable,qu,train,test){
  
  test<-test%>%
    mutate(a={{variable}}-predict(model_evt,test,type="quantile",prob=qu),b={{variable}}-predicted_1%>%as.vector())%>%
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    group_by(fac_lon,fac_lat)%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    dplyr::select(fac_lat,fac_lon,b,a,b_rho,a_rho,R2)%>%
    rename(lon=fac_lon,lat=fac_lat)
  
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
    mutate(a_rho=rho(a,tau=qu),b_rho=rho(b,tau=qu))%>%
    mutate(R2=1-sum(b_rho,na.rm=TRUE)/sum(a_rho,na.rm=TRUE))%>%
    distinct(fac_lon,fac_lat,.keep_all = TRUE)%>%
    ungroup()%>%
    rename(lon=fac_lon,lat=fac_lat)
  
  return(test)
  
  
}


Plot_R2<-function(dat,qu,opt="B",direc=1,low_lim=0,up_lim=0.75){
  
  dat<-dat%>%
    mutate(lon=ifelse(as.numeric(lon)>300,as.numeric(lon)-360,lon))%>%
    mutate(R2=as.numeric(R2))%>%
    mutate(pos=ifelse(R2>0,1,0)%>%as.character())
  
  world<-map_data("world")
  
  ggplot()+
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgray", size = 0.5
    )+
    geom_point(data=dat,aes(lon,lat,fill=R2),size=4,shape=21,stroke=NA)+
    coord_sf(xlim = c(-10, 1), ylim = c(35, 52.5), expand = FALSE)+
    theme_bw()+
    theme(axis.title = element_blank(),legend.title = element_text(size = 12))+
    scale_fill_gradient2(limits=c(low_lim,up_lim),oob=squish,
                         low = muted("red"),
                         mid = "#F6F3E7",
                         high = muted("green"))+
    labs(fill=expression(paste("Pseudo ", R^{2})))+
    scale_y_continuous(breaks = c(seq(30,60,by=2.5)))
  
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
    annotate("text",x=1,y=0.9, label=label,size=4)+
    theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))
}



fig_1<-function(dat,var,longitude,latitude,xlim_west,xlim_east,y_lim_south,y_lim_north){
  
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



pick_and_mean<-function(dat,var){
  
  res<-dat%>%
    sample_n(dim(dat)[1],replace=TRUE)%>%
    select(starts_with(var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  return(res[,1])
}

under_pick_and_mean_and_quant<-function(dat,var,var_in,qu){
  
  
  res<-dat%>%
    sample_n(dim(dat)[1],replace=TRUE)%>%
    mutate(thresh=quantile({{var_in}},qu,na.rm=TRUE))%>%
    rename(key_var={{var_in}})%>%
    filter(key_var<=thresh)%>%
    mutate(date=as.Date.character(date))%>%
    arrange(date)%>%
    group_by(lat,lon)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(id2,lat,lon)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(starts_with(var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  return(res[,1])
}

over_pick_and_mean_and_quant<-function(dat,var,var_in,qu){
  
  res<-dat%>%
    sample_n(dim(dat)[1],replace=TRUE)%>%
    mutate(thresh=quantile({{var_in}},qu,na.rm=TRUE))%>%
    rename(key_var={{var_in}})%>%
    filter(key_var>=thresh)%>%
    group_by(date)%>%
    filter(key_var==max(key_var))%>%
    ungroup()%>%
    group_by(lat,lon)%>%
    arrange(date)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(id2,lat,lon)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(starts_with(var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  return(res[,1])
}

POT_pick_and_dec<-function(data,in_var,declust_n=4,qu=0.05){ 
  
  data%>%
    rename(key_var={{in_var}})%>%
    mutate(thresh=quantile(key_var,qu,na.rm=TRUE))%>%
    group_by(date)%>%
    filter(key_var==min(key_var,na.rm=TRUE))%>%
    ungroup()%>%
    filter(key_var<=thresh)%>%
    mutate(date=as.Date(date))%>%
    arrange(date)%>%
    mutate(id2 = cumsum(c(T, diff(date) > declust_n)))%>%
    group_by(id2)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(-id2,-row)
}

under_pick_and_mean_fix_th<-function(dat,var,var_in,qu,thresho){
  
  
  res<-dat%>%
    sample_n(dim(dat)[1],replace=TRUE)%>%
    rename(key_var={{var_in}})%>%
    filter(key_var<=thresho)%>%
    mutate(date=as.Date.character(date))%>%
    arrange(date)%>%
    group_by(lat,lon)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(id2,lat,lon)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(starts_with(var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  return(res[,1])
}

over_pick_and_mean_fix_th<-function(dat,var,var_in,qu,thresho){
  
  res<-dat%>%
    sample_n(dim(dat)[1],replace=TRUE)%>%
    rename(key_var={{var_in}})%>%
    filter(key_var>=thresho)%>%
    group_by(date)%>%
    filter(key_var==max(key_var))%>%
    ungroup()%>%
    group_by(lat,lon)%>%
    arrange(date)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(id2,lat,lon)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(starts_with(var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  return(res[,1])
}

##################PERMUTATION TESTS###############

under_pick_and_mean_fix_monte<-function(dat_tot,var_in,var_out,len){
  
  #Shuffles the variable of interest 
  das<-dat_tot%>%
    select(date,{{var_in}})%>%
    slice_sample(n=len)
  
  
  dat_tot%>%
    select(date,{{var_out}},lon,lat)%>%
    right_join(das)%>%
    summarise(mean=mean({{var_out}},na.rm=TRUE))%>%
    pull(mean)
}

POT_cold_spell_monte_mean<-function(dat,in_var,out_var,qu=0.05,B=5,seed=3){
  
  quote_var=enquo(out_var)%>%rlang::as_label()
  dat<-dat%>%
    mutate(var_in={{in_var}},var_out={{out_var}})
  qu<-qu
  
  dat_1<-dat%>%
    rename(key_var={{in_var}})%>%
    mutate(thresh=quantile(key_var,qu,na.rm=TRUE))%>%
    filter(key_var<=thresh)%>%
    mutate(date=as.Date(date))%>%
    group_by(lat,lon)%>%
    arrange(date)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(lat,lon,id2)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(-id2,-row)%>%
    mutate(date=as.Date(date))
  
  len<-nrow(dat_1)/
    dat%>%
    distinct(lat,lon)%>%
    nrow()
  c<-replicate(B,under_pick_and_mean_fix_monte(dat=dat, var_in=var_in,
                                               var_out=var_out,len=len))
  
  pl_dat<-dat_1%>%
    select(starts_with(quote_var))%>%
    summarise(across(everything(), ~ mean(.x, na.rm=TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = quote_var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  
  high<-quantile(c,0.975,na.rm=TRUE)
  
  low<-quantile(c,0.025,na.rm=TRUE)
  
  pl<-ggplot(dat=pl_dat)+
    geom_line(aes(lag,{{out_var}}),color="blue",size=1.5)+
    geom_hline(yintercept = high,color="black",linetype="dashed")+
    geom_hline(yintercept=low,color="black",linetype="dashed")+
    geom_hline(yintercept=0)+
    theme_minimal()
  
  #returns the plot, the original data and the summary of the resampling
  l=list(dat=dat_1,pl_dat=pl_dat,pl=pl,high=high,low=low,len=len)
  
  return(l)
}


###################QUANTILE#################

under_pick_and_qu_fix_monte<-function(dat_tot,var_in,var_out,len,qu_out=qu_out){
  
  #Shuffles the variable of interest 
  das<-dat_tot%>%
    select(date,{{var_in}})%>%
    slice_sample(n=len)
  
  
  dat_tot%>%
    select(date,{{var_out}},lon,lat)%>%
    right_join(das)%>%
    summarise(qua=quantile({{var_out}},qu_out,na.rm=TRUE))%>%
    pull(qua)
}

POT_cold_spell_monte_qu<-function(dat,in_var,out_var,qu_in=0.05,qu_out=0.95,B=5,seed=3){
  
  
  quote_var=enquo(out_var)%>%rlang::as_label()
  dat<-dat%>%
    mutate(var_in={{in_var}},var_out={{out_var}})
  qu_in<-qu_in
  qu_out<-qu_out
  
  dat_1<-dat%>%
    rename(key_var={{in_var}})%>%
    mutate(thresh=quantile(key_var,qu_in,na.rm=TRUE))%>%
    filter(key_var<=thresh)%>%
    mutate(date=as.Date(date))%>%
    group_by(lat,lon)%>%
    arrange(date)%>%
    mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
    group_by(lat,lon,id2)%>%
    mutate(row=row_number())%>%
    filter(row==1)%>%
    ungroup()%>%
    select(-id2,-row)%>%
    mutate(date=as.Date(date))
  
  len<-nrow(dat_1)/
    dat%>%
    distinct(lat,lon)%>%
    nrow()
  
  c<-replicate(B,under_pick_and_qu_fix_monte(dat=dat, var_in=var_in,
                                             var_out=var_out,len=len,qu_out=qu_out))
  
  pl_dat<-dat_1%>%
    select(starts_with(quote_var))%>%
    summarise(across(everything(), ~ quantile(.x,qu_out, na.rm=TRUE)))%>%
    pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = quote_var)%>%
    mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
    mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
    select(-va,-type)
  
  high<-quantile(c,0.975)
  
  low<-quantile(c,0.025)
  
  center<-dat%>%
    pull({{out_var}})%>%
    quantile(qu_out,na.rm=TRUE)
  
  pl<-ggplot(dat=pl_dat)+
    geom_line(aes(lag,{{out_var}}),color="blue",size=1.5)+
    geom_hline(yintercept = high,color="black",linetype="dashed")+
    geom_hline(yintercept=low,color="black",linetype="dashed")+
    geom_hline(yintercept=center)+
    theme_minimal()
  
  #returns the plot, the original data and the summary of the resampling
  l=list(dat=dat_1,pl_dat=pl_dat,pl=pl,high=high,low=low,len=len)
  
  return(l)
}


POT_declust_out<-function(dat,out_var,qu,over=TRUE,span=length(year_samp),declust=FALSE){
  
  
  thresh<-dat%>%
    summarise(thresh=quantile({{out_var}},qu,na.rm=TRUE))%>%
    pull(thresh)
  
  
  dat<-dat%>%
    rename(key_var={{out_var}})%>%
    filter(key_var>=thresh)%>%
    mutate(date=as.Date(date))
  
  if(declust==TRUE){
    
    i<-dat%>%
      distinct(date,.keep_all = TRUE)%>%
      dplyr::select(date)%>%
      mutate(date=as.Date(date))%>%
      #filter(big_jet==1)%>%
      arrange(date)%>%
      mutate(id2 = cumsum(c(T, diff(date) <= 1)))%>%
      group_by(id2)%>%
      mutate(row=row_number())%>%
      filter(row==1)%>%
      ungroup()%>%
      dplyr::select(-id2,-row)%>%
      pull(date)
    
    dat<-dat%>%
      filter(!(date%in%i))
  }else{}
  
  a<-list(dat=dat,threshold=thresh)
  return(a)
}

