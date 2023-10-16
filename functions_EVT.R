density_comp<-function(out_var,ext_var,dat_out,dat_ext){
  ggplot()+
    geom_density(aes(x={{out_var}},y=..scaled..),dat=dat_out,fill="light blue",alpha=1.2,size=0.8)+
    geom_density(aes(x={{ext_var}},y=..scaled..),dat=dat_ext,fill="red",color="black",alpha=0.4,size=0.8)+
    theme_minimal()
}



############################COLD SPELL COMPOSITES#################
POT_cold_spell<-function(out_dat=Y,in_var,qu,over=TRUE){
  
  dat<-left_join(out_dat,X,by="date")
  
  thresh<-dat%>%
    summarise(thresh=quantile({{in_var}},qu,na.rm=TRUE))%>%
    pull(thresh)
  
 if(over==TRUE) 

{dat<-dat%>%
  rename(key_var={{in_var}})%>%
  filter(key_var>=thresh)%>%
  #group_by(date)%>%
  #filter(key_var==max(key_var))%>%
  #ungroup()%>%
  group_by(lat,lon)%>%
  arrange(date)%>%
  mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
  group_by(id2,lat,lon)%>%
  mutate(row=row_number())%>%
  filter(row==1)%>%
  ungroup()%>%
  select(-id2,-row)
 }else{
   dat<-dat%>%
     rename(key_var={{in_var}})%>%
     filter(key_var<=thresh)%>%
     #group_by(date)%>%
     #filter(key_var==min(key_var))%>%
     #ungroup()%>%  group_by(lat,lon)%>%
     group_by(lat,lon)%>%
     arrange(date)%>%
     mutate(id2 = cumsum(c(T, diff(date) > 4)))%>%
     group_by(id2,lat,lon)%>%
     mutate(row=row_number())%>%
     filter(row==1)%>%
     ungroup()%>%
     select(-id2,-row)
 }

pl_wind<-dat%>%
  #distinct(date,.keep_all = TRUE)%>%
  select(starts_with("wind"))%>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
  pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = "wind")%>%
  mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
  mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
  select(-va,-type)%>%
  ggplot()+
  geom_line(aes(lag,wind),color="dark blue",size=1.5)+
  theme_minimal()+
  ylim(-1.2,1.2)+
  geom_hline(yintercept = 0)

pl_temp<-dat%>%
  #distinct(date,.keep_all = TRUE)%>%
  select(starts_with("temp"))%>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
  pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = "temp")%>%
  mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
  mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
  select(-va,-type)%>%
  ggplot()+
  geom_line(aes(lag,temp),color="dark blue",size=1.5)+
  theme_minimal()+
  ylim(-1.2,1.2)+
  geom_hline(yintercept = 0)

pl_prec<-dat%>%
  #distinct(date,.keep_all = TRUE)%>%
  select(starts_with("prec"))%>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
  pivot_longer(cols=everything(),names_sep="_", names_to = c("va","type"), values_to = "prec")%>%
  mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
  mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
  select(-va,-type)%>%
  ggplot()+
  geom_line(aes(lag,prec),color="dark blue",size=1.5)+
  theme_minimal()+
  ylim(-1.2,1.2)+
  geom_hline(yintercept = 0)

pl_jet_strength<-dat%>%
  #distinct(date,.keep_all = TRUE)%>%
  select(starts_with("jet_strength"))%>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))%>%
  pivot_longer(cols=everything(),names_sep="_strength_", names_to = c("va","type"), values_to = "jet_strength")%>%
  mutate(lag=ifelse(str_detect(type,"lag")==TRUE,parse_number(type)*(-1),parse_number(type)))%>%
  mutate(lag=ifelse(is.na(lag)==TRUE,0,lag))%>%
  select(-va,-type)%>%
  ggplot()+
  geom_line(aes(lag,jet_strength),color="dark blue",size=1.5)+
  theme_minimal()+
  ylim(-1.2,1.2)+
  geom_hline(yintercept = 0)
  
l=list(pl_wind,pl_temp,pl_prec,pl_jet_strength,threshold=dat$thresh[1],dat,thresh)
return(l)
}


pick_and_mean<-function(dat,var){
  
  #var_quot=enquo(var)%>%rlang::as_label()
  
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


#lat - latitude of the point in Western Europe
#lon - longitude of the point in Western Europe
#prec - precipitation levels, anomaly in mm per day at the surface
#wind - wind speed anomaly at 10m height, meters per second
#temp - temperature anomaly at 10m height, kelvin 


#tot<-left_join(out,dat)%>%
#select(date,lat,lon,US_temp,wind,prec,temp)



# These are some of the functions I use to perform the monte carlo perform test. 
# You can use this function, or build some of your own . 
# If you use those functions it is good if you anyway try to understand how they work, 
# so that you understand the algorithm you follow to perform the test. 
# In addition, I wrote them quite fast - with high probability they contain at least a couple of small mistakes! :)

under_pick_and_mean_fix_monte<-function(dat_tot,var_in,var_out,len){
  
  #Shuffles the variable of interest 
  das<-dat_tot%>%
    select(date,{{var_in}})%>%
    slice_sample(n=len)#it is possible that a smaller sample could also work, feel free to experiment! 
  
  
  dat_tot%>%
    select(date,{{var_out}},lon,lat)%>%
    right_join(das)%>%
    summarise(mean=mean({{var_out}},na.rm=TRUE))%>%
    pull(mean)
  
  
}

#under_pick_and_mean_fix_monte(tot,US_temp_lag2,wind,186)



POT_cold_spell_monte_mean<-function(dat,in_var,out_var,qu=0.05,B=5,seed=3){
  
  # Data cleaning/preparation for following functions
  #dat<-left_join(out_dat,in_dat)
  
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
  
  #set.seed(seed)
  #Takes B samples to empirically determine 5% significance levels around H0. 
  #Based on the above function
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

#POT_cold_spell_monte_mean(tot,US_temp,wind,B=50)


###################QUANTILE#################



under_pick_and_qu_fix_monte<-function(dat_tot,var_in,var_out,len,qu_out=qu_out){
  
  #Shuffles the variable of interest 
  das<-dat_tot%>%
    select(date,{{var_in}})%>%
    slice_sample(n=len)#it is possible that a smaller sample could also work, feel free to experiment! 
  
  
  dat_tot%>%
    select(date,{{var_out}},lon,lat)%>%
    right_join(das)%>%
    summarise(qua=quantile({{var_out}},qu_out,na.rm=TRUE))%>%
    pull(qua)
  
  
}

#under_pick_and_qu_fix_monte(tot,US_temp_lag2,wind,186,0.95)


POT_cold_spell_monte_qu<-function(dat,in_var,out_var,qu_in=0.05,qu_out=0.95,B=5,seed=3){
  
  # Data cleaning/preparation for following functions
  #dat<-left_join(out_dat,in_dat)
  
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
  
  #set.seed(seed)
  #Takes B samples to empirically determine 5% significance levels around H0. 
  #Based on the above function
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
    #mutate(big_jet=ifelse(jet_strength_lag1>thresh_jet,1,0))%>%
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
      


POT_declust_jet<-function(dat,out_var,qu,over=TRUE,span=length(year_samp)){
  
  
  thresh<-dat%>%
    summarise(thresh=quantile({{out_var}},qu,na.rm=TRUE))%>%
    pull(thresh)
  
  thresh_jet<-dat%>%
    summarise(thresh_jet=quantile(jet_strength_lag1,qu,na.rm=TRUE))%>%
    pull(thresh_jet)
  
  
  dat<-dat%>%
    rename(key_var={{out_var}})%>%
    filter(key_var>=thresh)%>%
    mutate(big_jet=ifelse(jet_strength_lag1>thresh_jet,1,0))%>%
    mutate(date=as.Date(date))
  
  
  i<-dat%>%
       distinct(date,.keep_all = TRUE)%>%
       dplyr::select(date,big_jet)%>%
       filter(big_jet==1)%>%
       arrange(date)%>%
       mutate(id2 = cumsum(c(T, diff(date) <= 4)))%>%
       group_by(id2)%>%
       mutate(row=row_number())%>%
       filter(row==1)%>%
       ungroup()%>%
       dplyr::select(-id2,-row)%>%
       pull(date)
     
     dati<-dat%>%
       filter(!(date%in%i))
     
     a<-list(dati=dati,threshold=thresh,jet_thresh=thresh_jet)
     return(a)
     
}
  
  #if(over==TRUE){
    #e<-fevd(dat$key_var,threshold = thresh,type="GP",span=span,na.action = na.exclude)
    #plot(e)
  #}else{  
    #e<-fevd(-dat$key_var,threshold = -thresh,type="GP",span=span,na.action = na.exclude)
    #plot(e)}


