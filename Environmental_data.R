rm(list=ls())


library(ggplot2);library(tidyr);library(reshape2);library(dplyr); library(mgcv);library(lubridate);library(magrittr);
library(knitr);library(scales); library(readxl);library(ggpubr)

options(mc.cores = parallel::detectCores())##optimise model run

##Load Summarised data
mywd<-"C:/Users/ccastros/OneDrive - Australian Institute of Marine Science/Documents/2023/Light"
myplots=paste(mywd,"/Plots", sep="")

load(paste(mywd,"/RCode-Github/Environmental.RData", sep=""))


sitecolor=c("#CC6666", "#9999CC")
se <- function(x) sd(x) / sqrt(length(x))
 
##------------------
## MS9-loggers PLOTS
##------------------

#PAR

p1<-DailyMean%>%
  ggplot(aes(x=as.POSIXct(Date), y=meanPAR,group=Site, color=Site, fill=Site))+
  geom_line(aes(y = meanPAR))+
  scale_color_manual(values=sitecolor)+
  geom_ribbon(aes(y=meanPAR, ymin = meanPAR - sePAR, ymax = meanPAR + sePAR), alpha=0.3, colour="NA") +
  scale_fill_manual(values=sitecolor, aesthetics=c("fill"))+facet_wrap(vars(Season),scales = "free_x")+
  scale_x_datetime(date_breaks = "3 weeks",date_labels = "%b %d")+
  theme_bw()+xlab("Date")+
  ylab("Daily averaged PAR\n (µmol m\u207B\u00b2 s\u207B\u00b9)")+
theme(legend.position = "top",text = element_text(size=10))


#Tides
p2<-Light_mean%>%filter(!Season == "Spring")%>%ggplot(aes(x=as.POSIXct(DateTime), y=DEPTH, colour=Site))+
  geom_line(alpha=0.5)+xlab("")+#ggtitle("Depth")+
  facet_wrap(vars(Season),scales = "free_x")+xlab("Date")+ylab("Depth (m)")+
  scale_x_datetime(date_breaks = "3 week",date_labels = "%b %d")+
  scale_color_manual(values=sitecolor)+theme_bw()+
  theme(legend.position = "top",text = element_text(size=10))
ggsave(file=paste(mywd,myplots,"Depth.jpg",sep = "/"), dpi=300, units = c("in"),width = 5,height = 3)

MMM=27##Regional based on NOAA Coral Reef Watch
summerMM=MMM#Light_mean$MM[Light_mean$Month=="March"]
WinterMM=22.22#Light_mean$MM[Light_mean$Month=="August"] 

# text to annotate graphs
graphLabels <- data.frame(Season = c("Winter","Winter","Summer"),Date=rep(as.POSIXct("2022-07-11"),3),
                          ref = c("MMM","SM",""), Site= c("Bundegi", "Tantabiddi","Tantabiddi"))
graphLabels$Season=factor(graphLabels$Season, levels=c("Winter","Summer"))

seasonalref=Light_mean%>%select(Site,Season,DateTime,MM)%>%filter(!Season=="Spring")%>%droplevels()
seasonalref$MM[seasonalref$Season=="Summer"]=summerMM
seasonalref$MM[seasonalref$Season=="Winter"]=WinterMM
seasonalref$Season=factor(seasonalref$Season, levels = c("Winter","Summer"))

##--------------------------------------
##predict daily PAR curves 
##--------------------------------------

#daily predictions per second
daily_preds=list()

site=unique(Light_mean$Site)

for (i in 1:length(site)){
  d_data=Light%>%filter(Site %in% c(site[i]))%>%droplevels()
  date=unique(d_data$Date)
  for (d in 1:length(date)){
    
    D=d_data%>%filter(Date==date[d])%>%select("calcPAR","Time","DateTime")%>%droplevels()
    Dmodel<-gam(calcPAR~s(Time), data=D)
    timesequence=seq(from=D$DateTime[1], to=D$DateTime[length(D$DateTime)], by="sec")
    timeseq=format(as.POSIXct(timesequence), format = "%H%M%S") 
    
    D_newd=data.frame("Time"=as.numeric(timeseq),"DateTime"=timesequence)
    pred_PAR<-predict.gam(Dmodel,D_newd)
    D_newd$PAR=pred_PAR
    
    ##---------------
    ##Convert to DLI
    ##---------------
    
    #SUM instantaneous PAR over 12h period  (i.e. predicted PAR umol.m2.s)/ 1million to convert to mol
    D_newd$PAR[D_newd$PAR<0]=0## no negative PAR values or remove 
    D_newd$DLI=D_newd$PAR/1000000 ##convert to mol
    D_newd$Date=date[d]
    D_newd$Site=site[i]
    
    daily_preds=c(daily_preds,list(D_newd))
  }
}

MS9.preds=do.call("rbind",daily_preds)
rm(D, Dmodel,d_data,date,timesequence,timeseq,D_newd,pred_PAR,daily_preds)

##
MS9.preds$Month=months(MS9.preds$DateTime)
MS9.preds$Season[MS9.preds$Month=="June"|MS9.preds$Month=="July"|MS9.preds$Month=="August"]<-"Winter"
MS9.preds$Season[MS9.preds$Month=="September"|MS9.preds$Month=="October"|MS9.preds$Month=="November"]<-"Spring"##Incomplete. Leave out for now
MS9.preds$Season[MS9.preds$Month=="December"|MS9.preds$Month=="January"|MS9.preds$Month=="February"]<-"Summer"
MS9.preds$Season=factor(MS9.preds$Season, levels = c("Winter","Summer","Spring"))

#Extract time and add identical date for plotting 
MS9.preds$time=strftime(MS9.preds$DateTime, format="%H:%M:%S")
MS9.preds$time=as.POSIXct(MS9.preds$time, format="%H:%M:%S")

y_lab=("Daily PAR cycle\n (µmol m\u207B\u00b2 s\u207B\u00b9)")

wsplot=MS9.preds%>%
  filter(!Season == "Spring" & PAR>0)%>%
  ggplot(aes(x=time,y=PAR,colour=Site))+
  facet_wrap(vars(Season))+
  ylim(c(0,1500))+geom_smooth()+
  scale_color_manual(values = sitecolor)+
  stat_summary(fun="min", geom="line",linetype="dotted")+
  stat_summary(fun="max", geom="line", linetype="dotted")+
  theme_bw()+labs(x="Time (hours)", y=y_lab)+
  theme(legend.position = "top",text = element_text(size=10))


##----------------------------------------
##Use Temp loggers for daily variability
##----------------------------------------

p3<-Tempdata%>%filter(!Season=="Spring")%>%#mutate(Season=factor(Season, levels=c("Winter","Summer")))%>%
  ggplot(aes(x=as.POSIXct(time), y=Temp, colour=Site, fill=Site))+
  geom_line(aes(y=Temp), alpha=0.6)+
  scale_fill_manual(values=sitecolor, aesthetics=c("fill"))+facet_wrap(vars(Season),scales = "free_x")+
  geom_hline(aes(yintercept = MMM), color="black",linetype="dashed")+
  geom_line(data=seasonalref,aes(x=as.POSIXct(Date), y= MM), colour="magenta", linetype="dotted",size=0.5)+
  geom_text(data=subset(graphLabels,ref=="MMM"),aes(y=MMM+0.6, x=as.POSIXct(Date),label=ref),color="black",size=3)+
  geom_text(data=subset(graphLabels,ref=="SM"),aes(y=WinterMM+1, x=as.POSIXct(Date),label=ref),color="magenta",size=3)+
  scale_color_manual(values=sitecolor)+xlab("Date")+
  facet_wrap(vars(Season),scales = "free_x")+
  scale_x_datetime(date_breaks = "3 week",date_labels = "%b %d")+
  theme_bw()+xlab("Date")+ylab("Temperature (C\u00b0)\n at depth ")+
  theme(legend.position = "top",text = element_text(size=10))


Temp_range=Tempdata%>%
  group_by(Site,Month,Season,factor,Date)%>%
  summarise(meanTemp=mean(Temp), minTemp=min(Temp), maxTemp=max(Temp))

Temp_range%<>%mutate(deltaTemp=maxTemp-minTemp)

p4<-Temp_range%>%filter(!Season == "Spring")%>%
  ggplot(aes(x=as.POSIXct(Date), y=deltaTemp, colour=Site))+
  geom_line()+
  scale_color_manual(values=sitecolor)+xlab("Date")+facet_wrap(vars(Season),scales = "free_x")+
  scale_x_datetime(date_breaks = "3 week",date_labels = "%b %d")+
  theme_bw()+xlab("Date")+ylab("Daily temperature\n range (C\u00b0)")+
  theme(legend.position = "top",text = element_text(size=10))

##-----------------------------
##Create Figure with all panels
##-----------------------------
    
  ggarrange(p1,wsplot,p3,p4, ncol=1, common.legend = TRUE, labels = c("a)","b)","c)", "d)"),align = "v")%>%
  ggsave(filename =paste(mywd,myplots, "Fig2_v2.png",sep = "/"), dpi=400,units = c("in"),width = 5,height = 8)
  

##----------------------------------------------
##Create light profile for WINTER
##----------------------------------------------

##Filter a single day
tanks_winter=MS9.preds%>%filter(Date=="2022-07-11")%>%select(Time,DateTime)%>%droplevels()
#create treatment
tanks_winter$treatment='Low light'
tanks_winter$cPAR=0
##Data measured at noon in Heat-stress Exp
tanks_winter$cPAR[tanks_winter$Time >= 91500 & tanks_winter$Time <= 123000]=128
tanks_winter$Season="Winter"
#Duplicate for highlight
hl=tanks_winter
hl$treatment='High light'
hl$cPAR[hl$Time >= 91500 & hl$Time <= 123000]=540


tanks_winter%<>%rbind(hl)
treat=unique(tanks_winter$treatment)
predtw=tanks_winter%>%filter(Time %in% c(60002, 91500:123000, 170002))%>%droplevels()%>%arrange(DateTime)

##----------------------------------------------
##---- Estimate DLI exposure during winter experiments
##----------------------------------------------

winter_preds=list()

for (i in 1:length(treat)){
  d_data=predtw%>%filter(treatment==treat[i])%>%unique()%>%dplyr::arrange(Time)
  
  model=gam(cPAR~s(Time),data=d_data)
  
  timesequence=seq(from=d_data$DateTime[1], to=d_data$DateTime[length(d_data$DateTime)], by="sec")
  timeseq=format(as.POSIXct(timesequence), format = "%H%M%S") 
  D_newd=data.frame("Time"=as.numeric(timeseq),"DateTime"=timesequence)
  
  pred_PAR<-predict.gam(model,D_newd)
  
  D_newd$PAR=pred_PAR
  D_newd$DLI=D_newd$PAR/1000000#convert predicted PAR umol.m2.s to mol units (DLI)
  D_newd$Treatment=treat[i]
  
  winter_preds=c(winter_preds,list(D_newd))
}


winter_preds=do.call("rbind",winter_preds)
winter_preds%<>%filter(PAR>0)

#SUM instantaneous PAR (PPFD) over 12h period 
WLL=sum(winter_preds$DLI[winter_preds$Treatment=='Low light'])## 4 DLI 
WHL=sum(winter_preds$DLI[winter_preds$Treatment=='High light'])## 17 DLI

##Add time for sampling
sampling<-data.frame(Time=rep(233002,2),##Time of lights off, back to ambient temp and measure of Fv/Fm
                     DateTime=rep("2022-07-11 23:30:02",2),
                     treatment=c('Low light','High light'),       
                     cPAR=rep(0,2), 
                     Season="Winter")

predtw%<>%rbind(sampling)%>%arrange(DateTime)

##Repeat sequence for plotting
dayafter<-predtw%>%filter(DateTime < "2022-07-11 12:20:02")
dayafter%<>%mutate(DateTime= gsub("-11","-12",DateTime))

predtw%<>%rbind(dayafter)


##Add Temp 
tempdat=tanks_winter%>%filter(Time %in% c(60002, 91500:123000,153002,180002))%>%droplevels()%>%arrange(DateTime)
tempdat$Temp[tempdat$Time >=60002 & tempdat$Time < 93000]=22
tempdat$Temp[tempdat$Time >=120000 & tempdat$Time < 160000]=31

ramping=length(tempdat$Temp[tempdat$DateTime <"2022-07-12 07:00:00" & tempdat$Time >93000 & tempdat$Time <120000])
tempdat$Temp[tempdat$DateTime < "2022-07-12 07:00:00" & tempdat$Time >093000 & tempdat$Time < 120000]=seq(from=22, to=31, length.out=ramping)
tempdat$Heat="Heated"
tempdat$Temp[is.na(tempdat$Temp)==T]=22
#tempdat$Heat[is.na(tempdat$Temp)==T]="Control"

dayafter$Temp=22
dayafter$Heat="Control"

tempdat%<>%rbind(dayafter) 


##Create combined graph with two axes 
ylim.prim <- range(tempdat$cPAR) 
ylim.sec <- range(tempdat$Temp)#[!is.na(test$temp)]

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1] 

graphLabels2 <- data.frame(Season = rep(c("Summer", "Winter"),2),
                          treatment=c(rep("Low light",2), rep("High light",2)),
                          Date=rep(c(as.POSIXct("2022-12-14 09:30:02"),as.POSIXct("2022-07-12 09:30:02")),2),
                          cPAR= rep(c(30,25),2),
                          ref2 =c("Color","Color","",""))
graphLabels3=graphLabels2
graphLabels3$Date=rep(c(as.POSIXct("2022-12-13 20:30:02"),as.POSIXct("2022-07-11 20:30:02")),2)
graphLabels3$ref2=c("FvFm","FvFm","","")
  
labels=c("Heated", "High light","Low light","Control")
#mycolors=c("#ED1F1F","#FAA11B","#D9C582","blue")
mycolors=c("blue","#ED1F1F","#FAA11B","#D9C582")

ggplot(predtw, aes(x=as.POSIXct(DateTime), y=cPAR, group=treatment))+
  geom_line(size=1, linetype="dashed",aes(group=treatment,color=treatment))+
  geom_line(data=tempdat,aes(y=a+Temp*b, group=Heat, color=Heat),size=1)+#, linetype="dotted")+
  geom_hline(aes(yintercept=0, linetype="Control"),color="blue",size=0.70, linetype="solid")+
  scale_y_continuous("PAR (µmol m\u207B\u00b2 s\u207B\u00b9)", sec.axis = sec_axis(~ (. - a)/b, name = "Temperature (\u00b0C)")) +
  scale_color_manual(values=mycolors)+
  scale_x_datetime(date_breaks = "2 hours",date_labels = "%H")+
  xlab("time")+
  geom_point(data=tanks_winter%>%filter(cPAR==0),aes(x=as.POSIXct("2022-07-11 20:00:02",y =0)), color= "black", size=3)+
  geom_point(data=tanks_winter%>%filter(cPAR==0),aes(x=as.POSIXct("2022-07-12 09:00:02",y =0)), color= "black", size=3)+
  geom_text(data=subset(graphLabels3, Season=="Winter"),aes(y=cPAR, x=as.POSIXct(Date),label=ref2),color="black")+
  geom_text(data=subset(graphLabels2, Season=="Winter"),aes(y=cPAR, x=as.POSIXct(Date),label=ref2),color="black")+
  xlab("Time (Day hours)")+ggtitle("")+labs(color="Treatment")+
  theme_bw()+theme(legend.position = "bottom", text = element_text(size=13), legend.text = element_text(size=12))+
  scale_linetype_manual(name="", values=c("dashed"))+guides(color = guide_legend(nrow = 2))

ggsave(filename =paste(myplots, "Winter_tankprofile.png",sep = "/"), dpi=400,units = c("in"),width = 5,height = 4)


##----------------------------------------------
##Create light profile for SUMMER and estimate DLI
##----------------------------------------------
##Filter a single day
tanks_summer=MS9.preds%>%filter(Date=="2022-12-13")%>%select(Time,DateTime)%>%droplevels()
#create treatment
tanks_summer$treatment='Low light'
tanks_summer$cPAR=0
##light levels measured at noon and 3pm before ramping down 
tanks_summer$cPAR[tanks_summer$Time >= 83000 & tanks_summer$Time <= 153000]=238
tanks_summer$Season="summer"

#Duplicate for highlight
hls=tanks_summer
hls$treatment='High light'
hls$cPAR[hls$Time >= 83000 & hls$Time <= 153000]=508

tanks_summer%<>%rbind(hls)%>%arrange(Time)
treat=unique(tanks_summer$treatment)
##Filter time to predict on data
predts=tanks_summer%>%filter(Time %in% c(53002, 60002, 83000:153000))%>%droplevels()%>%arrange(DateTime)

#Adjust with obs data
##The average PAR measured at 6am: High: 133 (80-150) and Low:  66 (77-27)
predts$cPAR[predts$treatment=='Low light' & predts$Time==60002]=66
predts$cPAR[predts$treatment=='High light' & predts$Time==60002]=133
predts$cPAR[predts$Time==53002]=0

nd<-data.frame(Time=rep(190002,2),#extend daylight
               DateTime=rep(as.POSIXct("2022-12-13 19:00:02"),2),
        treatment= c('Low light','High light'),
        cPAR=rep(0,2),
        Season=rep("summer",2))       
predts%<>%rbind(nd)%>%arrange(Time)

##----------------------------------------------
##---- Estimate DLI exposure during Summer experiments
##----------------------------------------------
summer_preds=list()

for (i in 1:length(treat)){
  d_data=predts%>%filter(treatment==treat[i])%>%unique()%>%dplyr::arrange(Time)
  
  model=gam(cPAR~s(Time),data=d_data)
  
  timesequence=seq(from=d_data$DateTime[1], to=d_data$DateTime[length(d_data$DateTime)], by="sec")
  timeseq=format(as.POSIXct(timesequence), format = "%H%M%S") 
  D_newd=data.frame("Time"=as.numeric(timeseq),"DateTime"=timesequence)
  
  pred_PAR<-predict.gam(model,D_newd)
  
  D_newd$PAR=pred_PAR
  D_newd$DLI=D_newd$PAR/1000000####convert to mol (DLI)
  D_newd$Treatment=treat[i]
  
  summer_preds=c(summer_preds,list(D_newd))
}


summer_preds=do.call("rbind",summer_preds)

#SUM instantaneous PAR over 12h period  (i.e. predicted PAR umol.m2.s)/ 1million to convert to mol
SLL=sum(summer_preds$DLI[summer_preds$Treatment=='Low light'])## 11
SHL=sum(summer_preds$DLI[summer_preds$Treatment=='High light'])## 24


##-------------------------------------
##Plot experimental DLI vs observed
##-------------------------------------

meanpreds=MS9.preds%>%group_by(Site,Season,factor,Date)%>%summarise(DLI=sum(DLI))

#Aggregate empirical data
  expdata<-data.frame(Season=rep(c("Winter","Summer"),2),
                      Site=rep("Bundegi",4),
                      DLI=c(WLL,SLL,WHL,SHL),
                      ref=c("","Low","","High"))
  
  quantile(meanpreds$DLI[meanpreds$Site=="Bundegi" & meanpreds$Season=="Summer"])
  # 0%      25%      50%      75%     100% 
  # 22.12662 27.53021 29.13176 30.62358 34.27847 
  quantile(meanpreds$DLI[meanpreds$Site=="Bundegi" & meanpreds$Season=="Winter"])
  # 0%       25%       50%       75%      100% 
  # 6.851163 13.441380 14.892591 15.586080 19.387702 

  
    meanpreds%>%filter(!Season == "Spring")%>%
    ggplot(aes(x=Site, y=DLI,fill=Site, shape=Site))+
      geom_violin(alpha=0.6)+
      scale_color_manual(values = sitecolor, aesthetics = "fill")+
    geom_boxplot(width = .1, fill = "black", outlier.colour = NA) +
    stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21, size = 2.5)+
      facet_wrap(vars(Season))+
      geom_point(data=expdata,aes(x=Site,y=DLI),colour=rep(c("gold","darkorange"),each=2), size=4,
      position = position_nudge(0.5))+
      geom_text(data=subset(expdata,ref=="Low"),aes(label=ref,y=SLL-2),colour="gold", size=5, position = position_nudge(0.5))+
      geom_text(data=subset(expdata,ref=="High"),aes(label=ref,y=SHL-4),colour="darkorange", size=5, position = position_nudge(0.5))+
      ylab("Daily Light Integral\n DLI (mol photons m\u207B\u00b2)")+
    theme_bw()+theme(legend.position = "bottom", text = element_text(size=13), legend.text = element_text(size=12))

    ggsave(filename =paste(mywd,myplots, "obs_vs_exp_DLI.png",sep = "/"), dpi=400,units = c("in"),width = 5,height = 4)


    
    ##--------------------------------------------------------------------------------------
   ##  CHeck DLI on days where daily SST exceeded MMM+1 compared to days without thermal anomalies
    ##--------------------------------------------------------------------------------------
    
    Temp_range$anomaly=0
    Temp_range$anomaly[Temp_range$maxTemp> MMM+1]=1 ##Summer only
    
    dhw.DLI=Temp_range%>%mutate(Date=as.Date(Date))%>%left_join(meanpreds)%>%filter(!is.na(DLI))%>%filter(Season=="Summer")

    #tapply(dhw.DLI$DLI, dhw.DLI$Season, mean)
    tapply(dhw.DLI$DLI, dhw.DLI$anomaly, quantile)
    
