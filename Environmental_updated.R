##Final Plots
rm(list=ls())

mywd<-"~/Documents/2024"
myplots=paste(mywd, "Logger Plots", sep="/")


library(ggplot2);library(tidyr);library(reshape2);library(dplyr); library(mgcv);library(lubridate);library(magrittr);
library(knitr);library(scales); library(readxl);library(ggpubr);library(ggridges);library(readr);library(stringr);
library(ggh4x) ##To color code the facets
library(pracma)##Find outliers in time-series
library(forecast)

options(mc.cores = parallel::detectCores())##optimise model run


####Load Script outputs*******************************************************
load(paste(mywd,"Data/Environmental.tables.Rdata",sep="/")) 


sitecolor=scales::viridis_pal(end=0.7)(length(unique(DLI.preds$Site)))
ribbon_colors <- c("Bundegi" = sitecolor[1], "Tantabiddi" = sitecolor[2])
se <- function(x) sd(x) / sqrt(length(x))


##-----------------------------
## Figure S1
##-----------------------------
mydata$MS9depth%>%
  group_by(Site,Season)%>%arrange(DateTime)%>%#Complete missing dates to avoid connecting lines
  complete(Date=seq.Date(min(Date), max(Date), by="day"))%>%
  ggplot(aes(x=as.POSIXct(Date), y=mDEPTH, color=Site))+
  #Identify Seasons for reference
  annotate("rect", xmin=as.POSIXct("2022-07-01"), xmax=as.POSIXct("2022-09-30"), 
           ymin=min(mydata$MS9depth$mDEPTH), ymax=max(mydata$MS9depth$mDEPTH), alpha=.15, fill="grey")+
  annotate("rect", xmin=as.POSIXct("2023-06-01"), xmax=as.POSIXct("2023-09-30"), 
           ymin=min(mydata$MS9depth$mDEPTH), ymax=max(mydata$MS9depth$mDEPTH), alpha=.15, fill="grey")+
  annotate("rect", xmin=as.POSIXct("2022-12-01"), xmax=as.POSIXct("2023-03-31"), 
           ymin=min(mydata$MS9depth$mDEPTH), ymax=max(mydata$MS9depth$mDEPTH), alpha=.1, fill="orange")+
  annotate("rect", xmin=as.POSIXct("2023-12-01"), xmax=as.POSIXct("2024-01-10"), 
           ymin=min(mydata$MS9depth$mDEPTH), ymax=max(mydata$MS9depth$mDEPTH), alpha=.1, fill="orange")+
  scale_color_viridis_d(name = "Site", alpha = 0.3, end=0.7)+
  #Add seasons
  annotate("text", x =as.POSIXct("2022-07-30"), y = 1.3, label = "Winter", size=3)+
  annotate("text", x =as.POSIXct("2023-07-30"), y = 1.3, label = "Winter", size=3)+
  annotate("text", x =as.POSIXct("2023-02-01"), y = 1.3, label = "Summer", size=3)+
  annotate("text", x =as.POSIXct("2023-12-31"), y = 1.3, label = "Summer", size=3)+
  annotate("text", x =as.POSIXct("2022-07-10"), y = 4.8, label = "Bundegi", size=3, color="#440154FF")+
  annotate("text", x =as.POSIXct("2022-10-01"), y = 4.8, label = "Tantabiddi", size=3, color="#43BF71FF")+
  geom_line(linewidth=0.1)+labs(x="Date", y=ylab)+
  scale_x_datetime(date_breaks = "10 week",date_labels = "%b-%d-%y")+
  ylab("Logger depth (m)")+
  theme_bw()+theme(legend.position = "none",text = element_text(size=6))

ggsave(file=paste(myplots,"Depth.jpg",sep = "/"), dpi=300, units = c("in"),width = 5,height = 2.5)



##------------
##FIGURE 1
##------------

# Mean daily PAR per Site

f2a<-MS9.mean%>%group_by(Site,Date,Month,Season)%>%
  summarise(across(calcPAR_mean,list(mean=mean, se=se)))%>%
  arrange(Date)%>%group_by(Site)%>%
  complete(Date=seq.Date(min(Date), max(Date), by="day"))%>%##Complete missing days for plotting
  ggplot(aes(x = as.Date(Date),colour=Site)) +
  geom_line(aes(y = calcPAR_mean_mean), alpha = 0.5, size = 1) +
  geom_ribbon(aes(ymin = calcPAR_mean_mean-calcPAR_mean_se, 
                  ymax = calcPAR_mean_mean+calcPAR_mean_se, fill=Site),alpha = 0.2, color=NA) +
  scale_y_continuous(expand = c(0, 0))+
  ylab(expression(atop("Daily averaged PAR",(µmol~m^-2~s^-1))))+
  xlab("")+
  scale_color_manual(values = sitecolor) +
  scale_fill_manual(values = ribbon_colors) +  
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 8),
        panel.spacing.x = unit(5, "mm"))+
  annotate(geom="text", x=summary_PAR$YrM[3], y=80, label="2022", size=2)+
  annotate(geom="text", x=summary_PAR$YrM[14], y=80, label="2023", size=2)+
  annotate("text", x =as.Date("2022-08-10"), y = 720, label = "Bundegi", size=3, color="#440154FF")+
  annotate("text", x =as.Date("2022-08-10"), y = 650, label = "Tantabiddi", size=3, color="#43BF71FF")


##-----------------------------
##ADD DLI scale to PAR plot
##-----------------------------

DLI.preds$Date<-as.Date(DLI.preds$Date, format="%d-%m-%Y")
site.DLI.preds<-DLI.preds%>%
  group_by(Site,Date,Month,Season)%>%
  summarise(across(DLI,list(mean=mean)))%>%
  arrange(Date)%>%group_by(Site)
site.DLI.preds$Site<-as.factor(site.DLI.preds$Site)
levels(site.DLI.preds$Site)<-c("Bundegi","Tantabiddi")

site.MS9.mean<-MS9.mean%>%group_by(Site,Date,Month,Season)%>%
  summarise(across(calcPAR_mean,list(mean=mean, se=se)))%>%
  arrange(Date)%>%group_by(Site)

site.MS9.mean%<>%left_join(site.DLI.preds)
# Calculate the range of temperature
dli_range <- range(site.MS9.mean$DLI_mean)
# Calculate the range of light
par_range <- range(site.MS9.mean$calcPAR_mean_mean)

# Calculate the scaling factor for dual axes
scaling_factor <- diff(dli_range) / diff(par_range)

f2a=f2a+scale_y_continuous(sec.axis = sec_axis(~.*scaling_factor, 
                            name=expression(atop("Daily Light Integrals",(mol~photons~m^-2~d^-1)))))



MMM=27
SM=22
##Mean Temp
f2b<-Temp%>%
  group_by(Site,as.POSIXct(Date))%>%
  ggplot(aes(x=as.POSIXct(Date),y=Temp, color=Site,fill=Site))+geom_line(alpha = .4)+
  stat_summary_bin(fun = "mean",
                   geom = "ribbon",
                   alpha = .3,
                   color=NA,
                   fun.max = min,
                   fun.min = max,breaks = seq(0, 100, 10))+
  geom_hline(aes(yintercept = MMM), color="black",linetype="dashed")+
  geom_hline(aes(yintercept = SM), color="black",linetype="dotted",alpha=0.7)+
  ylab(expression(atop("Daily temperature", "at depth (\u00b0C)"))) +
  xlab("") +
  scale_color_manual(values = sitecolor) +
  scale_fill_manual(values = ribbon_colors) + 
  scale_x_datetime(date_breaks = "1 months", date_labels = "%b",expand = c(0, 0))+
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 8),
        panel.spacing.x = unit(5, "mm"))+
  annotate(geom="text", x=as.POSIXct("2022-09-01"), y=18.2, label="2022", size=2)+
  annotate(geom="text", x=as.POSIXct("2023-08-01"), y=18.2, label="2023", size=2)


##Temp range
f2c<-Temp_range%>%
  group_by(Site)%>%complete(Date=seq.Date(min(Date), max(Date), by="day"))%>%
  ggplot(aes(x=as.POSIXct(Date), y=deltaTemp, colour=Site))+
  geom_line()+
  scale_color_manual(values=sitecolor)+xlab("Date")+
  scale_x_datetime(date_breaks = "1 months", date_labels = "%b",expand = c(0, 0))+
  xlab("")+
  ylab(expression(atop("Daily temperature", "range (\u00b0C)")))+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 8),
        panel.spacing.x = unit(5, "mm"))+
  annotate(geom="text", x=as.POSIXct("2022-09-01"), y=0.02, label="2022", size=2)+
  annotate(geom="text", x=as.POSIXct("2023-08-01"), y=0.02, label="2023", size=2)



Temp_range%<>%left_join(Temp%>%select(-Temp))%>%unique()


##----------------------------------------------
##Plot combined monthly predicted time-series
##----------------------------------------------

variables<-c("Temp_mean","mean_PAR")
variablecolor=scales::viridis_pal(option="C",begin=0.9,end=0.6)(length(unique(summary_PAR$Site)))
variablecolor[1]="#FEBA2CFF"

# Calculate the range of temperature
temp_range <- range(summary_PAR$Temp_mean)
# Calculate the range of light
light_range <- range(summary_PAR$mean_PAR)

# Calculate the scaling factor for dual axes
scaling_factor <- diff(temp_range) / diff(light_range)
scaling_factor<-0.048## adjusted for plotting

strip <- strip_themed(background_x = elem_list_rect(fill = alpha(sitecolor, 0.3)))

tempbreaks=c(14,17,20,23,26,29,32)
lightbreaks=c(380,480,580)

f2d<-summary_PAR%>%
  ggplot(aes(x = YrM)) +
  geom_ribbon(aes(ymin = Temp_mean - Temp_se, ymax = Temp_mean + Temp_se, fill = "Temperature"), alpha = 0.2) +
  geom_ribbon(aes(ymin = mean_PAR * scaling_factor - se_PAR * scaling_factor, 
                  ymax = mean_PAR * scaling_factor + se_PAR * scaling_factor, fill = "Light"), alpha = 0.2) +
  geom_line(aes(y = Temp_mean, color = "Temperature"), alpha = 0.5, size = 2) +
  geom_line(aes(y = mean_PAR * scaling_factor, color = "Light"), alpha = 0.5, size = 2) +
  scale_y_continuous(sec.axis = sec_axis(~./scaling_factor,
                                         name = expression(atop("Monthly averaged PAR",(µmol~m^-2~s^-1))),
                           breaks=lightbreaks), breaks=tempbreaks)+
  ylab(expression(atop("Monthly averaged"," temperature (C\u00b0)")))+
  xlab("Date") +
  scale_color_manual(values = variablecolor) +
  scale_fill_manual(values = variablecolor) +  # Add this line to adjust fill colors
  scale_x_date(date_breaks = "2 month", date_labels = "%b", expand = c(0, 0)) +
  facet_grid2(~Site, strip=strip)+
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 8),
        panel.spacing.x = unit(3, "mm"),##Avoid overlapping axes
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color = variablecolor[2]), # Adjust colors here
        axis.title.y = element_text(color = variablecolor[2]), # Adjust colors here
        axis.text.y.right = element_text(color = variablecolor[1]),
        axis.title.y.right = element_text(color = variablecolor[1]))+
  annotate(geom="text", x=summary_PAR$YrM[3], y=15, label="2022", size=2)+
  annotate(geom="text", x=summary_PAR$YrM[14], y=15, label="2023", size=2)


##Create Figure 2 with all panels

ggarrange(f2a,f2b,f2c,f2d, ncol=1, labels = c("a","b","c","d"))%>%
  ggsave(filename =paste(myplots, "Fig2.png",sep = "/"), dpi=400,units = c("in"),width = 6,height = 8)


##----------
## FIGURE 3
##----------

# --------------------------
# #Aggregate empirical data
expdata<-data.frame(Season=rep(c("Winter","Summer"),2),
                    Site=rep("BU",4),
                    DLI=c(4,11,16,24),
                    ref=c("","Low","","High"))
expdata$Season<-factor(expdata$Season, levels=c("Winter","Summer"))

##Add Site Label
ann_text <- data.frame(Site = "BU",DLI = 31,lab = "Text",
                       Season = factor("Winter",levels = c("Winter","Summer")))


##Order for plotting
DLI.preds$Season<-factor(DLI.preds$Season, levels=c("Autumn","Winter","Spring","Summer"))


dli1<-DLI.preds%>%filter(Season %in% c("Summer","Winter"))%>%
  ggplot(aes(x=Site, y=DLI,fill=Site, shape=Site))+
  geom_text(data=ann_text, label="Bundegi",colour="#440154FF")+
  geom_text(data=ann_text, label="Tantabiddi",colour="#43BF71FF",nudge_x = 1)+
  geom_violin(alpha=0.6)+
  scale_color_manual(values = sitecolor, aesthetics = "fill")+
  geom_boxplot(width = .1, fill = "black", outlier.colour = NA) +
  stat_summary(fun.y = mean, geom = "point", fill = "white", shape = 21, size = 2.5)+
  facet_wrap(vars(Season))+xlab("")+
  ylab(expression(atop("Daily Light Integrals","DLI (mol photons m\u207B\u00b2)")))+
  geom_point(data=expdata,aes(x=Site,y=DLI),colour=rep(c("gold","darkorange"),each=2), size=3,
             position = position_nudge(0.5))+
  geom_text(data=subset(expdata,ref=="Low"),aes(label=ref,y=DLI-2),colour="gold", size=4, position = position_nudge(0.5))+
  geom_text(data=subset(expdata,ref=="High"),aes(label=ref,y=DLI-4),colour="darkorange", size=4, position = position_nudge(0.5))+
  theme_bw()+theme(legend.position = "none", text = element_text(size=8),
                   axis.text.x = element_blank(), axis.ticks.x = element_line(colour = "white"))


# #Aggregate empirical data
expdata2<-
  data.frame(Month=c(rep("December",2),rep("March",2),rep("July",2)),
             Site=rep(c("BU","TN"),3),
             DLI=c(11,24,rep(NA,2),4,16),
             ref=rep(c("Low","High"),3),
             Season=c(rep("Summer",4),rep("Winter",2)))


dli2<-DLI.preds%>%filter(Month%in%c(month.name[c(3,7,12)]))%>%
  ggplot(aes(x=Month, y=DLI, fill=Site, Shape=Site))+
  geom_violin(alpha=0.6,position = position_dodge(0.5),width=1)+
  scale_color_manual(values = sitecolor, aesthetics = "fill")+
  geom_boxplot(width = .1, fill = "black", outlier.colour = NA,position = position_dodge(0.5)) +
  stat_summary(fun.y = mean, geom = "point", fill = "white", shape = 21, size = 2.5,position = position_dodge(0.5))+
  facet_wrap(~factor(Season, levels=c("Winter","Summer")),scales="free_x")+
  ylab(expression(atop("Daily Light Integrals","DLI (mol photons m\u207B\u00b2)")))+
  xlab("")+
  geom_point(data=expdata2,aes(x=Month,y=DLI),colour=rep(c("gold","darkorange"),3), size=3,
             position = position_nudge(0.1))+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=8),
        axis.text.x = element_text(margin = margin(-20, 0, 0, 0), size=8),##MOve axis to inside of the plot
        axis.ticks.x = element_line(colour = NA),plot.margin = margin(10, 10, 10, 10, "pt"), 
        panel.background = element_rect(fill = "transparent", color = NA))



ggarrange(dli1,dli2, ncol=1, labels = c("a","b"), align = "hv")%>%
  ggsave(filename =paste(myplots, "Fig3.png",sep = "/"), dpi=500,units = c("in"),width = 4,height = 5)



##-------------------------------------------------------------------------
## Summary stats
MS9.mean%<>%na.omit()
##-------------------------------------------------------------------------

