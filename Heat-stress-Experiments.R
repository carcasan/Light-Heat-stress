rm(list=ls())

library(tidyverse)
library(ggplot2)
library(magrittr)

library(colorspace)
library(ggridges)
library(lubridate)
library(ggpubr)
library(rstatix)

library(RColorBrewer)
library(viridis)
library(wesanderson)
library(brms)
library(cmdstanr)
set_cmdstan_path("C:/cmdstan")
options(brms.backend = "cmdstanr")
check_cmdstan_toolchain()
## If error 
##dir.create(tempdir())
##check_cmdstan_toolchain()

library(ggmcmc)##to convert outputs to tibble for plots
library(ggthemes)
library(gridExtra)
library(stringr)
library(patchwork)

options(mc.cores = parallel::detectCores())##optimise model run


##Load Data
mywd="C:/Users/ccastros/OneDrive - Australian Institute of Marine Science/Documents/2023/Light/RCode-Github"
load(paste(mywd,"/Experiments_Tables.R", sep=("")))


##***************************************************
##***************************************************
## Bayesian Generalised Random-effects Model -- brms
##***************************************************
##***************************************************


model_winter_3f<-brm(FvFm~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), data=Winterdata, family = Beta(),warmup = 1000, iter=3000,adapt_delta=0.9)
ppcheck_wm1<-pp_check(model_winter_3f, ndraws = NULL)##estimated predictions vs observed data -OK 
##describe priors
brms::prior_summary(model_winter_3f)

#predicted probabilities- main effects
mycolors=rep(c("#EBF0AF", "#E0B336","#F08622", "#D11450"),2)## in order of treatments


##New data
nd <- Winterdata |>
  dplyr::distinct(Site,Temp, Light, Treatment)

##Predict
mwnd2 <- brms::posterior_epred(##already in response scale
  model_winter_3f, re_formula = .~ 1,##same as re_formula = NA,
  newdata = nd
) |>
  t() %>%
  cbind(nd, .)# 

mwinnd2<-mwnd2%>%dplyr::select(-c(Temp,Light))%>%tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment)%>%left_join(unique(mwnd2[c(2:4)]))%>%#droplevels()%>%
  ggplot(aes(y = Treatment, x = effect,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  geom_point(data=Winterdata%>%group_by(Tank,Genotype,Treatment,Light,Temp)%>%
  #geom_point(data=Winterdata%>%group_by(Site,Tank,Genotype,Treatment,Light,Temp)%>%
               summarise(FvFm=mean(FvFm)), aes(y=Light, x=FvFm, fill=Treatment, colour=Treatment), alpha=0.3, size=0.7)+
  ggdist::stat_halfeye(aes(y = Light, x = effect, group=Temp),color="grey",.width = c(0.8, 0.95), 
                       point_interval = "median_hdi",alpha=0.6,size=0.1)+
  facet_grid(vars(Light),scales = "free")+xlab("Posterior FvFm")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  scale_x_continuous(trans="logit", limits = c(0.1, 0.75))+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=7),text = element_text(size=8), axis.text = element_text(size=8), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

##-------------------
##Estimate contrasts
##-------------------
wcontrasts=model_winter_3f %>%
  emmeans::emmeans(~ Treatment)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>% as.data.frame(wcontrasts)
wcontrasts$model="winter.FvFm"

wcontrasts_sites=model_winter_3f %>%
  emmeans::emmeans(~ Site)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame(wcontrasts_sites)
wcontrasts_sites$model="winter.FvFm"

Contrasts=rbind(wcontrasts,wcontrasts_sites)


##---------------------------------------------------------------------------------
##compute posterior distribution of genotype-level response -predicted for new data
##---------------------------------------------------------------------------------

##Create New data  and Predict
gen<- Winterdata |>
  dplyr::distinct(Temp, Light, Site, Genotype,Treatment) 


pred_m3b <- brms::posterior_epred(
  model_winter_3f, re_formula = .~ 1 + (Temp * Light || gr(Genotype, by = Site)),
  newdata = gen,
) |>
  t() %>%
  cbind(gen, .)


##--------------------------------------------------------------------
##################-----SUMMER---------################################
##--------------------------------------------------------------------

##How many zeros
Summerdata %>% 
  count(FvFm == 0) %>% 
  mutate(prop = n / sum(n))## 22 (or 3% of the data)

##***************************************************
## Bayesian Generalised Random-effects Model -- brms
##***************************************************

model_summer_1b<-brm(bf(FvFm~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), zi~1), data=Summerdata, family = zero_inflated_beta(link = "logit"), warmup = 1000, iter=3000)##is this best?
ppcheck_sm1<-pp_check(model_summer_1b, type = "dens_overlay", ndraws = NULL)
##describe priors
brms::prior_summary(model_summer_1b)


##Predict
msnd <- brms::posterior_epred(
  model_summer_1b, re_formula = .~ 1,
  newdata = nd#,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)

##Create New data for Genotypes
gen2<- Summerdata |>
  dplyr::distinct(Temp, Light, Site, Genotype,Treatment) 


msfig2<-msnd%>%dplyr::select(-c(Temp,Light))%>%tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment)%>%left_join(msnd[c(2:4)])%>%
  ggplot(aes(y = Treatment, x = effect,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  geom_point(data=Summerdata%>%group_by(Tank,Genotype,Treatment,Light,Temp)%>%
  #geom_point(data=Summerdata%>%group_by(Site,Tank,Genotype,Treatment,Light,Temp)%>%
               summarise(FvFm=mean(FvFm)), aes(y=Light, x=FvFm, fill=Treatment, colour=Treatment), 
             alpha=0.3, size=0.7)+
  ggdist::stat_halfeye(aes(y = Light, x = effect, group=Temp),color="grey",.width = c(0.95), 
                       point_interval = "median_hdi",alpha=0.6,size=0.1)+
  facet_grid(vars(Light), scales = "free")+xlab("Posterior FvFm")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  scale_x_continuous(trans="logit", limits = c(0.1, 0.75))+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=7),text = element_text(size=8), axis.text = element_text(size=8), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

ggarrange(mwinnd2,msfig2+ylab(""), common.legend = T, legend = "bottom",labels = "auto")%>%
 ggsave(filename =paste(myplots, "Experiments/FvFm_Seasons_2.png",sep = "/"), dpi=400, width = 6, height = 3)

##-------------------
##Estimate contrasts
##-------------------
Scontrasts=model_summer_1b %>%
  emmeans::emmeans(~ Treatment)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame(Scontrasts)

Scontrasts$model="Summer.FvFm"


Scontrasts_site=model_summer_1b %>%
  emmeans::emmeans(~ Site)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame(Scontrasts_site)

Scontrasts_site$model="Summer.FvFm"
Contrasts%<>%rbind(Scontrasts,Scontrasts_site)


##compute posterior draws of genotype-level effects predicted for new data
pred_smod <- brms::posterior_epred(
  model_summer_1b, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen2#, transform = TRUE
) |>
  t() %>%
  cbind(gen2, .)


###########################
##RedChannel Response 
###########################

##PLOT  RAW data 

Summerdata$Season="Summer"
Winterdata$Season="Winter"

alldata= rbind(Summerdata,Winterdata)

alldata$SeasonTemp=paste(alldata$Season,alldata$Temp, sep=".")

ggplot(alldata, aes(y=Light, x=RedChannel*255, fill=SeasonTemp))+
  geom_density_ridges(rel_min_height = 0.005,scale = 1.2)+
  scale_fill_discrete_sequential(palette="Heat2", aesthetics = c("color","fill"), alpha=0.9)+#+facet_grid(~Site)
  xlab("Red Channel intensity")+ylab("Light")+#facet_wrap(~Site)+
  theme(legend.position = "bottom",text = element_text(size=8))+
  scale_x_continuous(breaks = c(0,50,100,150,200,250),label = c(0,50,100,150,200,250))

#ggsave(file=paste(myplots,"Experiments/Seasonal_Redchannel.jpeg",sep="/"),width = 5, height = 3.5, units = c("in"),dpi = 300)

##***************************************************
## Bayesian Generalised Random-effects Model -- brms
##***************************************************

model_winter_color_1<-brm(RedChannel~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), data=Winterdata, family = Beta(),warmup = 1000, iter=3000)
mwc_check<-pp_check(model_winter_color_1, ndraws=NULL)
##describe priors
brms::prior_summary(model_winter_color_1)

## Predict 
mwc <- brms::posterior_linpred(
  model_winter_color_1, re_formula = .~ 1,
  newdata = nd,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)# 


mwcfig2<-mwc%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment, group=Site)%>%left_join(mwc[c(2:4)])%>%
  ggplot(aes(y = Treatment, x = effect*255,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  geom_point(data=Winterdata%>%group_by(Tank,Genotype,Treatment,Light,Temp)%>%
  #geom_point(data=Winterdata%>%group_by(Site,Tank,Genotype,Treatment,Light,Temp)%>%
               summarise(RedChannel=mean(RedChannel)), aes(y=Light, x=RedChannel*255, fill=Treatment, colour=Treatment), alpha=0.3, size=0.9)+
  ggdist::stat_halfeye(aes(y = Light, x = effect*255),color="grey",.width = c(0.8, 0.95), point_interval = "median_hdi",alpha=0.6,size=2)+
  facet_grid(vars(Light),scales="free")+xlab("Posterior \nRed channel intensity")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  theme_bw()+
  theme(legend.position = "bottom",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=7),text = element_text(size=8), axis.text = element_text(size=8), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

##-------------------
##Estimate contrasts
##-------------------
mwc_contrasts=model_winter_color_1 %>%
  emmeans::emmeans(~ Treatment)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame()
mwc_contrasts$model="Winter.Color"

mwc_contrasts_sites=model_winter_color_1 %>%
  emmeans::emmeans(~ Site)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame()
mwc_contrasts_sites$model="Winter.Color"

Contrasts%<>%rbind(mwc_contrasts,mwc_contrasts_sites)

##compute posterior draws of genotype-level effects predicted for new data
pred_wrc_mod <- brms::posterior_linpred(
  model_winter_color_1, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen, transform = TRUE
) |>
  t() %>%
  cbind(gen, .)


########### SUMMER DATA###################

model_summer_color_1<-brm(RedChannel~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), data=Summerdata, family = Beta(),warmup = 1000, iter=3000)
model_summer_color_1
msccheck<-pp_check(model_summer_color_1, ndraws = NULL)##OK
##describe priors
brms::prior_summary(model_summer_color_1)


##Predict
csmod <- brms::posterior_linpred(
  model_summer_color_1, re_formula = .~ 1,
  newdata = nd,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)# 


mcsfig2<-csmod%>%dplyr::select(-c(Temp,Light))%>%tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment)%>%left_join(csmod[c(2:4)])%>%
  ggplot(aes(y = Light, x = effect*255,color=Treatment,fill=Treatment)) +
  geom_point(data=Summerdata%>%group_by(Tank,Genotype,Treatment,Light,Temp)%>%
               summarise(RedChannel=mean(RedChannel)), aes(y=Light, x=RedChannel*255, fill=Treatment, colour=Treatment), alpha=0.3, size=0.9)+
  ggdist::stat_halfeye(aes(y = Light, x = effect*255, group=Treatment),color="grey",.width = c(0.8, 0.95), point_interval = "median_hdi",alpha=0.6,size=2)+
  facet_grid(vars(Light), scales="free")+
  xlab("Posterior \nRed Channel intensity")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  theme_bw()+
theme(legend.position = "bottom",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
      legend.text =element_text(size=7),text = element_text(size=8), axis.text = element_text(size=8), 
      axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

##--PLOT
 ggarrange(mwcfig2,mcsfig2+ylab(""), common.legend = T, legend = "bottom",labels = "auto", align = "hv")%>%
   ggsave(filename =paste(myplots, "Experiments/Pop_Color2.png",sep = "/"), dpi=400, width = 6, height = 3)

##-------------------
##Estimate contrasts
##-------------------
SC_contrasts=model_summer_color_1 %>%
  emmeans::emmeans(~ Treatment)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame()
SC_contrasts$model="Summer.Color"

SC_contrasts_site=model_summer_color_1 %>%
  emmeans::emmeans(~ Site)%>%
  emmeans::regrid(transform="log")%>%## give percentage change response
  pairs%>%
  emmeans::regrid()%>%as.data.frame()
SC_contrasts_site$model="Summer.Color"

Contrasts%<>%rbind(SC_contrasts,SC_contrasts_site)

write.table(Contrasts, file=paste(mywd,"/Contrasts_table.csv", sep=""), sep=",",quote=FALSE, row.names=T, col.names=T)

##Plot Model fit for all
ggarrange(ppcheck_wm1,ppcheck_sm1,mwc_check,msccheck, labels = "auto")%>%
  ggsave(file=paste(myplots,"/Model/Modelfit.jpg", sep=""),width = 5.3, height = 6, units = c("in"),dpi = 300)

##------------------------------------------------------------------------------
##compute posterior draws of genotype-level effects predicted for new data
##------------------------------------------------------------------------------
pred_src_mod <- brms::posterior_linpred(
  model_summer_color_1, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen2,transform = TRUE
) |>
  t() %>%
  cbind(gen2, .)


##----------------------------------------------------------------------------
##Calculate magnitude of predicted shift and compared with change in FvFm response
##----------------------------------------------------------------------------

meanpreds<-pred_src_mod%>%tidyr::pivot_longer(cols=6:4005,names_to = "Iteration", values_to = "effect")%>%
group_by(Genotype,Treatment)%>%
tidybayes::mean_qi(Genotype_mean = effect*255)#, .width = c(.95, .66)) 

##Relative Change in Color
summercolor<-meanpreds%>%select(Genotype, Treatment, Genotype_mean)%>%spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

##Order Relative Change in Color
descolor<-summercolor%>%arrange(Change)%>%filter(Treatment=="High light")%>%select(color)
descolor=descolor$color
##Fix G2/G20
descolor[3]="#C75050"
descolor[4]="#3E45A8"

sc<-summercolor%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,Change), fill=Treatment))+
  xlab("Relative Change in Color\n (Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#FADD3C", "#F5970A"))+
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor
  glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "bottom", legend.title = element_blank(),text = element_text(size=8),legend.text = element_text(size=8))
  
##Predicted FvFm
meanpreds_fvfm=pred_smod%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Genotype,Treatment)%>%
  tidybayes::mean_qi(Genotype_mean = effect)

# Relative Change 
summerfvfm<-meanpreds_fvfm%>%select(Genotype, Treatment, Genotype_mean)%>%spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

##Order treatment levels
summerfvfm$Treatment<-factor(summerfvfm$Treatment,levels=c("Low light","High light"))

sp<-summerfvfm%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,summercolor$Change), fill=Treatment))+
  xlab("Relative Change in FvFm \n(Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#FADD3C", "#F5970A"))+
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  geom_vline(xintercept=-0.5,col = "gray60", linetype="dashed")+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor
    glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "bottom", legend.title = element_blank(),text = element_text(size=8),legend.text = element_text(size=8))+
  geom_text(x=-0.4, y="G11", label="Bundegi", color="#C75050", size=2)+
  geom_text(x=-0.4, y="G06", label="Tantabiddi", color="#3E45A8", size=2)

ggarrange(sp,sc+ylab("")+theme(axis.text.y = element_blank()), common.legend = T, legend = "bottom",labels = "auto", align = "hv")%>%
 ggsave(filename =paste(myplots, "Model/Figure5.png",sep = "/"), dpi=400, width = 4, height = 5)

#save(Summerdata, Winterdata, summercolor,summerfvfm,wintercolor,winfvfm,file="Experiments_Tables.R")
