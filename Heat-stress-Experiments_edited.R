rm(list=ls())

library(tidyverse)
library(ggplot2)
library(magrittr)
library(tidyr)

library(colorspace)
library(ggridges)
library(lubridate)
library(ggpubr)
library(rstatix)
library(ggh4x)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(ggsci)
library(brms)
library(cmdstanr)
library(ggeffects)
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
mywd="C:/Users/ccastros/OneDrive - Australian Institute of Marine Science/Documents/2023/Light/RCode-Github/Light-Heat-stress"
#load(paste(mywd,"/Experiments_Tables.R", sep=("")))
#predicted probabilities- main effects
load("bmrs model outputs.R")
myplots="C:/Users/ccastros/OneDrive - Australian Institute of Marine Science/Documents/2024/Logger Plots"

mycolors=rep(c("#EBF0AF", "#E0B336","#F08622", "#D11450"),2)## in order of treatments


##***************************************************
##***************************************************
## Bayesian Generalised Random-effects Model -- brms
##***************************************************
##***************************************************


model_winter_fvfm<-brm(FvFm~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), 
                  data=Winterdata, family = Beta(),warmup = 1000, iter=3000,adapt_delta=0.9)

ppcheck_wm1<-pp_check(model_winter_fvfm, ndraws = NULL)##estimated predictions vs observed data -OK 
##describe priors
#brms::prior_summary(model_winter_fvfm)


##New data for predicting responses (Global-no gens)
nd <- Winterdata |>
  dplyr::distinct(Site,Temp, Light, Treatment)##Excluding Genotypes



##Predict Marginal effects
predict.global.winter.fvfm <- brms::posterior_epred(##already in response scale
  model_winter_fvfm, re_formula = .~ 1,##same as re_formula = NA, (ignoring Genotype effects)
  newdata = nd
) |>
  t() %>%
  cbind(nd, .)# 

##AddSite Label
ann_text <- data.frame(Site = "Bundegi",Light = c("Low","High"), Treatment= c("Low.Control","High.Control"),
                       Season = c("Winter","Summer"), effect= 0.2)
ann_text$Light<-factor(ann_text$Light, levels=c("Low","High"))

##Plot posterior predictions
plot.winter.fvfm<-predict.global.winter.fvfm%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment)%>%left_join(unique(predict.global.winter.fvfm[c(2:4)]))%>%#droplevels()%>%
  ggplot(aes(y = Treatment, x = effect,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  ggdist::stat_halfeye(aes(y = Light, x = effect, group=Temp),color="grey",.width = c(0.8, 0.95), 
                       point_interval = "median_hdi",alpha=0.6,size=0.1)+
  facet_grid(vars(Light),scales = "free")+xlab("Posterior FvFm")+
  xlab("Posterior FvFm")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  scale_x_continuous(trans="logit")+
  geom_text(data=subset(ann_text,Light=="Low"),
            aes(label="Control",y=1.5, x=0.64),colour="gold", size=3, 
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="Low"),
            aes(label="Heated",y=1.2, x=0.64),colour="#F08622", size=3, 
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="High"),
            aes(label="Control",y=1.5, x=0.72),colour="#F08622", size=3,  
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="High"),
            aes(label="Heated",y=1.2, x=0.72),colour="#D11450", size=3, 
            position = position_nudge(0.0))+
  theme_bw()+
  theme(legend.position = "none",#legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        #legend.text =element_text(size=7),
        text = element_text(size=6), axis.text = element_text(size=6), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())


##--------------------------------------------------
##Estimate contrasts: mean marginal effect of light
##--------------------------------------------------

   require(tidybayes)
#   
# winter.fvfm.contrast<-model_winter_fvfm %>% 
#   emmeans::emmeans(~ Treatment,
#               epred = TRUE)%>%emmeans::contrast(method = "revpairwise")%>% 
#   gather_emmeans_draws()#%>% median_hdi()

##Site effects
winter.fvfm.site.contrast=
  model_winter_fvfm %>% 
  emmeans::emmeans(~ Site,
                   epred = TRUE)%>%emmeans::contrast(method = "revpairwise")%>% 
  gather_emmeans_draws()

## What prop is below zero?
winter.fvfm.site.prop=winter.fvfm.site.contrast%>%group_by(contrast)%>%count(.value<0)
winter.fvfm.site.prop$prop=NA
winter.fvfm.site.prop$prop[winter.fvfm.site.prop$`.value < 0`=="TRUE"]= (winter.fvfm.site.prop$n/8000)[winter.fvfm.site.prop$`.value < 0`=="TRUE"]
winter.fvfm.site.prop%<>%filter(`.value < 0`=="TRUE")

winter.fvfm.site.contrast%<>% median_hdi()
winter.fvfm.site.contrast$model="winter.FvFm"

winter.fvfm.site.contrast%<>%left_join(winter.fvfm.site.prop)

winter.fvfm.site.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop))


##---------------------------------------------------------------------------------
##compute posterior distribution of genotype-level response -predicted for new data
##---------------------------------------------------------------------------------

##Create New data  and Predict
gen<- Winterdata |>
  dplyr::distinct(Temp, Light, Site, Genotype,Treatment) 

##Predict Marginal effects for Genotypes
predict.winter.gen.fvfm <- brms::posterior_epred(##already in response scale
  model_winter_fvfm, re_formula = ~(1|Genotype),
  allow_new_levels=TRUE,
  sample_new_levels="uncertainty",#NULL,# predict on genotypes
  newdata = gen
) |>
  t() %>%
  cbind(gen, .)


Winter.fvfm.Genotype.contrast=
  model_winter_fvfm %>%
  emmeans::emmeans(~ Treatment+Genotype,
                   at = list(Genotype = "G00", Site="Bundegi"),
                   epred = TRUE, re_formula=NULL,
                   allow_new_levels = TRUE, sample_new_levels = "uncertainty") %>% 
  emmeans::contrast(method = "revpairwise")%>%
  gather_emmeans_draws()

Winter.fvfm.Genotype.contrast%>%median_hdi()


winter.effect<-Winter.fvfm.Genotype.contrast%>%
  filter(contrast %in% c("Low.Heated G00 - Low.Control G00","High.Heated G00 - High.Control G00"))%>%
  ggplot(aes(x = .value, fill=contrast, alpha=0.8)) + 
  ylab("Density")+
  ggdist::stat_halfeye(point_interval = "median_hdi",.width = c(0.8,0.95))+
  xlab("")+
  ylab("Density")+
  scale_color_manual(values = c("#D11450","#E0B336"))+
  scale_fill_manual(values = c("#D11450","#E0B336"))+
  annotate("text", x =-0.03, y = 0.9, label = "High light", size=3, color="#D11450")+
  annotate("text", x =-0.03, y = 0.76, label = "Low light", size=3, color="#E0B336")+
theme_bw()+
  theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=6),text = element_text(size=6), axis.text = element_text(size=6), 
        axis.ticks.y = element_blank())


## What prop is below zero?
winter.genfvfm.prop=Winter.fvfm.Genotype.contrast%>%group_by(contrast)%>%count(.value<0)
winter.genfvfm.prop$prop=NA
winter.genfvfm.prop$prop[winter.genfvfm.prop$`.value < 0`=="TRUE"]= (winter.genfvfm.prop$n/8000)[winter.genfvfm.prop$`.value < 0`=="TRUE"]
winter.genfvfm.prop%<>%filter(`.value < 0`=="TRUE")


Winter.fvfm.Genotype.contrast%<>% median_hdi()
Winter.fvfm.Genotype.contrast$model="winter.FvFm"


Winter.fvfm.Genotype.contrast%<>%left_join(winter.genfvfm.prop[c(1,4)])
Winter.fvfm.Genotype.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(Winter.fvfm.Genotype.contrast)


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

model_summer_fvfm<-brm(bf(FvFm~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), zi~1), 
                     data=Summerdata, family = zero_inflated_beta(link = "logit"), warmup = 1000, iter=3000)##is this best?
ppcheck_sm1<-pp_check(model_summer_fvfm, type = "dens_overlay", ndraws = NULL)
##describe priors
brms::prior_summary(model_summer_fvfm)


##Predict
predict.global.summer.fvfm <- brms::posterior_epred(
  model_summer_fvfm, re_formula = .~ 1,
  newdata = nd#,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)


##Plot predicted effects
plot.summer.fvfm<-predict.global.summer.fvfm%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment)%>%left_join(unique(predict.global.summer.fvfm[c(2:4)]))%>%#droplevels()%>%
  ggplot(aes(y = Treatment, x = effect,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  ggdist::stat_halfeye(aes(y = Light, x = effect, group=Temp),color="grey",.width = c(0.8, 0.95), 
                       point_interval = "median_hdi",alpha=0.6,size=0.1)+
  facet_grid(vars(Light),scales = "free")+xlab("Posterior FvFm")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  scale_x_continuous(trans="logit")+
  geom_text(data=subset(ann_text,Light=="Low"),
            aes(label="Control",y=1.6, x=0.54),colour="gold", size=3, 
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="Low"),
            aes(label="Heated",y=1.2, x=0.54),colour="#F08622", size=3, 
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="High"),
            aes(label="Control",y=1.6, x=0.54),colour="#F08622", size=3,  
            position = position_nudge(0.0))+
  geom_text(data=subset(ann_text,Light=="High"),
            aes(label="Heated",y=1.2, x=0.54),colour="#D11450", size=3, 
            position = position_nudge(0.0))+#, limits = c(0.1, 0.75))+
  theme_bw()+
  theme(legend.position = "none",
        #legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        #legend.text =element_text(size=7),
        text = element_text(size=6), axis.text = element_text(size=6), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())


##--------------------------------------------------
##Estimate contrasts: mean marginal effect of light
##--------------------------------------------------

summer.fvfm.site.contrast=
  model_summer_fvfm %>% 
  emmeans::emmeans(~ Site,
                   epred = TRUE)%>%emmeans::contrast(method = "revpairwise")%>% 
  gather_emmeans_draws()

## What prop is below zero?
summer.fvfm.site.prop=summer.fvfm.site.contrast%>%group_by(contrast)%>%count(.value<0)
summer.fvfm.site.prop$prop=NA
summer.fvfm.site.prop$prop[summer.fvfm.site.prop$`.value < 0`=="TRUE"]= (summer.fvfm.site.prop$n/8000)[summer.fvfm.site.prop$`.value < 0`=="TRUE"]
summer.fvfm.site.prop%<>%filter(`.value < 0`=="TRUE")


summer.fvfm.site.contrast%<>% median_hdi()
summer.fvfm.site.contrast$model="summer.FvFm"

summer.fvfm.site.contrast%<>%left_join(summer.fvfm.site.prop)

summer.fvfm.site.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop))


##---------------------------------------------------------------------------------
##compute posterior distribution of genotype-level response -predicted for new data
##---------------------------------------------------------------------------------

##Predict Marginal effects for Genotypes
predict.summer.gen.fvfm <- brms::posterior_epred(##already in response scale
  model_summer_fvfm, re_formula = ~(1|Genotype),
  allow_new_levels=TRUE,
  sample_new_levels="uncertainty",#NULL,# predict on genotypes
  newdata = gen
) |>
  t() %>%
  cbind(gen, .)

summer.fvfm.Genotype.contrast=
  model_summer_fvfm %>%
  emmeans::emmeans(~ Treatment+Genotype,
                   at = list(Genotype = "G00", Site="Bundegi"),
                   epred = TRUE, re_formula=NULL,
                   allow_new_levels = TRUE, sample_new_levels = "uncertainty") %>% 
  emmeans::contrast(method = "revpairwise")%>%
  gather_emmeans_draws()

summer.fvfm.Genotype.contrast%>%median_hdi()


summer.effect<-summer.fvfm.Genotype.contrast%>%
  filter(contrast %in% c("Low.Heated G00 - Low.Control G00","High.Heated G00 - High.Control G00"))%>%
  ggplot(aes(x = .value, fill=contrast, alpha=0.8)) + 
  ylab("Density")+
  ggdist::stat_halfeye(point_interval = "median_hdi",.width = c(0.8,0.95))+
  #xlab("Average marginal effect of Light on Heat stress")+
  ylab("Density")+
  xlab("")+
  scale_color_manual(values = c("#D11450","#E0B336"))+
  scale_fill_manual(values = c("#D11450","#E0B336"))+
  annotate("text", x =-0.07, y = 0.93, label = "High light", size=3, color="#D11450")+
  annotate("text", x =-0.07, y = 0.80, label = "Low light", size=3, color="#E0B336")+
  theme_bw()+
  theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=6),text = element_text(size=6), axis.text = element_text(size=6), 
        axis.ticks.y = element_blank())


## What prop is below zero?
summer.genfvfm.prop=summer.fvfm.Genotype.contrast%>%group_by(contrast)%>%count(.value<0)
summer.genfvfm.prop$prop=NA
summer.genfvfm.prop$prop[summer.genfvfm.prop$`.value < 0`=="TRUE"]= (summer.genfvfm.prop$n/8000)[summer.genfvfm.prop$`.value < 0`=="TRUE"]
summer.genfvfm.prop%<>%filter(`.value < 0`=="TRUE")


summer.fvfm.Genotype.contrast%<>% median_hdi()
summer.fvfm.Genotype.contrast$model="summer.FvFm"


summer.fvfm.Genotype.contrast%<>%left_join(summer.genfvfm.prop[c(1,4)])
summer.fvfm.Genotype.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(summer.fvfm.Genotype.contrast)


FvFm.Contrasts=rbind(winter.fvfm.site.contrast,Winter.fvfm.Genotype.contrast,
                     summer.fvfm.site.contrast,summer.fvfm.Genotype.contrast)


###########################
##RedChannel Response 
###########################

##PLOT  RAW data 

Summerdata$Season="Summer"
Winterdata$Season="Winter"

alldata= rbind(Summerdata,Winterdata)

alldata$SeasonTemp=paste(alldata$Season,alldata$Temp, sep=".")
# 
# ggplot(alldata, aes(y=Light, x=RedChannel*255, fill=SeasonTemp))+
#   geom_density_ridges(rel_min_height = 0.005,scale = 1.2)+
#   scale_fill_discrete_sequential(palette="Heat2", aesthetics = c("color","fill"), alpha=0.9)+#+facet_grid(~Site)
#   xlab("Red Channel intensity")+ylab("Light")+#facet_wrap(~Site)+
#   theme(legend.position = "bottom",text = element_text(size=6))+
#   scale_x_continuous(breaks = c(0,50,100,150,200,250),label = c(0,50,100,150,200,250))
# 
# #ggsave(file=paste(myplots,"Experiments/Seasonal_Redchannel.jpeg",sep="/"),width = 5, height = 3.5, units = c("in"),dpi = 300)

##***************************************************
## Bayesian Generalised Random-effects Model -- brms
##***************************************************

model_winter_color<-brm(RedChannel~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), data=Winterdata, family = Beta(),warmup = 1000, iter=3000)
mwc_check<-pp_check(model_winter_color, ndraws=NULL)
##describe priors
brms::prior_summary(model_winter_color)

## Predict 
winter_global_color_preds <- brms::posterior_linpred(
  model_winter_color, re_formula = .~ 1,
  newdata = nd,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)# 


plot.winter.color.preds<-winter_global_color_preds%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment, group=Site)%>%left_join(winter_global_color_preds[c(2:4)])%>%
  ggplot(aes(y = Treatment, x = effect*255,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
   ggdist::stat_halfeye(aes(y = Light, x = effect*255),color="grey",.width = c(0.8, 0.95), point_interval = "median_hdi",alpha=0.6,size=2)+
  facet_grid(vars(Light),scales="free")+xlab("Posterior Red channel intensity")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  theme_bw()+
  theme(legend.position = "none",#legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        #legend.text =element_text(size=7),
        text = element_text(size=6), axis.text = element_text(size=6), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

##-------------------
##Estimate contrasts
##-------------------

winter.color.site.contrast=
  model_winter_color %>% 
  emmeans::emmeans(~ Site,
                   epred = TRUE)%>%emmeans::contrast(method = "revpairwise")%>% 
  gather_emmeans_draws()

winter.color.site.contrast%>%median_hdi()
## What prop is ABOVE zero?
winter.sitecolor.prop=winter.color.site.contrast%>%group_by(contrast)%>%count(.value>0)
winter.sitecolor.prop$prop=NA
winter.sitecolor.prop$prop[winter.sitecolor.prop$`.value > 0`=="TRUE"]= (winter.sitecolor.prop$n/8000)[winter.sitecolor.prop$`.value > 0`=="TRUE"]
winter.sitecolor.prop%<>%filter(`.value > 0`=="TRUE")


winter.color.site.contrast%<>% median_hdi()
winter.color.site.contrast$model="winter.color"


winter.color.site.contrast%<>%left_join(winter.sitecolor.prop[c(1,4)])
winter.color.site.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(winter.color.site.contrast)



##Predict Marginal effects for Genotypes
predict.winter.gen.color <- brms::posterior_epred(##already in response scale
  model_winter_color, re_formula = ~(1|Genotype),
  allow_new_levels=TRUE,
  sample_new_levels="uncertainty",#NULL,# predict on genotypes
  newdata = gen
) |>
  t() %>%
  cbind(gen, .)


Winter.color.Genotype.contrast=
  model_winter_color %>%
  emmeans::emmeans(~ Treatment+Genotype,
                   at = list(Genotype = "G00", Site="Bundegi"),
                   epred = TRUE, re_formula=NULL,
                   allow_new_levels = TRUE, sample_new_levels = "uncertainty") %>% 
  emmeans::contrast(method = "revpairwise")%>%
  gather_emmeans_draws()

Winter.color.Genotype.contrast%>%median_hdi()


winter.color.effect<-Winter.color.Genotype.contrast%>%
  filter(contrast %in% c("Low.Heated G00 - Low.Control G00","High.Heated G00 - High.Control G00"))%>%
  ggplot(aes(x = .value, fill=contrast, alpha=0.8)) + 
  ylab("Density")+
  ggdist::stat_halfeye(point_interval = "median_hdi",.width = c(0.8,0.95))+
  xlab("")+
  ylab("Density")+
  scale_color_manual(values = c("#D11450","#E0B336"))+
  scale_fill_manual(values = c("#D11450","#E0B336"))+
  theme_bw()+
  theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=7),text = element_text(size=6), axis.text = element_text(size=6), 
        axis.ticks.y = element_blank())

## What prop is ABOVE zero?- an effect results in grater red chanel values 
winter.gencolor.prop=Winter.color.Genotype.contrast%>%group_by(contrast)%>%count(.value>0)
winter.gencolor.prop$prop=NA
winter.gencolor.prop$prop[winter.gencolor.prop$`.value > 0`=="TRUE"]= (winter.gencolor.prop$n/8000)[winter.gencolor.prop$`.value > 0`=="TRUE"]
winter.gencolor.prop%<>%filter(`.value > 0`=="TRUE")


Winter.color.Genotype.contrast%<>% median_hdi()
Winter.color.Genotype.contrast$model="winter.color"


Winter.color.Genotype.contrast%<>%left_join(winter.gencolor.prop[c(1,4)])
Winter.color.Genotype.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(Winter.color.Genotype.contrast)


########### SUMMER DATA###################

model_summer_color<-brm(RedChannel~-1+Treatment+Site+(1|Tank)+(Temp * Light||gr(Genotype, by = Site)), 
                        data=Summerdata, family = Beta(),warmup = 1000, iter=3000)
model_summer_color
msccheck<-pp_check(model_summer_color, ndraws = NULL)##OK
##describe priors
brms::prior_summary(model_summer_color)


## Predict 
summer_global_color_preds <- brms::posterior_linpred(
  model_summer_color, re_formula = .~ 1,
  newdata = nd,transform = TRUE
) |>
  t() %>%
  cbind(nd, .)# 


plot.summer.color.preds<-summer_global_color_preds%>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Treatment, group=Site)%>%left_join(summer_global_color_preds[c(2:4)])%>%
  ggplot(aes(y = Treatment, x = effect*255,color=Treatment,fill=Treatment)) + 
  ylab("Light treatment")+
  ggdist::stat_halfeye(aes(y = Light, x = effect*255),color="grey",.width = c(0.8, 0.95), point_interval = "median_hdi",alpha=0.6,size=2)+
  facet_grid(vars(Light),scales="free")+xlab("Posterior Red channel intensity")+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size=6), axis.text = element_text(size=6), 
        axis.text.y.left = element_blank(), axis.ticks.y = element_blank())

##-------------------
##Estimate contrasts
##-------------------

summer.color.site.contrast=
  model_summer_color %>% 
  emmeans::emmeans(~ Site,
                   epred = TRUE)%>%emmeans::contrast(method = "revpairwise")%>% 
  gather_emmeans_draws()

summer.color.site.contrast%>%median_hdi()
## What prop is ABOVE zero?
summer.sitecolor.prop=summer.color.site.contrast%>%group_by(contrast)%>%count(.value>0)
summer.sitecolor.prop$prop=NA
summer.sitecolor.prop$prop[summer.sitecolor.prop$`.value > 0`=="TRUE"]= (summer.sitecolor.prop$n/8000)[summer.sitecolor.prop$`.value > 0`=="TRUE"]
summer.sitecolor.prop%<>%filter(`.value > 0`=="TRUE")
summer.sitecolor.prop$prop[is.na(summer.sitecolor.prop$prop)]<-0


summer.color.site.contrast%<>% median_hdi()
summer.color.site.contrast$model="summer.color"


summer.color.site.contrast%<>%left_join(summer.sitecolor.prop[c(1,4)])
summer.color.site.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(summer.color.site.contrast)



##Predict Marginal effects for Genotypes
predict.summer.gen.color <- brms::posterior_epred(##already in response scale
  model_summer_color, 
  re_formula = ~(1|Genotype),
  #re_formula =~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  allow_new_levels=TRUE,
  sample_new_levels="uncertainty",
  newdata = gen
) |>
  t() %>%
  cbind(gen, .)


summer.color.Genotype.contrast=
  model_summer_color %>%
  emmeans::emmeans(~ Treatment+Genotype,
                   at = list(Genotype = "G00", Site="Bundegi"),
                   epred = TRUE, re_formula=NULL,
                   allow_new_levels = TRUE, sample_new_levels = "uncertainty") %>% 
  emmeans::contrast(method = "revpairwise")%>%
  gather_emmeans_draws()

summer.color.Genotype.contrast%>%median_hdi()


summer.color.effect<-summer.color.Genotype.contrast%>%
  filter(contrast %in% c("Low.Heated G00 - Low.Control G00","High.Heated G00 - High.Control G00"))%>%
  ggplot(aes(x = .value, fill=contrast, alpha=0.8)) + 
  ylab("Density")+
  ggdist::stat_halfeye(point_interval = "median_hdi",.width = c(0.8,0.95))+
  ylab("Density")+
  xlab("")+
  scale_color_manual(values = c("#D11450","#E0B336"))+
  scale_fill_manual(values = c("#D11450","#E0B336"))+
  theme_bw()+
  theme(legend.position = "none",legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'),
        legend.text =element_text(size=7),text = element_text(size=6), axis.text = element_text(size=6), 
        axis.ticks.y = element_blank())


## What prop is ABOVE zero?
summer.gencolor.prop=summer.color.Genotype.contrast%>%group_by(contrast)%>%count(.value>0)
summer.gencolor.prop$prop=NA
summer.gencolor.prop$prop[summer.gencolor.prop$`.value > 0`=="TRUE"]= (summer.gencolor.prop$n/8000)[summer.gencolor.prop$`.value > 0`=="TRUE"]
summer.gencolor.prop%<>%filter(`.value > 0`=="TRUE")


summer.color.Genotype.contrast%<>% median_hdi()
summer.color.Genotype.contrast$model="summer.color"

summer.color.Genotype.contrast%<>%left_join(summer.gencolor.prop[c(1,4)])
summer.color.Genotype.contrast%<>%select(c(model,contrast,.value,.lower,.upper,prop));head(summer.color.Genotype.contrast)
summer.color.Genotype.contrast$prop[is.na(summer.color.Genotype.contrast$prop)]<-0



Color.Contrasts=rbind(winter.color.site.contrast,Winter.color.Genotype.contrast,
                     summer.color.site.contrast,summer.color.Genotype.contrast)

##---------------
##Update Table
##---------------

Table.Contrasts= rbind(FvFm.Contrasts,Color.Contrasts) 
Table.Contrasts$contrast=gsub("G00","",Table.Contrasts$contrast)

write.table(Table.Contrasts, file=paste(mywd,"/Contrasts_table.csv", sep=""), sep=",",quote=FALSE, row.names=T, col.names=T)


##------------------------------------------------------------------------------
##compute posterior draws of genotype-level effects predicted for new data
##------------------------------------------------------------------------------
##----------------------------------------------------------------------------
##Calculate magnitude of predicted shift and compared with change in FvFm response
##----------------------------------------------------------------------------

predict.winter.gen.color2 <- brms::posterior_linpred(
  model_winter_color, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen,transform = TRUE
) |>
  t() %>%
  cbind(gen, .)

##WINTER
winter.color.gen.preds<-predict.winter.gen.color2%>%tidyr::pivot_longer(cols=6:4005,names_to = "Iteration", values_to = "effect")%>%
  group_by(Genotype,Treatment)%>%
  tidybayes::mean_qi(Genotype_mean = effect*255)#, .width = c(.95, .66)) 

##Relative Change in Color
winter.pred.color<-winter.color.gen.preds%>%select(Genotype, Treatment, Genotype_mean)%>%
  spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

gensite<-unique(Winterdata[c(1,2)])##Add site based on GenotypeID
winter.pred.color%<>%left_join(gensite)

winter.pred.color$color="#440154FF"##Bundegi
winter.pred.color$color[winter.pred.color$Site=="Tantabiddi"]<-"#43BF71FF"

##Order Relative Change in Color
descolor<-winter.pred.color%>%arrange(Change)%>%filter(Treatment=="High light")%>%select(color,Genotype)
myorder=forcats::fct_reorder(winter.pred.color$Genotype,winter.pred.color$Change)
descolor$Genotype=factor(descolor$Genotype,levels = levels(myorder))
descolor%<>%arrange(Genotype)


winter.color.Rchange<-winter.pred.color%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,Change), fill=Treatment))+
  xlab("Relative Change in Color\n (Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#F5970A","#FADD3C"))+##Order of factors
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor$color
    glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  geom_text(x=-0.05, y="G02", label="Low light", color="gold", size=2)+
  geom_text(x=-0.05, y="G10", label="High light", color="#F5970A", size=2)+
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "bottom", legend.title = element_blank(),text = element_text(size=6),legend.text = element_text(size=6))

##Predicted FvFm
predict.winter.fvfm.change <- brms::posterior_epred(
  model_winter_fvfm, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen#, transform = TRUE
) |>
  t() %>%
  cbind(gen, .)

predict.winter.fvfm.change%<>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Genotype,Treatment)%>%
  tidybayes::mean_qi(Genotype_mean = effect)

# Relative Change 
winter.R.fvfm<-predict.winter.fvfm.change%>%select(Genotype, Treatment, Genotype_mean)%>%spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

##Order treatment levels
winter.R.fvfm$Treatment<-factor(winter.R.fvfm$Treatment,levels=c("Low light","High light"))

W.fvfm.change<-winter.R.fvfm%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,wintercolor$Change), fill=Treatment))+
  xlab("Relative Change in FvFm \n(Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#FADD3C", "#F5970A"))+
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor$color
    glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "bottom", legend.title = element_blank(),text = element_text(size=6),legend.text = element_text(size=6))+
  geom_text(x=min(winter.R.fvfm$Change)+0.01, y="G13", label="Bundegi", color="#440154FF", size=2)+
  geom_text(x=min(winter.R.fvfm$Change)+0.01, y="G16", label="Tantabiddi", color="#43BF71FF", size=2)



##---------
##SUMMER

predict.summer.gen.color2 <- brms::posterior_linpred(
  model_summer_color, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen,transform = TRUE
) |>
  t() %>%
  cbind(gen, .)

summer.color.gen.preds<-predict.summer.gen.color2%>%tidyr::pivot_longer(cols=6:4005,names_to = "Iteration", values_to = "effect")%>%
group_by(Genotype,Treatment)%>%
tidybayes::mean_qi(Genotype_mean = effect*255)#, .width = c(.95, .66)) 

##Relative Change in Color
summer.pred.color<-summer.color.gen.preds%>%select(Genotype, Treatment, Genotype_mean)%>%spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

gensite<-unique(Summerdata[c(1,2)])##Add site based on GenotypeID
summer.pred.color%<>%left_join(gensite)

summer.pred.color$color="#440154FF"##Bundegi
summer.pred.color$color[summercolor$Site=="Tantabiddi"]<-"#43BF71FF"

##Order Relative Change in Color
descolor<-summer.pred.color%>%arrange(Change)%>%filter(Treatment=="High light")%>%select(color,Genotype)
myorder=forcats::fct_reorder(summer.pred.color$Genotype,summer.pred.color$Change)
descolor$Genotype=factor(descolor$Genotype,levels = levels(myorder))
descolor%<>%arrange(Genotype)


summer.color.Rchange<-summer.pred.color%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,Change), fill=Treatment))+
  xlab("Relative Change in Color\n (Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#F5970A","#FADD3C"))+##Order of factors
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor$color
  glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  geom_text(x=0.4, y="G20", label="Low light", color="gold", size=2)+
  geom_text(x=0.4, y="G11", label="High light", color="#F5970A", size=2)+
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "none", legend.title = element_blank(),text = element_text(size=6),legend.text = element_text(size=6))
  
##Predicted FvFm
predict.summer.fvfm.change <- brms::posterior_epred(
  model_summer_fvfm, re_formula = .~ 1 + (Temp * Light||gr(Genotype, by = Site)),
  newdata = gen#, transform = TRUE
) |>
  t() %>%
  cbind(gen, .)

predict.summer.fvfm.change%<>%dplyr::select(-c(Temp,Light))%>%
  tidyr::pivot_longer(cols=4:4003,names_to = "Iteration", values_to = "effect")%>%
  group_by(Genotype,Treatment)%>%
  tidybayes::mean_qi(Genotype_mean = effect)

# Relative Change 
S.fvfm.change<-predict.summer.fvfm.change%>%select(Genotype, Treatment, Genotype_mean)%>%spread(Treatment,Genotype_mean)%>%
  mutate(`Low light`=((Low.Heated-Low.Control)/Low.Control))%>%
  mutate(`High light`=((High.Heated-High.Control)/High.Control))%>%select(Genotype,`Low light`,`High light`)%>%
  gather(key="Treatment",value="Change",`Low light`,`High light`)%>%droplevels()

##Order treatment levels
S.fvfm.change$Treatment<-factor(S.fvfm.change$Treatment,levels=c("Low light","High light"))

s.fvfm.Rchange<-S.fvfm.change%>%ggplot(aes(x=Change, y=forcats::fct_reorder(Genotype,summercolor$Change), fill=Treatment))+
  xlab("Relative Change in FvFm \n(Posterior mean)")+ylab("Genotype")+
  geom_point(aes(colour=Treatment),size=1.5, position=position_dodge(width = 0.3))+
  scale_color_manual(values=c("#FADD3C", "#F5970A"))+
  geom_linerange(aes(xmin=0, xmax=Change, y=Genotype, yend=Genotype,colour=Treatment),position=position_dodge(width = 0.3))+
  geom_vline(xintercept=-0.5,col = "gray60", linetype="dashed")+
  theme_bw()+scale_y_discrete(labels = function(x) {
    colours <-descolor$color
    glue::glue("<span style='color:{colours}'>{x}</span>")
  }) +
  theme(axis.text.y = ggtext::element_markdown())+
  theme(legend.position = "bottom", legend.title = element_blank(),text = element_text(size=6),legend.text = element_text(size=6))+
  geom_text(x=-0.4, y="G18", label="Bundegi", color="#440154FF", size=2)+
  geom_text(x=-0.4, y="G10", label="Tantabiddi", color="#43BF71FF", size=2)


##-------------------------------------
### Correlate physiological parameters
##-------------------------------------

head(alldata)

require(ggpmisc)


cor_phenotype<-alldata%>%filter(FvFm>0 & Season=="Summer")%>%
  ggplot(aes(y=FvFm, x=RedChannel,grp.label = Season))+
  geom_point(aes(colour=Temp), alpha=0.3)+
  geom_smooth(method = "lm",color=c("grey40"))+
  stat_poly_eq(use_label(c("grp","adj.R2")), label.y = c(0.95, 0.95), label.x = c(0.1,0.9), size=2)+
  labs(x="Red channel intensity", y= "Fv/Fm")+
  scale_color_manual(values = rep(c("#8491B4B2","#D11450"),10))+
  annotate("text", x =0.60, y = 0.38, label = "Control", size=2, color="#8491B4B2")+
  annotate("text", x =0.60, y = 0.33, label = "Heated", size=2, color="#D11450")+
  theme_bw()+
  theme(legend.position = "none",text = element_text(size=6), axis.text = element_text(size=6))


alldata$Season=factor(alldata$Season, levels=c("Winter","Summer"))
gencolors<-ggsci::pal_d3("category20b")(20)

##Genotype-level  
summer_cor_phenotype<-alldata%>%filter(FvFm>0 & Season=="Summer")%>%
  ggplot(aes(y=FvFm, x=RedChannel, color = Genotype))+
  geom_point(aes(colour=Genotype), alpha=0.3)+
  geom_smooth(method = "lm")+
ggsci::scale_color_d3("category20b")+
    labs(x="Red channel intensity", y= "Fv/Fm")+
  theme_bw()+
  annotate("text", x =0.55, y = 0.45, label = "G01", size=2, color=gencolors[1])+
  annotate("text", x =0.55, y = 0.42, label = "G02", size=2, color=gencolors[2])+
  annotate("text", x =0.55, y = 0.39, label = "G03", size=2, color=gencolors[3])+
  annotate("text", x =0.55, y = 0.36, label = "G04", size=2, color=gencolors[4])+
  annotate("text", x =0.55, y = 0.33, label = "G05", size=2, color=gencolors[5])+
  annotate("text", x =0.55, y = 0.30, label = "G06", size=2, color=gencolors[6])+
  annotate("text", x =0.55, y = 0.27, label = "G07", size=2, color=gencolors[7])+
  annotate("text", x =0.55, y = 0.24, label = "G08", size=2, color=gencolors[8])+
  annotate("text", x =0.55, y = 0.21, label = "G09", size=2, color=gencolors[9])+
  annotate("text", x =0.55, y = 0.18, label = "G10", size=2, color=gencolors[10])+
  annotate("text", x =0.62, y = 0.45, label = "G11", size=2, color=gencolors[11])+
  annotate("text", x =0.62, y = 0.42, label = "G12", size=2, color=gencolors[12])+
  annotate("text", x =0.62, y = 0.39, label = "G13", size=2, color=gencolors[13])+
  annotate("text", x =0.62, y = 0.36, label = "G14", size=2, color=gencolors[14])+
  annotate("text", x =0.62, y = 0.33, label = "G15", size=2, color=gencolors[15])+
  annotate("text", x =0.62, y = 0.30, label = "G16", size=2, color=gencolors[16])+
  annotate("text", x =0.62, y = 0.27, label = "G17", size=2, color=gencolors[17])+
  annotate("text", x =0.62, y = 0.24, label = "G18", size=2, color=gencolors[18])+
  annotate("text", x =0.62, y = 0.21, label = "G19", size=2, color=gencolors[19])+
  annotate("text", x =0.62, y = 0.18, label = "G20", size=2, color=gencolors[20])+
  theme(legend.position = "none",text = element_text(size=6), axis.text = element_text(size=6))

##-----------------------
### FINAL PLOTS ######
##-----------------------
sitecolor=scales::viridis_pal(end=0.7)(length(unique(Winterdata$Site)))
strip <- strip_themed(background_x = elem_list_rect(fill = alpha(sitecolor, 0.3)))

#plot observed response per treatment
w1<-Winterdata%>%group_by(Site,Tank,Treatment)%>%summarise(FvFm=mean(FvFm))%>%
  ggplot(aes(x=Treatment, y=FvFm, fill=Treatment))+
  ylab("Photochemical efficiency (Fv/Fm)")+
  scale_x_discrete(labels=c('Low\nControl', 'Low\nHeated', 'High\nControl', 'High\nHeated'))+
  xlab("")+
  geom_boxplot()+
  scale_fill_discrete_sequential(palette="Heat2",labels=c("Low.Control (LC)","Low.Heated (LH)","High.Control (HC)", "High.Heated (HH)"))+
  facet_grid2(~Site, strip=strip)+theme_bw()+
  theme(legend.position = "none",legend.key.size = unit(0.5,"cm"),text = element_text(size=6))#, panel.background = element_rect(fill="transparent"))

##Add text to plot
treats=c("Low.Control","Low.Heated","High.Control","High.Heated")
sdata<-data.frame(Treatment=treats,
             Site=rep("Bundegi",4),
             FvFm=rep(0.6,4))

##raw points for average per treatments
w1a<-w1+ggtitle("Observed Fv/Fm response")+
  geom_point(data=Winterdata%>%group_by(Site,Genotype,Treatment)%>%
               summarise(FvFm=mean(FvFm)), aes(x=Treatment, y=FvFm, fill=Treatment), alpha=0.2, size=0.3)+
  geom_jitter(alpha=0.1, size=0.3)+
  ylim(c(0.6, 0.75))

##Color
w1rc<-Winterdata%>%group_by(Site,Tank,Treatment)%>%summarise(RedChannel=mean(RedChannel)*255)%>%
  ggplot(aes(x=Treatment, y=RedChannel, fill=Treatment))+
  ylab("Red channel intensity")+
  scale_x_discrete(labels=c('Low\nControl', 'Low\nHeated', 'High\nControl', 'High\nHeated'))+
  xlab("")+
  geom_boxplot()+
  scale_fill_discrete_sequential(palette="Heat2",labels=c("Low.Control (LC)","Low.Heated (LH)","High.Control (HC)", "High.Heated (HH)"))+
  facet_grid2(~Site, strip=strip)+
  theme_bw()+
  ggtitle("Observed color response")+
  geom_point(data=Winterdata%>%group_by(Site,Genotype,Treatment)%>%
  summarise(RedChannel=mean(RedChannel)*255), aes(x=Treatment, y=RedChannel, fill=Treatment), alpha=0.2, size=0.3)+
  geom_jitter(alpha=0.1, size=0.3)+
  theme(legend.position = "none",legend.key.size = unit(0.5,"cm"),text = element_text(size=6),
  axis.text.x = element_text(margin = margin(0, 0, 0, 0), size=5))


##---------------------------------------
##Figure 4- Winter results
##---------------------------------------

ggarrange(ncol=2, nrow=3, w1a,w1rc, 
          plot.winter.fvfm+ggtitle("Predicted Fv/Fm response"),
          plot.winter.color.preds+ggtitle("Predicted color response")+ylab(""),
          winter.effect+ggtitle("Marginal effect on Fv/Fm stress response"),
          winter.color.effect+ggtitle("Marginal effect on color stress response")+ylab(""),
          labels = c("a","","b","","c","d"), font.label = list(size=12, face="bold"))%>%
ggsave(file=paste(myplots,"Fig 4.jpg", sep="/"),width = 5.4, height = 6.2, units = c("in"),dpi = 500)
    
##---------------------------------------
##Figure 5- Summer results
##---------------------------------------

#Observed FvFm
Sa<-Summerdata%>%group_by(Site,Tank,Treatment)%>%summarise(FvFm=mean(FvFm))%>%
  ggplot(aes(x=Treatment, y=FvFm, fill=Treatment))+
  ylab("Photochemical efficiency (Fv/Fm)")+
  scale_x_discrete(labels=c('Low\nControl', 'Low\nHeated', 'High\nControl', 'High\nHeated'))+
  xlab("")+
  geom_boxplot()+
  scale_fill_discrete_sequential(palette="Heat2",labels=c("Low.Control (LC)","Low.Heated (LH)","High.Control (HC)", "High.Heated (HH)"))+
  facet_grid2(~Site, strip=strip)+theme_bw()+
  theme(legend.position = "none",legend.key.size = unit(0.5,"cm"),text = element_text(size=6))+
  ggtitle("Observed Fv/Fm response")+
  geom_point(data=Summerdata%>%group_by(Site,Genotype,Treatment)%>%
               summarise(FvFm=mean(FvFm)), aes(x=Treatment, y=FvFm, fill=Treatment), alpha=0.2, size=0.3)+
  geom_jitter(alpha=0.1, size=0.3)

##Color
Sb<-Summerdata%>%group_by(Site,Tank,Treatment)%>%summarise(RedChannel=mean(RedChannel)*255)%>%
  ggplot(aes(x=Treatment, y=RedChannel, fill=Treatment))+
  ylab("Red channel intensity")+
  scale_x_discrete(labels=c('Low\nControl', 'Low\nHeated', 'High\nControl', 'High\nHeated'))+
  xlab("")+
  geom_boxplot()+
  scale_fill_discrete_sequential(palette="Heat2",labels=c("Low.Control (LC)","Low.Heated (LH)","High.Control (HC)", "High.Heated (HH)"))+
  facet_grid2(~Site, strip=strip)+
  theme_bw()+
  ggtitle("Observed color response")+
  geom_point(data=Summerdata%>%group_by(Site,Genotype,Treatment)%>%
               summarise(RedChannel=mean(RedChannel)*255), aes(x=Treatment, y=RedChannel, fill=Treatment), alpha=0.2, size=0.3)+
  geom_jitter(alpha=0.1, size=0.3)+
  theme(legend.position = "none",legend.key.size = unit(0.5,"cm"),text = element_text(size=6),
        axis.text.x = element_text(margin = margin(0, 0, 0, 0), size=5))


ggarrange(ncol=2, nrow=3, Sa,Sb, 
          plot.summer.fvfm+ggtitle("Predicted Fv/Fm response"),
          plot.summer.color.preds+ggtitle("Predicted color response")+ylab(""),
          summer.effect+ggtitle("Marginal effect on Fv/Fm stress response"),
          summer.color.effect+ggtitle("Marginal effect on color stress response")+ylab(""),
          labels = c("a","","b","","c","d"), font.label = list(size=12, face="bold"))%>%
  ggsave(file=paste(myplots,"Fig 5.jpg", sep="/"),width = 5.4, height = 6.2, units = c("in"),dpi = 500)

##---------------------------------------
##Figure 6- Summer results
##---------------------------------------

###Summer only
ggarrange(s.fvfm.Rchange,summer.color.Rchange+ylab("")+theme(axis.text.y = element_blank()), 
          common.legend = F, legend = "none",
          labels = c("a","b","c",""), font.label = list(size=12, face="bold"),
          align = "hv")%>%
  ggsave(filename =paste(myplots, "Figure 6.png",sep = "/"), dpi=500, width = 3.5, height = 5)


ggarrange(s.fvfm.Rchange+theme(legend.position = "none"),
          summer.color.Rchange+ylab("")+theme(axis.text.y = element_blank()), 
          cor_phenotype,
          summer_cor_phenotype+ylab(""),
          labels = c("a","b","c",""), font.label = list(size=12, face="bold"),
          align = "hv")%>%
  ggsave(file=paste(myplots,"Figure 6v2.jpg", sep="/"),
         width = 4, height = 6.2, units = c("in"),dpi = 500)


## Winter
ggarrange(W.fvfm.change,winter.color.Rchange+ylab("")+theme(axis.text.y = element_blank()),
          common.legend = F, legend = "none",labels = "auto", align = "hv")%>%
  ggsave(filename =paste(myplots, "Figure S4.png",sep = "/"), dpi=400, width = 4, height = 5)


#save(Summerdata, Winterdata, summercolor,summerfvfm,wintercolor,winfvfm,file="Experiments_Tables.R")
save.image("bmrs model outputs.R")

