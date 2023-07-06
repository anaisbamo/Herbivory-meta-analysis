## Title: Herbivory meta-analysis
## Author: Anaïs Barrera Montoro
## Contact: anais.barrera@ceab.csic.es
## Date: 04/04/2023



##1. Loading packages

library(ggplot2)
library(maptools)
library(rgdal)
library(metafor)
library(readxl)
library(easypackages)
libraries("ggplot2","patchwork","ggthemes", "tidyverse", 
          "dplyr", "ggsci", "tidyr", "metafor", "readr", 
          "readxl", "ggpubr", "performance", "car", "MuMIn", 
          "lmtest", "sandwich", "DHARMa", "effects")



predictvals <- function(model, xvar, yvar, xrange = NULL, samples = 100, ...) {
  
  # If xrange isn't passed in, determine xrange from the models.
  # Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm")))
      xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess"))
      xrange <- range(model$x)
  }
  
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  newdata
}



##2. Set the working directory
setwd("~/Documents/Documents/FEINA/CEAB/PROJECTES/METANALISI/TREBALL FI DE GRAU/TERCERA CERCA/Herbivory-meta-analysis/Data")


##3. Loading the data for temperature (ºC) and consumption rates.

data_T <- read_excel( "TERCERA_CERCA_R.xlsx", sheet = "data_temperature")


##4. Map the location of the studies included

map.data <-data_T[,c(2, 8:12)] 

buffer <- 0.11

geo_bounds <- c(left = min(map.data$Long)-buffer, 
                bottom = min(map.data$Lat)-buffer, 
                right = max(map.data$Long)+buffer, 
                top = max(map.data$Lat)+buffer)

Sites.grid <- expand.grid(lon_bound = c(geo_bounds[1], geo_bounds[3]), 
                          lat_bound = c(geo_bounds[2], geo_bounds[4]))

coordinates(Sites.grid) <- ~ lon_bound + lat_bound

world <- readOGR("~/Documents/Documents/FEINA/CEAB/PROJECTES/METANALISI/TREBALL FI DE GRAU/TERCERA CERCA/Herbivory-meta-analysis1/Scripts/TM_WORLD_BORDERS-0/TM_WORLD_BORDERS-0.3.shp")

plot(world)
ggplot() + 
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill="grey", colour="black") +
  coord_equal() +
  geom_point(data=map.data, aes(x=Long, y=Lat), colour="#4883B4FF") +
  labs(x="Longitude", y="Latitude") +
  theme_classic()



######### 5. Meta-analysis: effect of T on grazing rates ########

# 5.1  Calculate effect size (ES)

## 5.1.1 with the log transformed ratio of means (LRR) and its variance (LRR_var) measure= ROM. 

T_ES <-escalc(measure="ROM",
              m1i=treat_mean,
              m2i=Control_mean,
              sd1i=treat_SD,
              sd2i=control_SD,
              n1i=N_treat,
              n2i=N_control, 
              data=data_T,var.names=c("LRR","LRR_var"),digits=4)

hist(T_ES$LRR, breaks=10)

# Include LRR and LRR_var produced in the dataset and see effect sizes values distribution (histogram)

data_T$EffectSize_LRR <- T_ES$LRR
data_T$ESVariance_LRR <- T_ES$LRR_var

# study number 43 and 47 produce +-infin values for effect size because the log ratio of the means log(m1i/m2i) have opposite signs
# (and some means are 0 = no plant material remaining) so, exclude these studies from data set. ROM measure can
# only be used when both values are positive.
# Any study with Na in the EffectSize_LRR column or the ESVariance_LRR column will be excluded
data_T2 <-data_T %>% filter(! is.na(EffectSize_LRR))
data_T3 <- data_T2 %>% filter(! is.na(ESVariance_LRR))

## 5.1.2 with the standarized mean difference (Hedges' g), measure= SMD.

T_ES_SMD <-escalc(measure="SMD", m1i=treat_mean, m2i=Control_mean, sd1i=treat_SD, sd2i=control_SD, n1i=N_treat, n2i=N_control,  data=data_T3,var.names=c("SMD","SMD_var"),digits=4)


data_T3$EffectSize_SMD <- T_ES_SMD$SMD
data_T3$ESVariance_SMD <- T_ES_SMD$SMD_var

hist(T_ES_SMD$SMD)


# 5.2 Meta-analysis for temperature effect on grazing rates. RANDOM EFFECT MODEL rma() with 
# restricted maximum likelihood estimator (REML).

# Exclude Sampaio et al., 2017
data_T4 <- data_T3 %>% filter(Authors != "Sampaio_et_al.") 
meta_T_random <-rma(yi=EffectSize_LRR, vi=ESVariance_LRR, slab=paste(Authors, Year, sep=", "), method="REML", data=data_T4)
funnel(meta_T_random)
forest(meta_T_random, 
       vi=data_T4$ESVariance_LRR,
       slab = data_T4$Authors,
       xlab = "LRR of consumption rates for Temperature",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)

# 5.2.2 subset for SEAGRASS studies
seagrass <- subset(data_T4, data_T4$Habitat_stucture=="erect_seagrass")
meta_T_randomS <- rma(yi=EffectSize_LRR, vi=ESVariance_LRR, data=seagrass,method="REML")
funnel(meta_T_randomS, xlab = "LRR of consumption rates for Temperature on seagrasses")
forest(meta_T_randomS,
       vi=seagrass$ESVariance_LRR,
       slab = seagrass$Authors,
       xlab = "LRR of consumption rates for Temperature on seagrasses",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)

# 5.2.3 subset for ALGAE studies 
algae <- subset(data_T4, data_T4$Habitat_stucture=="erect_macroalge")
meta_T_randomA <- rma(yi=EffectSize_LRR, vi=ESVariance_LRR, data=algae,method="REML")
funnel(meta_T_randomA)
forest(meta_T_randomA,
       vi=algae$ESVariance_LRR,
       slab = algae$Authors,
       xlab = "LRR of consumption rates for Temperature on algae",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)


# 5.3 Meta-analysis for temperature effect on grazing rates. MIXED-EFFECT MODEL rma.mv() with
# restricted maximum likelihood estimator (REML).

meta_T_mixed <- rma.mv(yi = EffectSize_LRR , V = ESVariance_LRR, 
                       method = "REML", 
                       random = ~1|Authors + specie,
                       level = 95, digits = 4, data = data_T4) 
funnel(meta_T_mixed)
forest(meta_T_mixed,
       vi=data_T4$ESVariance_LRR,
       slab = data_T4$Authors,
       xlab = "Mixed effects model for grazing rates and temperature",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)

# 5.3.1 mixed effect model for seagrasses

meta_T_mixedS <- rma.mv(yi = EffectSize_LRR , V = ESVariance_LRR, 
                        method = "REML", 
                        random = ~1|Authors + specie,
                        level = 95, digits = 4, data = seagrass) 
funnel(meta_T_mixedS)
forest(meta_T_mixedS,
       vi=seagrass$ESVariance_LRR,
       slab = seagrass$Authors,
       xlab = "Mixed effects model for grazing rates and temperature on seagrass",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)

# 5.3.2 mixed effect model for algae

meta_T_mixedA <- rma.mv(yi = EffectSize_LRR , V = ESVariance_LRR, 
                        method = "REML", 
                        random = ~1|Authors + specie,
                        level = 95, digits = 4, data = algae) 
funnel(meta_T_mixedA)
forest(meta_T_mixedA,
       vi=algae$ESVariance_LRR,
       slab = algae$Authors,
       xlab = "Mixed effects model for grazing rates and temperature on algae",
       xlim = c(-15,23),
       showweights = TRUE,
       col = "#84BD00FF", cex=0.8)


# 6.Test between Random and Mixed effects model.

#6.1 for all macrophytes

AIC(meta_T_random, meta_T_mixed)
# meta_T_random  2  167.5051
# meta_T_mixed   3 1520.5925


#6.2 for seagrass
AIC(meta_T_randomS, meta_T_mixedS)

#df      AIC
#meta_T_randomS  2 21.10014
#meta_T_mixedS   3 62.23343

#6.3 for algae

AIC(meta_T_randomA, meta_T_mixedA)
#               df       AIC
# meta_T_randomA  2  148.0821
# meta_T_mixedA   3 1461.3835



# 7. Calculate the % of increment in temperature between control and treatment for each study. 
# Plot the response (effect size as LRR) as a function of the % increment in temp. 
data_T4$perc.temp <- ((data_T4$Temp_mean-data_T4$TempC_mean)/(data_T4$TempC_mean))*100
plot(data_T4$EffectSize_LRR~data_T4$perc.temp)
# 7.1 Do subsets for seagrass and algae again to include % increment in temp. in data_4
algae <- subset(data_T4, data_T4$Habitat_stucture=="erect_macroalge")
seagrass <- subset(data_T4, data_T4$Habitat_stucture=="erect_seagrass")
# plots
plot(seagrass$EffectSize_LRR~seagrass$perc.temp)
plot(algae$EffectSize_LRR~algae$perc.temp)



# 8. plot the response (effect size as LRR) by latitude
plot(data_T4$EffectSize_LRR~data_T4$Lat, xlim=c(-60, 60))
# 8.1 do it for the algae
plot(algae$Lat~algae$EffectSize_LRR, ylim=c(-60, 60))
# 8.2 do it for seagrass
plot(seagrass$Lat~seagrass$EffectSize_LRR, ylim=c(-60, 60))


# 9. Calculate the degrees of difference in temp. between control and treatment
data_T4$dif.temp <- data_T4$Temp_mean-data_T4$TempC_mean
# modeling response of all herbivores to degrees of difference in tmeperature
model_performance_T <- lm(EffectSize_LRR~I(dif.temp^2), data=data_T4)
# plot
ggplot(data_T4, aes(x=dif.temp, y=EffectSize_LRR, colour= herb_group))+ 
  geom_point() + stat_smooth(method=lm, se = TRUE) + xlim(0,15)
labs(title="Effect of temperature on urchins grazing rates",x = "Increase in temperature vs. control (degrees)", y = "Effect size")+
  theme(title = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black" ))+
  theme(axis.title.x = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(axis.title.y = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line= element_line(colour = "black",size = 0.5))+
  theme(axis.text.x = element_text(family = "Helvetica", face = "plain",colour = "black", size = 30))+
  theme(axis.text.y = element_text( size = 30))+
  theme(legend.title = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.text = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.key.size = unit(2, 'cm')) 



# 10. Do subsets for herbivore groups and plot the EffectSize_LRR as a response of the difference in temperature degrees.Also, 
#.   quadratic models for each herb group.

# 10.1 for urchin
urchin <- subset(data_T4, data_T4$herb_group=="urchin")
hist(urchin$EffectSize_LRR)
urchin_model <- lm(EffectSize_LRR~I(dif.temp^2), data=urchin)
urchin_model <- lm(log(EffectSize_LRR+abs(min(EffectSize_LRR))+1)~I(dif.temp^2), data=urchin)
hist(log(urchin$EffectSize_LRR+abs(min(urchin$EffectSize_LRR))+1))
plot(urchin_model)
summary(urchin_model)
#    Estimate Std. Error t value Pr(>|t|)  
#(Intercept)    0.202468   0.225679   0.897   0.3815  
#I(dif.temp^2) -0.009451   0.003434  -2.752   0.0131 *

car::Anova(urchin_model)
##Response: EffectSize_LRR
#Sum Sq Df F value  Pr(>F)  
#I(dif.temp^2) 3.9380  1  7.5736 0.01312 *
#Residuals     9.3593 18                  

library(DHARMa)
plot(urchin_model)
shapiro.test(resid(urchin_model))
hist(resid(urchin_model))
ggplot(urchin, aes(x=dif.temp, y=EffectSize_LRR))+ 
  geom_point(data= urchin, aes(colour=Authors)) + stat_smooth() + xlim(0,14)+
  labs(title="Effect of temperature on urchins grazing rates",x = "Increase in temperature vs. control (degrees)", y = "Effect size")+
  theme(title = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black" ))+
  theme(axis.title.x = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(axis.title.y = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line= element_line(colour = "black",size = 0.5))+
  theme(axis.text.x = element_text(family = "Helvetica", face = "plain",colour = "black", size = 30))+
  theme(axis.text.y = element_text( size = 30))+
  theme(legend.title = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.text = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.key.size = unit(2, 'cm'))

# 10.2 for gastropods
gastropods <- subset(data_T4, data_T4$herb_group=="gastropod")
gastropods_model <- lm(EffectSize_LRR~I(dif.temp^2), data= gastropods)
summary(gastropods_model)
plot(gastropods_model)
ggplot(gastropods, aes(x=dif.temp, y=EffectSize_LRR))+ 
  geom_point(data= gastropods, aes(colour=Authors)) + stat_smooth(method= lm, se=TRUE) + xlim(0,14)+
  labs(title="Effect of temperature on gastropods grazing rates",x = "Increase in temperature vs. control (degrees)", y = "Effect size")+
  ylim(-3,3)+
  theme(title = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black" ))+
  theme(axis.title.x = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(axis.title.y = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line= element_line(colour = "black",size = 0.5))+
  theme(axis.text.x = element_text(family = "Helvetica", face = "plain",colour = "black", size = 30))+
  theme(axis.text.y = element_text( size = 30))




# 10.3 for isopods
isopods <- subset(data_T4, data_T4$herb_group=="isopod")
isopods_model <- lm(EffectSize_LRR~I(dif.temp^2), data= isopods)
summary(isopods_model)
plot(isopods_model)
ggplot(isopods, aes(x=dif.temp, y=EffectSize_LRR))+ 
  geom_point(data= isopods, aes(colour=Authors)) + stat_smooth(method=lm, se = TRUE) + xlim(0,15) + xlim(0,15) +
  labs(title="Effect of temperature on isopods grazing rates",x = "Increase in temperature vs. control (degrees)", y = "Effect size")+
  ylim(-5,5)+
  theme(title = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black" ))+
  theme(axis.title.x = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(axis.title.y = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line= element_line(colour = "black",size = 0.5))+
  theme(axis.text.x = element_text(family = "Helvetica", face = "plain",colour = "black", size = 30))+
  theme(axis.text.y = element_text( size = 30))+
  theme(legend.title = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.text = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.key.size = unit(2, 'cm'))

# 10.4 for amphipods
amphipods <-subset(data_T4, data_T4$herb_group=="amphipod")
amphipods_model <-  lm(EffectSize_LRR~I(dif.temp^2), data= amphipods)
summary(amphipods_model)
plot(amphipods_model)
ggplot(amphipods, aes(x=dif.temp, y=EffectSize_LRR))+ 
  geom_point(data= amphipods, aes(colour=Authors)) + geom_smooth(method="lm", formula= y~x) + xlim(0, 7)+
  labs(title="Effect of temperature on amphipod grazing rates",x = "Increase in temperature vs. control (degrees)", y = "Effect size")+
  theme(title = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black" ))+
  theme(axis.title.x = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(axis.title.y = element_text(size = 30, lineheight = .9,family = "Helvetica", face = "plain", colour = "black"))+
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line= element_line(colour = "black",size = 0.5))+
  theme(axis.text.x = element_text(family = "Helvetica", face = "plain",colour = "black", size = 30))+
  theme(axis.text.y = element_text( size = 30))+
  theme(legend.title = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.text = element_text(face = "plain",family = "Helvetica",colour = "black",size = 20))+
  theme(legend.key.size = unit(2, 'cm'))

# 10.5 plot different herb groups on the same plot fitted to the quadratic models.
library(dplyr)
#quadratic_models <- data_T4 %>%
group_by(herb_group) %>%
  do(model = lm (EffectSize_LRR~I(dif.temp^2), .)) %>%
  ungroup()

#quadratic_models$model

# get predicted values from each model
#predvals_qmodels <- quadratic_models %>%  group_by(herb_group) %>% do(predictvals(.$model[1], xvar="dif.temp", yvar="EffectSize_LRR"))

# plot all quadratic models
#ggplot(data_T4, aes(x=dif.temp, y=EffectSize_LRR, colour= herb_group)) + 
geom_point() + theme_bw() + geom_line(data = predvals_qmodels)



#11. Overall effect size (LRR) against temperature  by temperature variation steps (+- degrees from control) for all studies. 
# build a meta-analysis model for each temperature step and extract the estimates from the model (overall effect as a matrix=b,
# lower confidenceinterval and upper confidenc einterval) and store it as a data frame. 

# plot all the studies according to the temperature "step" in which they're placed. Coloured by the herbivore group.

ggplot(data_T4, aes(x=dif.temp, y=EffectSize_LRR, colour= herb_group))+ 
  geom_point() +  theme_bw() +
  facet_wrap(~dif.temp) + theme(legend.position="bottom")

# create categories for temperature steps
#data_T4_Ov.est <- data_T4 %>% mutate("t_steps"= case_when (dif.temp < 0 ~ "< -4",
#dif.temp >= 0 & dif.temp <= 2 ~ "2",
#dif.temp == 3 ~ "3",
#dif.temp > 3 & dif.temp <= 4 ~ "4",
#dif.temp > 4 & dif.temp <= 6 ~ "5",
#dif.temp == 7 ~ "8",
#dif.temp == 8 ~ "8",
#dif.temp == 10 ~"10",
#dif.temp >=11 ~ ">12"))




#11.1. Model for < -4 temperature step (includes -10, -5 i -4)

Tstep_-4 <- subset()


#overall_estimate <- round(meta_T_mixed$b[[1,1]],4)

data_frame_estimates<- data.frame(estimate=meta_T_mixed$b[[1,1]],
                                  lb= meta_T_mixed$ci.lb,
                                  ub=meta_T_mixed$ci.ub)

































########### 12. Meta-analysis: effect of N on grazing rates ########

##12.1. Loading the data for nutrients and grazing rates.

data_N <- read_excel( "TERCERA_CERCA_R.xlsx", sheet = "data_nutrients")

## Map the location of the studies included

map.data <-data_N[,c(2, 8:12)] 

buffer <- 0.11

geo_bounds <- c(left = min(map.data$Long)-buffer, 
                bottom = min(map.data$Lat)-buffer, 
                right = max(map.data$Long)+buffer, 
                top = max(map.data$Lat)+buffer)

Sites.grid <- expand.grid(lon_bound = c(geo_bounds[1], geo_bounds[3]), 
                          lat_bound = c(geo_bounds[2], geo_bounds[4]))

coordinates(Sites.grid) <- ~ lon_bound + lat_bound

world <- readOGR("~/Documents/Documents/FEINA/CEAB/PROJECTES/METANALISI/TREBALL FI DE GRAU/TERCERA CERCA/Herbivory-meta-analysis1/Scripts/TM_WORLD_BORDERS-0/TM_WORLD_BORDERS-0.3.shp")

plot(world)
ggplot() + 
  geom_polygon(data = world, aes(x=long, y=lat, group=group), fill="grey", colour="black") +
  coord_equal() +
  geom_point(data=map.data, aes(x=Long, y=Lat), colour= "#3FD092", size = 3) +
  labs(x="Longitude", y="Latitude") +
  theme_classic() + labs(title= "Studies testing the effects of nutrient concentration on grazing rates")


## 12.2  Calculate effect size (ES) for nitrogen addition (N)

## 5.1.1 with the log transformed ratio of means (LRR) and its variance (LRR_var) measure= ROM. 
N_ES <-escalc(measure="ROM",
              m1i=Enriched_mean,
              m2i=Control_mean,
              sd1i=Enriched_SD,
              sd2i=Control_SD,
              n1i=N_Treat,
              n2i=N_Contr, 
              data=data_N,var.names=c("EffectSize_N","Variance_N"),digits=4)

hist(N_ES$Effectsize_N)

data_N$EffectSize_N <- N_ES$EffectSize_N
data_N$ESVariance_N <- N_ES$Variance_N
