library(lme4)
library(DHARMa)
library(partR2)
library(car)
library(effects)
library(MuMIn)
library(dplyr)
library(effectsize)
options(na.action='na.fail') # for partR2


## change names in fig

# load data
data <- read.table("data_hb2021.csv", header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
survival <- read.table("wintermortality2021_binom.csv", header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
survival <- survival[-c(101:104),] # delete N1145 (had to be removed from the site per land-owner's request)


# scale variables
data$SNH750 <- data$SNH750*100/(pi*750^2/10000)
data$Org_crop2000 <- data$Org_crop2000*100/(pi*2000^2/10000)
data$Annual_fl500 <- data$Annual_fl500*100/(pi*500^2/10000)
data$OSR500 <- data$OSR500*100/(pi*500^2/10000)
data$Z.p_rich1 <- scale(data$p_rich1)
data$Z.p_rich2 <- scale(data$p_rich2)
data$Z.Annual_fl500 <- scale(data$Annual_fl500)
data$Z.Org_crop2000 <- scale(data$Org_crop2000)
data$Z.SNH750 <- scale(data$SNH750)
data$Z.OSR500 <- scale(data$OSR500)
data$Varroa_p <- data$varroa_board/data$no_bees*100
data$Z.Varroa <- scale(data$Varroa_p)


###################################
### 1. PARASITE  RICHNESS MODEL ###
###################################

# trying to fit parasite richness (p_rich2) to different distributions
#prich.M <- glmer(cbind(p_rich2, 11-p_rich2) ~ Z.Varroa + Z.Org_crop2000 + (1|site), family = 'binomial', data) ### doesn't work well, assumes all parasites have equal probability of being contracted - not true
#prich.M <- glmer(p_rich2 ~ Z.Varroa + Z.Org_crop2000 + (1|site), family = 'poisson', data) ### not good either - underdispersed
#prich.M <- glmer(p_rich2 ~ Z.Varroa + Z.Org_crop2000 + (1|site), family = 'nbinom1', data) ### not good either - underdispersed
#data$p_rich2_r <- data$p_rich2/11
#library(glmmTMB)
#data$p_rich2_P <- data$p_rich2/11
#prich.M <- glmmTMB(p_rich2_r ~ Z.Varroa + Z.Org_crop2000 + (1|site), family = beta_family(), data) ### works well but cannot be used in piecewiseSEM



#using logarithm of parasite richness

prich.M <- glmer(log(p_rich2) ~ Z.p_rich1 + Z.Varroa + (Z.Org_crop2000 + Z.Annual_fl500 +  Z.SNH750 + Z.OSR500)^2 + (1|site), family = Gamma(link='log'), data)

# model selection
b<-dredge(prich.M, m.lim = c(0, 5), rank="AIC") 
b

# best model
prich.M <- lmer(log(p_rich2) ~ Varroa_p + Z.Org_crop2000 + (1|site), REML=T, data) #works ok but not recommended
vif(prich.M)
summary(prich.M)
car::Anova(prich.M, type=3)
plot(allEffects(prich.M))
testDispersion(prich.M)
plot(so <- simulateResiduals(fittedModel = prich.M, plot = F))
r.squaredGLMM(prich.M)


#######################
### 2. VARROA MODEL ###
#######################


#trying to fit Varroa
#data$Varroa_L <- log(data$Varroa_p)
#var.M <- lmer(Varroa_L ~ Z.SNH750 + Z.Annual_fl500 + (1|site), REML=T, data) ### works alright, but recommended only as a last resort
#var.M <- glmer(Varroa_p ~ Z.SNH750 + Z.Annual_fl500 + (1|site), family = Gamma(link='log'), data) #same result as logarithm, but  underdispersed (though not significantly)
#var.M <- glmmTMB(Varroa_p ~ Z.SNH750 + Z.Annual_fl500 + (1|site), family = Gamma(link='log'), REML=T, data) #same result, better dispersion but cannot be used in piecewiseSEM
#var.M <- glmer(varroa_board ~  Z.Annual_fl500 + Z.SNH750 + (1|site), family = 'poisson', weights=no_bees, data) #poisson with weights as no_bees (each bee had an equal prob. of aquiring a varroa mite)


var.M <- glmer(Varroa_p ~ (Z.Org_crop2000 +  Z.SNH750 + Z.OSR500 + Z.Annual_fl500)^2 +  (1|site), family = Gamma(link='log'), data)
vif(var.M)

# model selection
b<-dredge(var.M, m.lim = c(0, 4), rank="AIC") 
b

# best model with interaction, but VIF is inflated - removed interaction
var.M <- glmer(Varroa_p ~ Z.SNH750 + Z.Annual_fl500 + (1|site), family = Gamma(link='log'), data)

vif(var.M)
summary(var.M)
car::Anova(var.M, type=3)
plot(allEffects(var.M))
testDispersion(var.M)
plot(so <- simulateResiduals(fittedModel = var.M, plot=T))
r.squaredGLMM(var.M)


##############################
### 3. COLONY GROWTH MODEL ###
##############################



col.M <- glmer.nb(no_bees ~ log(p_rich2) + Z.Varroa + (Z.OSR500 + Z.Org_crop2000 + Z.SNH750 + Z.Annual_fl500)^2 + (1|site), family='nbinom2', data) 

# model selection
b<-dredge(col.M, m.lim = c(0, 3), rank="AIC") 
b

# best model 
col.M <- glmer.nb(no_bees ~ Z.Org_crop2000 + SNH750 + log(p_rich2)  + (1|site), family='nbinom2',data) # a bit underdispersed
#col.M <- lmer(no_bees ~ Z.Org_crop2000 + SNH750 + L.p_rich2 + (1|site), data)  - normal distribution less underdispersed but the same result

vif(col.M)
summary(col.M)

car::Anova(col.M, type=3)
plot(allEffects(col.M))
testDispersion(col.M)
plot(so <- simulateResiduals(fittedModel = col.M, plot = F))
r.squaredGLMM(col.M)


##### individual parasites instead of parasite richness, model 3b:

col.M <- glmer.nb(no_bees ~ dwvb + tryp + bqcv + Z.Varroa + (Z.OSR500 + Z.Org_crop2000 + Z.SNH750 + Z.Annual_fl500)^2 + (1|site), family='nbinom2', data) 


b<-dredge(col.M, m.lim = c(0, 3), rank="AIC") 
b

col.M <- glmer.nb(no_bees ~ tryp  + Z.Org_crop2000 + Z.SNH750 +  (1|site), control = glmerControl(optimizer = 'bobyqa'), family='nbinom2',data)
vif(col.M)
summary(col.M)

car::Anova(col.M, type=3)
plot(allEffects(col.M))
testDispersion(col.M)
plot(so <- simulateResiduals(fittedModel = col.M, plot = F))
r.squaredGLMM(col.M)


#################################
### 4. COLONY SURVIVAL MODEL  ###
#################################


data_S <- merge(data, survival, by='sample') #merging with survival - 2 colonies less, because they had to been removed from the site in the fall
data_S$Z.no_bees <- scale(data_S$no_bees)

# not going through as one model - splitting into two
surv.M <- glmer(status ~ Z.Org_crop2000 + Z.Annual_fl500 + Z.SNH750 * Z.OSR500 + (1|site), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial, data_S)
surv.M <- glmer(status ~ log(p_rich2) + Z.Varroa + Z.No_bees (1|site), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), family=binomial, data_S)

# model selection - first the landscape one, then the parasite one
b<-dredge(surv.M, m.lim = c(0, 3), rank="AICc") 
b

# best landscape model = NULL
# best parasite model
surv.M <- glmer(status ~  Z.no_bees + (1|site), family=binomial, data_S)

summary(surv.M)
car::Anova(surv.M, type=3)
plot(allEffects(surv.M))
testDispersion(surv.M)
plot(so <- simulateResiduals(fittedModel = surv.M, plot = F))
r.squaredGLMM(surv.M)



###############
### 5. SEM  ###
###############

library(piecewiseSEM)
data$L.p_rich2 <- log(data$p_rich2)

# best models - skipped colony survival
m1 <- glmer(Varroa_p ~ Z.Annual_fl500 + Z.SNH750 + (1|site), family = Gamma(link='log'), data) 
m2 <- lmer(L.p_rich2 ~ Varroa_p + Z.Org_crop2000 + (1|site), REML=T, data) 
m3 <- glmer.nb(no_bees ~ L.p_rich2 + Z.SNH750 + Z.Org_crop2000 + (1|site),  family = 'nbinom2', data) 

sem <- psem(m1, m2, m3)

summary(sem)

plot(sem) #not working since R update


########### Standarization of the estimates for Gamma distribution 
########### (not implemented in piecewiseSEM) via observation error approach

summary(m1)
beta_annfl <- summary(m1)$coefficients[2, 1]
beta_SNH <- summary(m1)$coefficients[3, 1]
# Extract predicted values on the link scale
preds <- predict(m1, type = "link")
r.squaredGLMM(m1)
# Compute sd of error variance using theoretical variances
sd.y <- sqrt(var(preds) + 0.52)
# Compute sd of x (which is 1 because the variables are standardized)
sd.x.annfl <- sd(data$Z.Annual_fl500)
sd.x.SNH <- sd(data$Z.SNH750)

# std. beta of annual flower
beta_annfl * sd.x / sd.y

#std. beta of SNH
beta_SNH * sd.x / sd.y



#############################################
### Recalculating coefficients from GLMMs ###
#############################################


## for Gamma and negative binomial (log link) - exponentiating:
#exp(X)

## colony growth model
## Org_crop = 1.17
## SNH = 0.9
## p_rich = 0.82

## for binomial models (logit link) - reversing log-odds:
#exp(X)/(1+exp(X))


#############################################
### INDIVIDUAL PARASITE PREVALENCE MODELS ###
#############################################


# loading data from first sampling
data1 <- read.table("data_hb2021_1sampling.csv", header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
data12 <- merge(data, data1, by='sample')

# DWV-A <10% positives

# DWV-B
dwvb.M <- glmer(dwvb ~ dwvb1 + Z.Varroa + (Z.Org_crop2000 + Z.Annual_fl500 + Z.SNH750 + Z.OSR500)^2 + (1 | site), family='binomial', data=data12) #excluding 2-way interactions becasue of VIF

b<-dredge(dwvb.M, m.lim = c(0, 4), rank="AIC") 
b

dwvb.M <- glmer(dwvb ~ dwvb1 + Z.OSR500 * Z.Org_crop2000 + (1 | site), family='binomial', data=data12) #third best model because of heterogeneity
vif(dwvb.M)
summary(dwvb.M)
plot(allEffects(dwvb.M))
testDispersion(dwvb.M)
plot(so <- simulateResiduals(fittedModel = dwvb.M, plot = F))

#BQCV - no interaction term because of the error
bqcv.M <- glmer(bqcv ~ bqcv1 + Z.Org_crop2000  + Z.OSR500 + Z.Annual_fl500 + Z.SNH750 + (1 | site), family='binomial', control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12), data=data12)

b<-dredge(bqcv.M, m.lim = c(0, 4), rank="AIC") 
b

bqcv.M <- glmer(bqcv ~ Z.Org_crop2000  + Z.OSR500 + Z.SNH750 + (1 | site), family='binomial', data=data12) #third best model because of convergence problems
summary(bqcv.M)
plot(allEffects(bqcv.M))
testDispersion(bqcv.M)
plot(so <- simulateResiduals(fittedModel = bqcv.M, plot = F))


# SBV <10% positives
# ABPV <10% positives
# CBPV <10% positives
# SBPV <10% positives

# Trypanosomes
tryp.M <- glmer(tryp ~ tryp1 + Z.Org_crop2000 + Z.Annual_fl500 + Z.SNH750 + Z.OSR500 + (1 | site), family='binomial', control=glmerControl(nAGQ0initStep=F, tolPwrss=1e-12), data=data12)

b<-dredge(tryp.M, m.lim = c(0, 5), rank="AIC") 
b

tryp.M <- glmer(tryp ~ tryp1  + (1 | site), family='binomial', data=data12)
summary(tryp.M)
plot(allEffects(tryp.M))
testDispersion(tryp.M)
plot(so <- simulateResiduals(fittedModel = tryp.M, plot = F))


# Neogregarines <10% positives
# Nosema bombi <10% positives
# Nosema ceranae <10% positives
  

############################################
### MORAN'S I - SPATIAL AUTOCORRELATION  ###
############################################


library(spdep)
coord <- read.table("C:/Users/patry/OneDrive/PhD/tralala/combee/RAW_Data_2021/prev_map.csv", 
                    header=TRUE, stringsAsFactors=TRUE, sep=",", na.strings="NA", dec=".", 
                    strip.white=TRUE)
coord <- coord %>% distinct()

# parasite  richness
plot(so <- simulateResiduals(fittedModel = m2, plot = F))
so2 <- recalculateResiduals(so, group = data$site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# colony growth
plot(so <- simulateResiduals(fittedModel = m3, plot = F))
so2 <- recalculateResiduals(so, group = data$site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# varroa infestation
plot(so <- simulateResiduals(fittedModel = m1, plot = F))
so2 <- recalculateResiduals(so, group = data$site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)
# colony survival
plot(so <- simulateResiduals(fittedModel = surv.M, plot = F))
so2 <- recalculateResiduals(so, group = data_S$site) #because of 2 data points from each location
testSpatialAutocorrelation(so2, coord$long, coord$lat, plot = T)