###############
### FIGURES ###
###############
library(ggeffects)
library(ggplot2)

###############################
# 1. parasite  richness model #
###############################

#unscaled model
prich.M <- lmer(log(p_rich2) ~ Varroa_p + Org_crop2000 + (1|site), data)

# effect plots
plot(allEffects(prich.M, residuals=TRUE),
     partial.residuals=list(span=0.001), selection=1)

plot(allEffects(prich.M, residuals=TRUE),
     partial.residuals=list(span=0.001), selection=2)

# predictions
plot(ggpredict(prich.M, terms = 'Varroa_p', back.transform = T)) + 
  theme_classic(base_size = 24)+
  theme(axis.text = element_text(color='black'))+
  ggtitle(NULL)+
  xlab('Daily mite fall per 100 bees') +
  ylim(0,7.8)+
  ylab('Predicted parasite richness')+
  ggtitle('B')
ggsave('par_var.png', device = 'png', height = 5.5, width=6)

  
plot(ggpredict(prich.M, terms = 'Org_crop2000', back.transform = T)) + 
  theme_classic(base_size = 24)+
  theme(axis.text = element_text(color='black'))+
  ggtitle(NULL)+
  xlab('Organic farming [%]') +
  ylim(0,7.8)+
  ylab('Predicted parasite richness')+
  ggtitle('A')
ggsave('par_org.png', device = 'png', height = 5.5, width=6)


#more plot editing
#p <- ggpredict(var.M, terms = 'SNH750[0:21 by=0.1]')
#p <- ggpredict(col.M, terms = 'Org_crop2000[0:11 by=0.1]')

#ggplot(p, aes(x, predicted))+
#  geom_smooth(color = 'black', se = F, size = 0.8)+
#  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = 'grey', alpha = 0.4)+
#  theme_classic(base_size = 16)+
#  xlab('Daily varroa mite fall per 100 bees') +
#  ylab('Predicted parasite richness')+
#  ggtitle('B')
#ggsave('par_var.png', device = 'png', height = 5.5, width=6)

###################
# 2. Varroa model #
###################

var.M <- glmer(Varroa_p ~ SNH750 + Annual_fl500 + (1|site), family = Gamma(link='log'), data)

plot(ggpredict(var.M, terms = 'SNH750[0:21 by=0.1]')) + 
  theme_classic(base_size = 24)+
  theme(axis.text = element_text(color='black'))+
  ggtitle(NULL)+
  xlab('Perennial SNH [%]') +
  ylab('Predicted daily mite fall')+
  ggtitle('B')
ggsave('var_snh.png', device = 'png', height = 5.5, width=6)


plot(ggpredict(var.M, terms = 'Annual_fl500[0:2.6 by=0.1]')) + 
  theme_classic(base_size = 24)+
  theme(axis.text = element_text(color='black'))+
  ggtitle(NULL)+
  xlab('Annual flower area [%]') +
  ylab('Predicted daily mite fall')+
  ggtitle('A')
ggsave('var_anf.png', device = 'png', height = 5.5, width=6)


##########################
# 3. colony growth model #
##########################

col.M <- glmer.nb(no_bees ~ Org_crop2000 + SNH750 + log(p_rich2)  + (1|site), family='nbinom2',data)

plot(ggpredict(col.M, terms = 'Org_crop2000[0:10.4 by=1]'))+ 
  theme_classic(base_size = 24)+
  ggtitle(NULL)+
  xlab('Organic farming [%]') +
  ylab('Predicted number of bees')+
  ggtitle('A')+
  ylim(4000,16000)
ggsave('col_org.png', device = 'png', height = 5.5, width=6)

plot(ggpredict(col.M, terms = 'p_rich2[1:5 by=0.2]', back.transform = T))+ 
  theme_classic(base_size = 24)+
  ggtitle(NULL)+
  xlab('Parasite richness') +
  ylab('Predicted number of bees')+
  ggtitle('B')+
  ylim(4000,16000)
ggsave('col_prich.png', device = 'png', height = 5.5, width=6)

plot(ggpredict(col.M, terms = 'SNH750'))+ 
  theme_classic(base_size = 16)+
  ggtitle(NULL)+
  xlab('Perennial SNH [%]') +
  ylab('Predicted number of bees')+
  ggtitle('')
ggsave('col_snh.pdf', device = 'pdf', height = 6, width=6)



############################
# 3. colony survival model #
############################

surv.M <- glmer(status ~  no_bees + (1|site), family=binomial, data_S)

plot(ggpredict(surv.M, terms = 'no_bees'))
