# analysis of association between rodenticide poisoning and zoonotic infection 

# run calculateSMI.R script before running this one

# Ben Bolker's function for calculating R2
# based on Nakagawa & Schielzeth 2013
# https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/inst/misc/rsqglmm.R
source("./R/rsqglmm.R")

# select columns for analysis---------------------------------------------------
clean <- data %>% 
  select(poisoned, lepto, ecoli, sex, ageClass, communityArea, SMI) %>% 
  mutate(across(lepto:communityArea, as.factor)) 

write.csv(clean, "./Data/cleanData.csv", row.names = FALSE)

# summary of tested rats for Table 1--------------------------------------------
xtabs(~communityArea + sex, data = clean)
xtabs(~communityArea + ageClass, data = clean)
xtabs(~communityArea + poisoned, data = clean)

# look more closely at poisoned rats
poisoned <- filter(clean, poisoned == 1)
summary(poisoned)

# lepto model-------------------------------------------------------------------

m.lepto <- glmmTMB(lepto ~ poisoned + sex + ageClass + SMI + (1|communityArea), 
                   data = clean, family = binomial(link = "logit"))

summary(m.lepto)
round(exp(confint(m.lepto)), 2)

# calculate R2
my_rsq(m.lepto)
# marginal R2 describes proportion of variance explained by only fixed effects
# conditional R2 describes proportion of variance explained by fixed and random

# ecoli model-------------------------------------------------------------------

m.ecoli <- glmmTMB(ecoli ~ poisoned + sex + ageClass + SMI + (1|communityArea), 
                    data = clean, family = binomial(link = "logit"))

summary(m.ecoli)
round(exp(confint(m.ecoli)), 2)

# calculate R2
my_rsq(m.ecoli)
