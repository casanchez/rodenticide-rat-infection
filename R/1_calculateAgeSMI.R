# code for calculating age and SMI

source("./R/loadPackages.R")

# read in data
data <- read.csv("./Data/rodenticideInfectionData.csv", 
                 na.strings = c("N/A", ""), fileEncoding = 'UTF-8-BOM') 

# estimate age based on weight--------------------------------------------------

# impute rat age from weight using von Bertalanffy equation
# weight = a * (1-exp(-r*(age - c)))
# can be re-written as age = c - log(1 - weight/age)/r

# using parameters from Minter et al. 2017: 
# "Evidence of multiple intraspecific transmission routes 
# for Leptospira acquisition in Norway rats (Rattus norvegicus)
c <- 23 # age at which maximum growth occurs (days)
a <- 562 # asymptote (days)
r <- 0.01337 # constant growth rate (grams per day)

data %<>% 
  # rename some columns for easier use
  rename(communityArea = Community.Area, mass = Mass..g., 
         tipTip = Tip.to.tip.length..cm., tailLength = Tail.length..cm.,
         bodyLength = Body.length..cm., earLength = Ear.length..mm., 
         footLength = Hind.foot.length..mm., sex = Sex, lepto = Leptospira, 
         ecoli = Ecoli, poisoned = AR_binary) %>% 
  # calculate age in days, rounding up to nearest day
  mutate(ageDays = ceiling(c - log(1 - mass/a)/r)) 

# look at age distribution to decide breakpoints
range(data$ageDays)
hist(data$ageDays, breaks = seq(30, 220, 5))
table(data$ageDays)

# seems like there's a natural break around 65

data %<>% 
  mutate(ageClass = case_when(
    ageDays %in% 30:65 ~ "30-65",
    ageDays >65 ~ ">65")) %>% 
  mutate(ageClass = fct_relevel(ageClass, "30-65"))

# SMI step 1: bivariate plots of mass and length observations-------------------

# Peig & Green "recommend the use of that single L variable which has 
# the strongest correlation with M on a log-log scale"

par(mfrow = c(2, 3))

plot(log(mass) ~ log(tipTip), data = data)
lm1 <- lm(log(mass) ~ log(tipTip), data = data)
abline(lm1)

plot(log(mass) ~ log(tailLength), data = data)
lm2 <- lm(log(mass) ~ log(tailLength), data = data)
abline(lm2)

plot(log(mass) ~ log(bodyLength), data = data)
lm3 <- lm(log(mass) ~ log(bodyLength), data = data)
abline(lm3)

plot(log(mass) ~ log(earLength), data = data)
lm4 <- lm(log(mass) ~ log(earLength), data = data)
abline(lm4)

plot(log(mass) ~ log(footLength), data = data)
lm5 <- lm(log(mass) ~ log(footLength), data = data)
abline(lm5)

summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)
summary(lm5)

# based on R2, tipTip is most correlated with mass (0.76)
# followed by tailLength (0.71), bodyLength (0.68), footLength(0.36), and
# earLength(0.21). footLength and earLength appear most influenced by outliers

# we can remove outliers to further refine estimates
# interactive plots are useful to ID outliers
p1 <- ggplot(data = data) +
  geom_point_interactive(aes(x = log(tipTip), y = log(mass), tooltip = Rat.ID)) +
  theme_minimal()
girafe(ggobj = p1)

noOut1 <- data %>% 
  filter(!(Rat.ID %in% c("BE01", "LV27")))


p2 <- ggplot(data = data) +
  geom_point_interactive(aes(x = log(tailLength), y = log(mass), 
                             tooltip = Rat.ID)) +
  theme_minimal()
girafe(ggobj = p2)

noOut2 <- data %>% 
  filter(!(Rat.ID %in% c("LV27")))


p3 <- ggplot(data = data) +
  geom_point_interactive(aes(x = log(bodyLength), y = log(mass), 
                             tooltip = Rat.ID)) +
  theme_minimal()
girafe(ggobj = p3)

noOut3 <- data %>% 
  filter(!(Rat.ID %in% c("BE01", "LS03")))


p4 <- ggplot(data = data) +
  geom_point_interactive(aes(x = log(earLength), y = log(mass), 
                             tooltip = Rat.ID)) +
  theme_minimal()
girafe(ggobj = p4)

noOut4 <- data %>% 
  filter(!(Rat.ID %in% c("SL15")))


p5 <- ggplot(data = data) +
  geom_point_interactive(aes(x = log(footLength), y = log(mass), 
                             tooltip = Rat.ID)) +
  theme_minimal()
girafe(ggobj = p5)

noOut5 <- data %>% 
  filter(!(Rat.ID %in% c("BE01", "NN14")))

# replot now that outliers have been removed
par(mfrow = c(2, 3))
plot(log(mass) ~ log(tipTip), data = noOut1)
lm1 <- lm(log(mass) ~ log(tipTip), data = noOut1)
abline(lm1)

plot(log(mass) ~ log(tailLength), data = noOut2)
lm2 <- lm(log(mass) ~ log(tailLength), data = noOut2)
abline(lm2)

plot(log(mass) ~ log(bodyLength), data = noOut3)
lm3 <- lm(log(mass) ~ log(bodyLength), data = noOut3)
abline(lm3)

plot(log(mass) ~ log(earLength), data = noOut4)
lm4 <- lm(log(mass) ~ log(earLength), data = noOut4)
abline(lm4)

plot(log(mass) ~ log(footLength), data = noOut5)
lm5 <- lm(log(mass) ~ log(footLength), data = noOut5)
abline(lm5)

summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)
summary(lm5)

# tipTip is still most correlated with mass
# R2 has also been slightly improved (0.788)

# last step, check to see if strength of fit improves when only
# rats from the reference population (ie not exposed to ARs) are included 
plot(log(mass) ~ log(tipTip), data = subset(noOut1, poisoned == 0))
lm1 <- lm(log(mass) ~ log(tipTip), data = subset(noOut1, poisoned == 0))
abline(lm1)

summary(lm1)
# R2 has been slightly improved again (0.792)

# one possible explanation is that animals exposed to ARs have inhibited
# or disrupted growth (disrupting bSMA)

# step 2: fitting an SMA regression line to ln-transformed data-----------------

# to get best estimate of slope, 
# we will use data from non-AR exposed rats with outliers excluded

reliableData <- subset(noOut1, poisoned == 0)

# first, let's check if the scaling relationship between M and L is the same
# between sexes
sma(log(mass) ~ log(tipTip)*sex, data = reliableData)
# P = 0.051, so we fail to reject null H0 that the slopes are equal

# and how about between age classes
sma(log(mass) ~ log(tipTip)*ageClass, data = reliableData)
# P = 0.469 so we fail to reject H0 again 

# this means we can use the same bSMA for all rats
# and we can make valid comparisons between sexes and age classes

SMAres <- sma(log(mass) ~ log(tipTip), data = reliableData)

bSMA <-  SMAres$coef[[1]][2, 1]
# 3.05 aligns with expected value of ~3 under assumption of isometric growth

# step 3: select L0 and compute SMI---------------------------------------------

# Peig and Green: "The scaled mass index computes the mass each individual 
# would have at a fixed body size as indicated by L0. 
# The constant L0 denotes a particular point along the M-L relationship, 
# i.e. a specific body size along the power curve. 
# The arithmetic mean or median of the linear body measure (L) 
# for all sampled individuals make suitable L0 values, 
# but we emphasize that any value within the range of L observations can be used.

# let's use mean of tipTip for all individuals (including outliers)
L0 <- mean(data$tipTip, na.rm = T)

# now can calculate SMI for all individuals (including outliers)
data %<>% 
  mutate(SMI = mass*(L0 / tipTip)^bSMA)


# according to Peig & Green:
# "Unless some factor affects the well-being of animals at a specific growth
# stage...mean condition scores should be equal for different age classes...
# Similarly, even if a species shows sexual size dimorphism...we should expect 
# no difference in CIS (i.e. well-being) between sex.
younger <- filter(data, ageClass == "30-65")
older <- filter(data, ageClass == ">65")

mean(younger$SMI, na.rm = T) # 167.14
mean(older$SMI, na.rm = T) # 162.16

m <- filter(data, sex == "Male")
f <- filter(data, sex == "Female")

mean(m$SMI, na.rm = T) # 164.86
mean(f$SMI, na.rm = T) # 166.83

# indeed, average SMIs are similar between groups 
