library(multcomp)
library(dplyr)
library(glmmTMB)
library(lme4)
library(car)

standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x))


############## WITHOUT WEIRD REPS ############## 

data <- read.csv('EXP4-table.csv', header=TRUE, sep=',')
data <- read.csv('no_stressors.csv', header=TRUE, sep=',')

data <- data %>%
  mutate(day = case_when(
    replicate %in% c(1, 2) ~ 1,
    replicate %in% c(3, 4) ~ 2,
    replicate %in% c(5, 6) ~ 3,
    replicate %in% c(7, 8) ~ 4,
    replicate %in% c(9, 10) ~ 5,
    replicate %in% c(11, 12) ~ 6,
    replicate %in% c(13, 14) ~ 7
  ))

means <- aggregate( nb_shakes ~ treatment, FUN=mean, data=data)
#sd <- aggregate( nb_shakes ~ treatment, FUN=sd, data=data)
error <- aggregate( nb_shakes ~ treatment, FUN=standard_error, data=data)

means$se <- error$nb_shakes
write.csv(means, 'no_stressors_means.csv', row.names = FALSE)


data$nb_shakes <- as.numeric(data$nb_shakes)
data$Observation.id <- as.factor(data$Observation.id)
data$replicate <- as.factor(data$replicate)
data$treatment <- as.factor(data$treatment)
data$day <- as.factor(data$day)



#### myrmica vs nestmate vs non nestmate ###

data_myrmica_nestmate <- subset(data, data$treatment=='nestmate' | data$treatment=='non-nestmate' | data$treatment=='myrmica')

data_mean <- aggregate(nb_shakes ~ treatment, FUN=mean, data=data_myrmica_nestmate)
data_sd <- aggregate(nb_shakes ~ treatment, FUN=sd, data=data_myrmica_nestmate)
data_se <- aggregate(nb_shakes ~ treatment, FUN=standard_error, data=data_myrmica_nestmate)

data_mean$data_sd <- data_sd$nb_shakes
data_mean$data_se <- data_se$nb_shakes


model_myrmica_nestmate <- glmer(nb_shakes ~ treatment + (1|replicate) + (1|day), family='poisson', data=data_myrmica_nestmate)
Anova(model_myrmica_nestmate, type='II')
post_hoc <- summary(glht(model_myrmica_nestmate, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
aggregate(nb_shakes ~ treatment, FUN=mean, data=data_myrmica_nestmate)
print(cld(post_hoc))



#### nestmate cadaver vs sporulating cadaver ###

data_cadavers <- subset(data, data$treatment=='nestmate-cad' | data$treatment=='sporulating-cad')
model_cadavers <- glmer(nb_shakes ~ treatment + (1|replicate) + (1|day), family='poisson', data=data_cadavers)
Anova(model_cadavers, type='II')
post_hoc <- summary(glht(model_cadavers, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
aggregate(nb_shakes ~ treatment, FUN=mean, data=data_cadavers)
print(cld(post_hoc))

data_mean <- aggregate(nb_shakes ~ treatment, FUN=mean, data=data_cadavers)
data_sd <- aggregate(nb_shakes ~ treatment, FUN=sd, data=data_cadavers)
data_se <- aggregate(nb_shakes ~ treatment, FUN=standard_error, data=data_cadavers)

data_mean$data_sd <- data_sd$nb_shakes
data_mean$data_se <- data_se$nb_shakes


#### food vs brood vs glass ball ###

data_positives <- subset(data,  data$treatment=='brood' | data$treatment=='food' | data$treatment=='glass-ball')
model_positives <- glmer(nb_shakes ~ treatment + (1|replicate) + (1|day), family='poisson', data=data_positives)
Anova(model_positives, type='II')
post_hoc <- summary(glht(model_positives, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
aggregate(nb_shakes ~ treatment, FUN=mean, data=data_positives)
print(cld(post_hoc))


data_mean <- aggregate(nb_shakes ~ treatment, FUN=mean, data=data_positives)
data_sd <- aggregate(nb_shakes ~ treatment, FUN=sd, data=data_positives)
data_se <- aggregate(nb_shakes ~ treatment, FUN=standard_error, data=data_positives)

data_mean$data_sd <- data_sd$nb_shakes
data_mean$data_se <- data_se$nb_shakes



#### first shakes ###


first_ant <- read.csv('first_ants.csv', header=TRUE, sep=',')

first_ant

data_mean <- aggregate(nb_shakes ~ treatment, FUN=mean, data=first_ant)
data_sd <- aggregate(nb_shakes ~ treatment, FUN=sd, data=first_ant)
data_se <- aggregate(nb_shakes ~ treatment, FUN=standard_error, data=first_ant)

data_mean$data_sd <- data_sd$nb_shakes
data_mean$data_se <- data_se$nb_shakes

write.csv(data_mean, 'first_ant_means_sds.csv', row.names=FALSE)

aggregate(nb_shakes ~ treatment, FUN=mean, data=first_ant)

aggregate(nb_shakes ~ treatment, FUN=standard_error, data=first_ant)

model_first <- glmer(nb_shakes ~ treatment + (1|replicate), family='poisson', data=first_ant)
Anova(model_first, type='II')

post_hoc <- summary(glht(model_first, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))


## pathogen ##

pathogen <- subset(first_ant, treatment=='sporulating-cad' | treatment == 'nestmate-cad')

pathogen_model <- glmer(nb_shakes ~ treatment + (1|replicate), family='poisson', data=pathogen)
Anova(pathogen_model, type='II')
post_hoc <- summary(glht(pathogen_model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
print(cld(post_hoc))

## intruder ##

intruder <- subset(first_ant, treatment=='nestmate' | treatment == 'non-nestmate' | treatment == 'myrmica')

means <- aggregate(nb_shakes ~ treatment, data = intruder, mean)

intruder$treatment <- factor(intruder$treatment, levels = means[order(means$nb_shakes), "treatment"])

intruder_model <- glmer(nb_shakes ~ treatment + (1|replicate), family='poisson', data=intruder)
Anova(intruder_model, type='II')
post_hoc <- summary(glht(intruder_model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
print(cld(post_hoc))

## positive ##

positive <- subset(first_ant, treatment=='brood' | treatment == 'food' | treatment == 'glass-ball')

positive_model <- glmer(nb_shakes ~ treatment + (1|replicate), family='poisson', data=positive)
Anova(positive_model, type='II')
post_hoc <- summary(glht(positive_model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
print(cld(post_hoc))




## shakes directed towards stressors

data <- read.csv('shakes_stressors_vs_ants.csv', header=TRUE, sep = ',')
data


data_means <- aggregate(nb_shakes ~ treatment + direction, FUN=mean, data=data)
data_sd <- aggregate(nb_shakes ~ treatment + direction, FUN=sd, data=data)
data_se <- aggregate(nb_shakes ~ treatment + direction, FUN=standard_error, data=data)

data_means$sd <- data_sd$nb_shakes
data_means$se <- data_se$nb_shakes


write.csv(data_means, 'stressor_means4plots.csv')


myrmica <- subset(data, treatment=='myrmica')
aggregate(nb_shakes ~ direction, FUN=mean, data=myrmica)

model <- glmer(nb_shakes ~ direction + (1|replicate), family='poisson', data=myrmica)
Anova(model, type='II')

post_hoc <- summary(glht(model, linfct=mcp(direction="Tukey")), test=adjusted("BH"))
print(cld(post_hoc))


stressor <- subset(data, direction=='stressor')
stressor

stressor <- stressor %>%
  mutate(day = case_when(
    replicate %in% c(1, 2) ~ 1,
    replicate %in% c(3, 4) ~ 2,
    replicate %in% c(5, 6) ~ 3,
    replicate %in% c(7, 8) ~ 4,
    replicate %in% c(9, 10) ~ 5,
    replicate %in% c(11, 12) ~ 6,
    replicate %in% c(13, 14) ~ 7
  ))


data_sums <- aggregate(nb_shakes ~ treatment, FUN=sum, data=stressor)
data_means <- aggregate(nb_shakes ~ treatment, FUN=mean, data=stressor)
data_sd <- aggregate(nb_shakes ~ treatment, FUN=sd, data=stressor)
data_se <- aggregate(nb_shakes ~ treatment, FUN=standard_error, data=stressor)

data_means$sum <- data_sums$nb_shakes
data_means$sd <- data_sd$nb_shakes
data_means$se <- data_se$nb_shakes




model_nb <- glmer.nb(nb_shakes ~ treatment + (1|replicate) + (1|day), data = stressor)
Anova(model_nb, type='II')
post_hoc <- summary(glht(model_nb, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))




#### INCREASE FOLD COMPARISONS ####


data <- read.csv('increase_fold.csv', header=TRUE, sep=',')

data <- data %>%
  mutate(day = case_when(
    replicate %in% c(1, 2) ~ 1,
    replicate %in% c(3, 4) ~ 2,
    replicate %in% c(5, 6) ~ 3,
    replicate %in% c(7, 8) ~ 4,
    replicate %in% c(9, 10) ~ 5,
    replicate %in% c(11, 12) ~ 6,
    replicate %in% c(13, 14) ~ 7
  ))


data$adjusted_nb_shakes <- as.numeric(data$adjusted_nb_shakes)
data$Observation.id <- as.factor(data$Observation.id)
data$replicate <- as.factor(data$replicate)
data$treatment <- as.factor(data$treatment)
data$day <- as.factor(data$day)


model <- lmer(sqrt(adjusted_nb_shakes) ~ treatment + (1|replicate) + (1|day), data=data)

shapiro.test(residuals(model))
hist(residuals(model))
Anova(model, type='II')
summary(glht(model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))

aggregate(adjusted_nb_shakes ~ treatment, FUN=mean, data=data)
