
rm(list=ls())


setwd("C:/home/rb17990/Documents/EXP1_code")
events <- read.csv(file="events.csv", header=T, sep=",")


# install.packages(
#   "ggplot2",
#   repos = c("http://rstudio.org/_packages",
#             "http://cran.rstudio.com")
# )
# install.packages("tidyverse")
#  install.packages("pscl")
#  install.packages("performance")
#  install.packages("ggbeeswarm")
#  install.packages("R2admb")
#  install.packages("glmmADMB",
#                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                         getOption("repos")),
#                 type="source")
# install.packages("glmmTMB")
# install.packages("blmeco")
# install.packages("ggplot2")
# install.packages("car")
# install.packages("multcomp")
# install.packages("MuMIn")
# install.packages("dplyr")
# install.packages("e1071")
# install.packages("plyr")

library(blmeco)
library(performance)
library(stringr)
library(lme4)
library(car)
library(pscl)
library(glmmADMB)
library(ggplot2)
library(multcomp)
library(glmmTMB)
library(ggbeeswarm)
library(numDeriv)
library(MuMIn)
library(dplyr)
library(e1071)
library(plyr)

#allows expand grid function to pass a data frame
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

ant_list <- paste("ANT",c(1:11),sep="_")

##removing spaces in events
##first part removes leading spaces 
events$Replicate <- trimws(events$Replicate, which=c("left"))
events$Duration_.s. <- trimws(events$Duration_.s., which=c("left"))
events$Start_.s. <- trimws(events$Start_.s., which=c("left"))
events$Stop_.s.<- trimws(events$Stop_.s., which=c("left"))

#adds in _ in spaces
events <- as.data.frame(apply(events, 2, function(x)
  str_replace_all(string=x, pattern=" ", repl="_")))


##adding unique identifiers for each trial
events$group_size <- 1+as.numeric(sapply(strsplit(events$Trial, ",_"), getElement, 1))
events$treatment <- sapply(strsplit(events$Trial, ",_"), getElement, 2)
events$ant_ID <- events$Subject

no_events <- events[ which(   events$ant_ID  == "NO_EVENTS"      )  ,   ]
discards  <- events[ which(   events$ant_ID  == "DISCARD"      )  ,   ]
events    <- events[ which( !  events$ant_ID   %in% c("DISCARD","NO_EVENTS")      )  ,   ]

###get non-discarded replicate list
replicate_list <- rbind(unique(events[c("Replicate","treatment","group_size")]),unique(no_events[c("Replicate","treatment","group_size")]))

#ant status
events$ant_status <- "NESTMATE"
events[ which (events$ant_ID=="ANT_1")  , "ant_status"    ] <- "FOCAL"

unique(events$ant_status)

events$time_of_year <- "WINTER"
events[which(grepl("\\/08\\/",events$Observation_date)),"time_of_year"] <- "SUMMER"


events1 <- events[ which(events$Behavior == "T")   ,   ]

write.csv(events1, 'events1.csv', row.names=FALSE)

##adding twitch modifiers into separate columns
#twitch size, repetition, proximity to brood and nestmates
events1$twitch_size <- sapply(strsplit(events1$Modifiers, ","), getElement, 1)
events1$twitch_size <- sapply(strsplit(events1$twitch_size, "_"), getElement, 1)
events1 [ which(events1$twitch_size == "Small") , "twitch_size" ] <- "small"
unique(events1$twitch_size)

events1$twitch_repetition <- sapply(strsplit(events1$Modifiers, ","), getElement, 2)

events1[ which(events1$twitch_repetition == "single_twitch")  , "twitch_repetition" ] <- "single"
events1[ which(events1$twitch_repetition == "part_of_a_repetitive_set")  , "twitch_repetition" ] <- "repetitive"
unique(events1$twitch_repetition)

events1$proximity_to_nestmates <- sapply(strsplit(events1$Modifiers, ","), getElement, 3)
events1[ which(events1$proximity_to_nestmates == "isolated_twitch/alone") , "proximity_to_nestmates" ] <- "alone"
events1[ which(events1$proximity_to_nestmates == "in_contact_with_nestmate") , "proximity_to_nestmates" ] <- "contact"
events1[ which(events1$proximity_to_nestmates == "in_contact_with_multiple_nestmates") , "proximity_to_nestmates" ] <- "multiple_contact"
unique(events1$proximity_to_nestmates)

events1$proximity_to_brood <- sapply(strsplit(events1$Modifiers, ","), getElement, 4 )
unique(events1$proximity_to_brood)

##adding in temporal effects
##if before half way point (600) then time == B, and then after == A
# events1$Start_.s. <- as.numeric(events1$Start_.s.)
# events1$time <- "B"
# events1 [ which ( events1$Start_.s. > 600 ) , "time" ] <- "A"



#each ant gets a bloc of all possible combinations of twitch modifiers
events1$N <- 1




ant_bloc <- expand.grid (  
  twitch_size = unique(events1$twitch_size)
  ,
  twitch_repetition = unique (events1$twitch_repetition)
  ,
  proximity_to_nestmates = unique (events1$proximity_to_nestmates)
  ,
  proximity_to_brood = unique (events1$proximity_to_brood)
  ,
  duration = 20*60
  , 
  time = unique ( events1$time )
  
)

twitch_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list[ 1:   replicate_list[repl,"group_size"] ])
  repl_twitches <- expand.grid.df(   data.frame( replicate_list[repl,] )
                                     ,
                                     ants
                                     ,
                                     ant_bloc
  )
  twitch_table <- rbind(twitch_table , repl_twitches)
}

events_Nb <- events1[c(names(twitch_table)[which(names(twitch_table)%in%names(events1))],"ant_status","Behavior", "N")]
events_Nb <- aggregate(N ~ . , FUN=sum, data= events_Nb)

twitch_table <- merge(twitch_table,events_Nb,all.x=T)

###fill in missing cases
twitch_table[ which(twitch_table$ant_ID=="ANT_1")   ,    "ant_status"] <- "FOCAL"
twitch_table[ which(twitch_table$ant_ID!="ANT_1")   ,    "ant_status"] <- "NESTMATE"

twitch_table$Behavior <- unique(events1$Behavior)

twitch_table  [ is.na(twitch_table$N)  ,"N"] <- 0


######add a full ant ID column
twitch_table$full_ant_ID <- with (twitch_table,paste(Replicate,treatment,group_size,ant_ID,sep="_"))

##deaths
deaths <- subset(events, Comment_start == "DEAD")
dead_ants <- with (deaths, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
twitch_table <- twitch_table [ which( ! twitch_table$full_ant_ID %in% dead_ants) ,   ]
actual_group_size <- aggregate (ant_ID ~ Replicate + treatment + group_size, function(x)length(unique(x)),data=twitch_table)
names(actual_group_size)[which(names(actual_group_size)=="ant_ID")] <- "actual_group_size"
twitch_table <- merge(twitch_table,actual_group_size,all.x=T)


###escapes
escapes <- subset(events, Comment_start == "ESCAPE")
escapes$full_ant_ID <-  with (escapes, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "duration"] <- twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "duration"] - as.numeric(escapes[match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )))],"Duration_.s."])

###total number of twitches
# total_twitches <-aggregate (cbind(frequency,N) ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
total_twitches <-aggregate (N ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
total_twitches$duration <- total_twitches$duration / 60
total_twitches$frequency <- total_twitches$N/total_twitches$duration


##add a replicate set
##to separate replicates from spring experiment and summer experiment
total_twitches$Replicate <- as.numeric(total_twitches$Replicate)
total_twitches$rep_set <- "1"
total_twitches [ which (total_twitches$Replicate > 16 ) , "rep_set" ] <- "2"


###remove summer focals
rep_set_one <- subset(total_twitches, rep_set=="1")
rep_set_two <- subset(total_twitches, rep_set=="2")

summer_nestmates <- subset(rep_set_two, ant_status=="NESTMATE")

total_twitches <- rbind(rep_set_one, summer_nestmates)


##temporal effects
total_twitches$time_frequency <- total_twitches$N / 600


total_twitches$treatment_full <- NA
total_twitches[which(total_twitches$treatment == "P" ), "treatment_full" ] <- "Pathogen"
total_twitches[which(total_twitches$treatment == "S" ), "treatment_full" ] <- "Sham"
total_twitches[which(total_twitches$treatment == "C" ), "treatment_full" ] <- "Control"
unique(total_twitches$treatment_full)


standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x))
mean_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
sd_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=sd, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
mean_status$se <- error_status$frequency
colnames(mean_status) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status$Group_size <- factor(mean_status$Group_size,
                                 levels=c("1", "2", "6", "11"))


write.csv(total_twitches, "total_twitches.csv", row.names = FALSE)


nestmate <- aggregate( frequency ~ treatment_full + ant_status, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
sd_status <- aggregate( frequency ~  treatment_full + ant_status, FUN=sd, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status <- aggregate( frequency ~ treatment_full + ant_status, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])

nestmate$se <- error_status$frequency
colnames(nestmate) <- c("Treatment", "ant_status", "Mean", "standard_error")

nestmate <- subset(nestmate, ant_status=="NESTMATE")
write.csv(nestmate, 'nestmate.csv', row.names = FALSE)








########### 5 MIN BIN ##############


##subset events to only include shake behaviours
events2 <- events[ which(events$Behavior == "T")   ,   ]

events2$petri_dish <- with ( events2, paste (Replicate, treatment, group_size, sep="_"))




##add replicate set for winter vs summer replicates
events2$rep_set <- "1"
events2 [ which (events2$Replicate > 16 ) , "rep_set" ] <- "2"

##remove summer focals from dataset
events2_focals <- subset(events2, ant_status=="FOCAL")
events2_nestmates <- subset(events2, ant_status=="NESTMATE")
events2_focals_first_rep <- subset(events2_focals, rep_set=="1")
events2 <- rbind(events2_nestmates, events2_focals_first_rep)


round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

##add time bins so events get rounded to the nearest minutes (60)
events2$Start_.s. <- as.numeric(events2$Start_.s.)
events2$time_bins <- round_any(events2$Start_.s., 60)

events2$N <- 1


events2$smins <- NA

events2[which(events2$Start_.s. <= 300 ), "smins" ] <- "0-5mins"
events2[which(events2$Start_.s. > 300 & events2$Start_.s. <= 600 ), "smins" ] <- "5-10mins"
events2[which(events2$Start_.s. > 600 & events2$Start_.s. <= 900 ), "smins" ] <- "10-15mins"
events2[which(events2$Start_.s. > 900 ), "smins" ] <- "15-20mins"






##sanity checks
shake_per_min <- aggregate(N ~ group_size + treatment + ant_status + smins, FUN=sum, data=events2)
shakes_min <- aggregate(N ~ group_size + treatment + ant_ID + time_bins + Replicate + rep_set + Observation_date,
                        FUN=sum, data=events2)



###create a data frame with all options to add zeros in
replicate_list <- rbind(unique(events[c("Replicate","treatment","group_size")]),unique(no_events[c("Replicate","treatment","group_size")]))

ant_list <- paste("ANT",c(1:11),sep="_")

shake_bloc <- expand.grid (
  smins = unique (events2$smins)
)


shake_time_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list [ 1 : replicate_list [repl, "group_size"] ] )
  repl_shake <- expand.grid.df( data.frame ( replicate_list[repl,])
                                , ants,
                                shake_bloc)
  shake_time_table <- rbind(shake_time_table, repl_shake)
}


shake_Nb <- aggregate(N ~ Replicate + treatment + group_size + ant_ID + smins + rep_set, data=events2, FUN=sum  )
shake_time_table1 <- merge(shake_time_table, shake_Nb, all.x=T)


##add in missing values
shake_time_table1 [is.na (shake_time_table1$N), "N" ] <- 0
shake_time_table1$ant_status <- "NESTMATE"
shake_time_table1 [ which ( shake_time_table1$ant_ID=="ANT_1"), "ant_status"] <- "FOCAL"
shake_time_table1$rep_set <- "1"
shake_time_table1 [ which (shake_time_table1$Replicate > 16 ) , "rep_set" ] <- "2"


##create a column to identify unique ants
shake_time_table1$full_ant_ID <-  with (shake_time_table1,paste(Replicate,treatment,group_size,ant_ID,sep="_"))
shake_time_table1$group_size <- factor(shake_time_table1$group_size,
                                       levels=c("1", "2", "6", "11"))

shake_time_table1$smins <- factor(shake_time_table1$smins, levels=c("0-5mins", "5-10mins", "10-15mins", "15-20mins"))


shake_time_table1$duration <- 1200
escapes <- subset(events, Comment_start == "ESCAPE")
escapes$full_ant_ID <-  with (escapes, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
shake_time_table1[which(!is.na(match   (   shake_time_table1$full_ant_ID,escapes$full_ant_ID  ))), "duration"] <- shake_time_table1[which(!is.na(match   (   shake_time_table1$full_ant_ID,escapes$full_ant_ID  ))), "duration"] - as.numeric(escapes[match   (   shake_time_table1$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   shake_time_table1$full_ant_ID,escapes$full_ant_ID  )))],"Duration_.s."])
shake_time_table1$duration <- shake_time_table1$duration / 60

shake_time_table1$frequency <- shake_time_table1$N / shake_time_table1$duration




mean_shake_mins <- aggregate(frequency  ~ treatment + ant_status + smins, FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~   treatment + ant_status + smins, FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c("Treatment", "ant_status", "time_bins", "mean_shakes", "se")


mean_shake_mins[which(mean_shake_mins$Treatment=="C"), "Treatment"] <- "Control"
mean_shake_mins[which(mean_shake_mins$Treatment=="S"), "Treatment"] <- "Sham"
mean_shake_mins[which(mean_shake_mins$Treatment=="P"), "Treatment"] <- "Pathogen"

write.csv(mean_shake_mins, 'mean_shake_mins.csv', row.names=FALSE)


test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}



bouts_df <- read.csv('full_table.csv', header=T, sep=',')


means = aggregate(number_of_bouts ~ group_size + ant_status, FUN=mean, data=bouts_df)

bouts_df$Replicate <- factor(bouts_df$Replicate)
bouts_df$group_size <- factor(bouts_df$group_size)
bouts_df$treatment <- factor(bouts_df$treatment)

model <- glmer(number_of_bouts ~ treatment+group_size + (1|Replicate), family = 'poisson', data=bouts_df[which(bouts_df$ant_status=="FOCAL"),])

#test_norm(residuals(model))
#hist(residuals(model))


Anova(model, type='III')

summary(glht(model, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))


model <- glmer(number_of_bouts ~ treatment+group_size + (1|Replicate), family = 'poisson', data=bouts_df[which(bouts_df$ant_status=="NESTMATE"),])

Anova(model, type='III')
summary(glht(model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
summary(glht(model, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
hist(residuals(model))


aggregate(number_of_bouts ~ treatment + ant_status, FUN=mean, data=bouts_df)




bouts <- read.csv('bouts_df.csv', header=T, sep=',')

bouts$ant_status <- "NESTMATE"
bouts [ which ( bouts$Subject=="ANT_1"), "ant_status"] <- "FOCAL"


bouts$Replicate <- factor(bouts$Replicate)
bouts$group_size <- factor(bouts$group_size)
bouts$treatment <- factor(bouts$treatment)

## focals
model <- glmer(shake_count ~ group_size*treatment + (1|Replicate), family = 'poisson', 
             data=bouts[which(bouts$ant_status=="FOCAL"),] )

#test_norm(residuals(model))
#hist(residuals(model))

Anova(model, type='III')

summary(glht(model, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))


means <- aggregate(shake_count ~ treatment  + ant_status, FUN=mean, data=bouts)
write.csv(means, 'treatment_means.csv', row.names=FALSE)



means <- aggregate(shake_count ~ treatment + group_size  + ant_status, FUN=mean, data=bouts)



write.csv(means, 'nb_bouts_treatment.csv', row.names=FALSE)


bouts

bouts$interac_size_treat <- with (bouts, interac_size_treat <- paste(group_size,treatment,sep="_"))

model <- glmer(shake_count ~ interac_size_treat + (1|Replicate), family = 'poisson', 
               data=bouts[which(bouts$ant_status=="FOCAL"),] )

Anova(model, type='II')

summary(glht(model, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))





model <- glmer(shake_count ~ group_size+treatment + (1|Replicate), family = 'poisson', 
               data=bouts[which(bouts$ant_status=="NESTMATE"),] )





standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x))
means = aggregate(number_of_bouts ~ treatment + ant_status, FUN=mean, data=bouts_df)
stdev <- aggregate(number_of_bouts ~ treatment   + ant_status, FUN=sd, data=bouts_df)
se <- aggregate(number_of_bouts ~ treatment   + ant_status, FUN=standard_error, data=bouts_df)
means$se <- se$number_of_bouts






### CHAIN TRANSMISSION ###

data = read.csv('final_chains_removal.csv', header=TRUE, sep=',')
data

data$N <- 1

data$Observation.id <- factor(data$Observation.id)
data$group_size <- factor(data$group_size)
data$Treatment <- factor(data$Treatment)



model <- glmer(Chain.Length ~ group_size+Treatment + (1|Observation.id), family=Gamma(link='log'),
               data=data )
summary(model)

#hist(residuals(model))

Anova(model, type='II')


summary(glht(model, linfct=mcp(Treatment="Tukey")), test=adjusted("BH"))

summary(glht(model, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))




nb_chains_data = read.csv('all_trials_nb_chains.csv', header=TRUE, sep=',')
nb_chains_data$Replicate <- factor(nb_chains_data$Replicate)
nb_chains_data$group_size <- factor(nb_chains_data$group_size)
nb_chains_data$Treatment <- factor(nb_chains_data$Treatment)


model <- glmer(nb_chains ~ group_size+Treatment + (1|Replicate), family='poisson', data=nb_chains_data)
model_nb <- glmer.nb(nb_chains ~ Treatment * group_size + (1|Replicate), data = nb_chains_data)

Anova(model, type='II')


summary(glht(model, linfct=mcp(Treatment="Tukey")), test=adjusted("BH"))

summary(glht(model, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))


aggregate(nb_chains ~ Treatment , FUN=mean, data=nb_chains_data)

hist(as.numeric(data$Chain.Length))


blue = read.csv('all_trials_nb_chains_BLUE.csv', header=TRUE, sep=',')
blue$Replicate <- factor(blue$Replicate)
blue$group_size <- factor(blue$group_size)
blue$Treatment <- factor(blue$Treatment)


model <- glmer(nb_chains ~ group_size+Treatment + (1|Replicate), family='poisson', data=blue)

Anova(model, type='III')
summary(glht(model, linfct=mcp(Treatment="Tukey")), test=adjusted("BH"))

summary(glht(model, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))


aggregate(nb_chains ~ Treatment , FUN=mean, data=blue)

hist(as.numeric(data$Chain.Length))




