rm(list=ls())


#setwd("C:/Users/bzniks/Dropbox/SeniorLectureship_Bristol/Students_postdocs/PhD_students/2021_Rachael_Brown/experiment1")
# setwd("C:/Users/racha/OneDrive/Documents/PhD/Winter 2022/Experiment/DATA ANALYSIS")
setwd("C:/Users/rb17990/OneDrive - University of Bristol/Documents")
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







#####GRAPHS#####

## package to have a font of choice in the plots
library(showtext)
font_add_google("Roboto") #choose a font from google fonts
showtext_auto() #must be called to indicate that showtext is going to be automatically invoked to draw text whenever a plot is created.



##focal vs nestmate mean twitches
standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x))
mean_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
sd_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=sd, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status <- aggregate( frequency ~ group_size + treatment_full + ant_status, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])

mean_status$se <- error_status$frequency
colnames(mean_status) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status$Group_size <- factor(mean_status$Group_size,
                                 levels=c("1", "2", "6", "11"))

# cld <- c( "a", "b", "a", "ab", "d", "c", "a", "e", "e", "a", "d", "c", "", "", "", "", "", "", "", "", "")
# cld <- c("a", "ab", "a", "a", "b", "d", "e", "d", "a", "c", "e", "c", "", "", "", "", "", "", "", "", "")
cld <- c("a", "b", "ab", "a", "c", "e", "f", "e", "ab", "d", "e", "e", "", "", "", "", "", "", "", "", "")

ggplot(mean_status, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), linewidth=1) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label= cld, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=5, family="Roboto") +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Mean number body shakes per ant per minute")+
  ylim(0,0.3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=18),
        legend.title=element_text(size=18), legend.text=element_text(size=16), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16), legend.position = "bottom", strip.background = element_blank(), strip.text.x = element_text(size=18, family="Roboto") )



nestmate <- aggregate( frequency ~ treatment_full + ant_status, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
sd_status <- aggregate( frequency ~  treatment_full + ant_status, FUN=sd, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status <- aggregate( frequency ~ treatment_full + ant_status, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])

nestmate$se <- error_status$frequency
colnames(nestmate) <- c("Treatment", "ant_status", "Mean", "standard_error")

nestmate <- subset(nestmate, ant_status=="NESTMATE")

letters <- c("a", "c", "b")


ggplot(nestmate, aes(fill= Treatment, y=Mean, x=Treatment))+
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.6, position=position_dodge(0.9), size=1, linewidth=1.4) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label= letters, y = Mean + standard_error), vjust= -1, position = position_dodge(0.9), size=6, family="Roboto") +
  labs(x="Treatment", y="Mean number body shakes per ant per minute", title="Number of body shakes performed by nestmates")+
  theme_bw()+
  guides(fill="none")+
  ylim(c(0,0.04 ))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=18),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=16, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )






# mean_status_control <- mean_status
# mean_status_control$control <- mean_status$Mean
# mean_status_control[which(mean_status_control$Treatment!="Control"), "control" ] <- 0
# mean_status_control$control_se <- mean_status$standard_error
# mean_status_control[which(mean_status_control$Treatment!="Control"), "control_se" ] <- 0
# 
# ggplot(mean_status_control, aes(fill=Treatment, y=control, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=control-control_se, ymax=control+control_se), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9), size=5) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")+
#   ylim(0,0.01)+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=17),
#         legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
#         axis.title.y=element_text(size=16) )
# 
# mean_status_control <- mean_status
# mean_status_control$control <- mean_status$Mean
# mean_status_control[which(mean_status_control$Treatment=="Pathogen"), "control" ] <- 0
# mean_status_control$control_se <- mean_status$standard_error
# mean_status_control[which(mean_status_control$Treatment=="Pathogen"), "control_se" ] <- 0
# 
# ggplot(mean_status_control, aes(fill=Treatment, y=control, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=control-control_se, ymax=control+control_se), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label= cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9), size=5) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")+
#   ylim(0,0.01)+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=17),
#         legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
#         axis.title.y=element_text(size=16) )
# 

##hist of shakes per trial
total_twitches$petri_dish <- with( total_twitches, paste (Replicate, treatment, group_size, sep="_"))

mean_status <- aggregate( N ~ petri_dish, FUN=sum, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
ggplot(mean_status,aes(x=N))+
  geom_histogram(binwidth=5, colour="black", fill="white")+
  labs(title="Number of body shakes per trial", x="total number of shakes", y="trial count")



events
##boxplot of frequency of twitches
boxplot1 <- total_twitches

boxplot1$group_size <- factor(boxplot1$group_size,
                              levels=c("1", "2", "6", "11"))

ggplot(boxplot1, aes(fill=treatment, y=N, x=group_size, )) +
  geom_boxplot(alpha=0.3) +
  labs ( x = "Group size", y = "Frequency of body shakes" , fill = "Treatment" ) +
  scale_fill_discrete(labels = c("Control", "Pathogen", "Sham"))+
  facet_wrap(~ant_status)

# ggplot(boxplot1, aes(x=N )) +
#   geom_histogram(binwidth = 2) +
#   labs ( x = "Group size", y = "Frequency of body shakes" , fill = "Treatment" ) +
#   scale_fill_discrete(labels = c("Control", "Pathogen", "Sham"))+
#   facet_wrap(~ant_status)

#frequency polygon of twitch decay
scatter <- events1
scatter$N <- 1
scatter1 <- cbind(as.numeric(scatter$N), as.numeric(scatter$Start_.s.), scatter$group_size, scatter$treatment, scatter$ant_status, scatter$ant_ID, scatter$Observation_id)
scatter1 <- data.frame(scatter1)
colnames(scatter1) <- c("N", "Time", "Group_size", "Treatment", "Ant_status", "Ant_ID", "Observation_ID")

scatter1$N <- as.numeric(scatter1$N)
scatter1$Time <- as.numeric(scatter1$Time)


scatter1$Group_size <- factor(scatter1$Group_size, 
                              levels = c("1", "2", "6", "11"))

scatter1_focal <- subset(scatter1, Ant_status=="FOCAL")
scatter1_nestmate <- subset(scatter1, Ant_status=="NESTMATE")

# ggplot(scatter1, aes(Time, colour=Treatment)) +
#   geom_freqpoly(binwidth=60, linewidth=1) +
#   xlab("Time (s)") + ylab("Body shake count") +
#   xlim(0, 1200) + facet_wrap( ~ Group_size)

# ggplot(scatter1, aes(Time, colour=Treatment)) +
#   geom_histogram(binwidth=60, linewidth=1) +
#   xlab("Time (s)") + ylab("Body shake count") +
#   xlim(0, 1200) + facet_wrap( ~ Group_size)

ggplot(scatter1_focal, aes(Time, colour=Treatment))+
  geom_freqpoly(binwidth=30, linewidth=1)+
  xlab("Time (s)") + ylab("Body shake count") + labs(subtitle = "Focal body shakes") +
  xlim(0,1200) + ylim(0,80) + facet_wrap ( ~ Group_size)


ggplot(scatter1_nestmate, aes(Time, colour=Treatment))+
  geom_freqpoly(binwidth=30, linewidth=1)+
  xlab("Time (s)") + ylab("Body shake count") + labs(subtitle="Nestmate body shakes")+
  xlim(0, 1200) + ylim(0, 80) + facet_wrap ( ~ Group_size)


##body shakes in first and second 10 minute periods
# mean_status_time <- aggregate( time_frequency ~ group_size + treatment + ant_status + time, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# error_status_time <- aggregate( frequency ~ group_size + treatment + ant_status + time, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# mean_status_time$se <- error_status_time$frequency
# colnames(mean_status_time) <- c("Group_size", "Treatment", "ant_status", "time", "Mean", "standard_error")
# mean_status_time$Group_size <- factor(mean_status_time$Group_size,
#                                       levels=c("1", "2", "6", "11"))
# 
# first <- subset(mean_status_time, time=="B")
# second <- subset(mean_status_time, time=="A")
# 
# 
# ggplot(first, aes(fill=Treatment, y=Mean, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label =x cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")+
#   ggtitle("body shakes in first 10 minutes")+
#   ylim(0, 0.03)
# 
# 
# 
# ggplot(second, aes(fill=Treatment, y=Mean, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label =x cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")+
#   ggtitle("body shakes in second 10 minutes") +
#   ylim(0, 0.03)


##scatter plots with regression line 
##number of shakes across trial
#]time bin width of 1 minute
events2 <- events1
events2$rep_set <- "1"
events2 [ which (events2$Replicate > 16 ) , "rep_set" ] <- "2"
events2_focals <- subset(events2, ant_status=="FOCAL")
events2_nestmates <- subset(events2, ant_status=="NESTMATE")
events2_focals_first_rep <- subset(events2_focals, rep_set=="1")
events2 <- rbind(events2_nestmates, events2_focals_first_rep)

events2$Start_.s. <- as.numeric(events2$Start_.s.)

events2$time_bins <- round_any(events2$Start_.s., 60)
shake_per_min <- aggregate(N ~ group_size + treatment + ant_status + time_bins, FUN=sum, data=events2)

ggplot(shake_per_min, aes(x=time_bins, y=N, color=treatment))+
  geom_point()+
  geom_smooth(data=subset(shake_per_min, treatment=="S"), method="lm", formula=y~x+I(x^2), colour="blue")+
  geom_smooth(data=subset(shake_per_min, treatment=="P"), method="lm", formula=y~x+I(x^2), colour="green")+
  geom_smooth(data=subset(shake_per_min, treatment=="C"), method="lm", formula=y~x+I(x^2), colour="red")+
  scale_x_continuous(name="time", breaks=seq(0,1200,300), limits=c(0,1200))+
  ylab("Number of body shakes") +
  facet_wrap(~ant_status)


### number of twitches of each type
twitch_table$twitch_type <- with ( twitch_table, paste (twitch_size, twitch_repetition, sep="_"))
twitches_per_type <-aggregate (N ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
twitches_per_type$frequency <- twitches_per_type$N/twitches_per_type$duration

# ##barplots of twitch types in group sizes and treatments
# ##some zeros in data
# mean_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = mean, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood"))   ] )
# standard_dev_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = sd, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood"))   ] )
# error_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = standard_error, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood"))   ] )
# 
# mean_type$se <- error_type$frequency
# colnames(mean_type) <- c("group_size", "treatment", "twitch_type", "mean", "standard_error")
# mean_type$group_size <- as.character(mean_type$group_size)
# mean_type$group_size <- factor(mean_type$group_size,
#                                levels=c("1", "2", "6", "11"))
# 
# ggplot(mean_type, aes(fill=treatment, y=mean, x=group_size))+
#   geom_bar(position="dodge", stat="identity")+
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
#   facet_wrap(~twitch_type)


##just repetitive
total_twitches$frequency
mean_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition + ant_status, FUN = mean, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )
standard_dev_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition + ant_status, FUN = sd, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )
error_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition + ant_status, FUN = standard_error, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )

mean_rep$se <- error_rep$frequency
colnames(mean_rep) <- c("group_size", "Treatment", "repetition", "ant_status", "mean", "standard_error")
mean_rep$group_size <- factor(mean_rep$group_size,
                              levels=c("1", "2", "6", "11"))


mean_rep$status_rep <- with (mean_rep, paste(repetition,ant_status,sep="_"))
mean_rep$repetition <- as.character(mean_rep$repetition)

mean_rep[which(mean_rep$Treatment == "C") , "Treatment" ] <- "Control"
mean_rep[which(mean_rep$Treatment == "S") , "Treatment" ] <- "Sham"
mean_rep[which(mean_rep$Treatment == "P") , "Treatment" ] <- "Pathogen"

mean_rep[which(mean_rep$status_rep == "single_FOCAL"), "status_rep"] <- "FOCAL: Single shakes"
mean_rep[which(mean_rep$status_rep == "repetitive_FOCAL"), "status_rep"] <- "FOCAL: Repetitive shakes"
mean_rep[which(mean_rep$status_rep == "single_NESTMATE"), "status_rep"] <- "NESTMATES: Single shakes"
mean_rep[which(mean_rep$status_rep == "repetitive_NESTMATE"), "status_rep"] <- "NESTMATES: Repetitive shakes"

mean_rep[which(mean_rep$repetition ==  "single"), "repetition" ] <- "SINGLE"
mean_rep[which(mean_rep$repetition ==  "repetitive"), "repetition" ] <- "REPETITIVE"

mean_rep$repetition <- factor(mean_rep$repetition,
                              levels = c("SINGLE", "REPETITIVE"))


ggplot(mean_rep, aes(fill=Treatment, y=mean, x=group_size))+
  geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.6, position=position_dodge(0.9), size=1)+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(repetition ~ ant_status) +
  theme_bw()+
  labs(x="Group size", y="Number of body shakes per ant per minute")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )


##just size
total_twitches$frequency
mean_size <- aggregate ( frequency ~ group_size + treatment + twitch_size + ant_status, FUN = mean, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_repetition"))   ] )
standard_dev_rep <- aggregate ( frequency ~ group_size + treatment + twitch_size + ant_status, FUN = sd, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_repetition"))   ] )
error_rep <- aggregate ( frequency ~ group_size + treatment + twitch_size + ant_status, FUN = standard_error, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_repetition"))   ] )

mean_size$se <- error_rep$frequency
colnames(mean_size) <- c("group_size", "Treatment", "size", "ant_status", "mean", "standard_error")
mean_size$group_size <- factor(mean_size$group_size,
                               levels=c("1", "2", "6", "11"))


mean_size$status_rep <- with (mean_size, paste(size,ant_status,sep="_"))
mean_size$size <- as.character(mean_size$size)

mean_size[which(mean_size$Treatment == "C") , "Treatment" ] <- "Control"
mean_size[which(mean_size$Treatment == "S") , "Treatment" ] <- "Sham"
mean_size[which(mean_size$Treatment == "P") , "Treatment" ] <- "Pathogen"

mean_size[which(mean_size$status_rep == "large_FOCAL"), "status_rep"] <- "FOCAL: large shakes"
mean_size[which(mean_size$status_rep == "small_FOCAL"), "status_rep"] <- "FOCAL: small shakes"
mean_size[which(mean_size$status_rep == "large_NESTMATE"), "status_rep"] <- "NESTMATES: large shakes"
mean_size[which(mean_size$status_rep == "small_NESTMATE"), "status_rep"] <- "NESTMATES: small shakes"

mean_size[which(mean_size$size ==  "large"), "size" ] <- "LARGE"
mean_size[which(mean_size$size ==  "small"), "size" ] <- "SMALL"

mean_size$size <- factor(mean_size$size,
                         levels = c("SMALL", "LARGE"))


ggplot(mean_size, aes(fill=Treatment, y=mean, x=group_size))+
  geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.6, position=position_dodge(0.9), size=1)+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(size ~ ant_status) +
  theme_bw()+
  labs(x="Group size", y="Number of body shakes per ant per minute")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )





# ###beeswarm plots
# sum_status <- aggregate(N ~ group_size + treatment + ant_status + full_ant_ID, FUN=sum, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# sum_status_int <- subset(sum_status, sum_status$N > 0)
# 
# ggplot(sum_status_int, aes(fill=treatment, y=N, x=group_size, color=treatment))+
#   geom_beeswarm(cex=1, priority = "density") +
#   labs(x="group size", y="Number of body shakes of each ant") +
#   facet_wrap(~ant_status)
# 
# 
# #just focal
# sum_status_int_focal <- subset(sum_status_int, sum_status_int$ant_status=="FOCAL")
# 
# ggplot(sum_status_int_focal, aes(y=N, x=treatment, color=treatment))+
#   geom_beeswarm(cex=2, priority = "density") +
#   labs(x="group size", y="Number of body shakes of each ant") +
#   facet_wrap(~group_size)
# 
# #just nestmates
# sum_status_int_nestmate <- subset(sum_status_int, sum_status_int$ant_status=="NESTMATE")
# 
# ggplot(sum_status_int_nestmate, aes(y=N, x=treatment, color=treatment))+
#   geom_beeswarm(cex=1, priority = "density") +
#   labs(x="group size", y="Number of body shakes of each ant") +
#   facet_wrap(~group_size)


##focal vs nestmate size and repetition
#focal
focal_mean_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = mean, data = twitches_per_type[ which(twitches_per_type$ant_status=="FOCAL"),  ] )
focal_error_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = standard_error, data =  twitches_per_type[ which(twitches_per_type$ant_status=="FOCAL"),  ]  )

focal_mean_type$se <- focal_error_type$frequency
colnames(focal_mean_type) <- c("group_size", "treatment", "twitch_type", "mean", "standard_error")
focal_mean_type$group_size <- as.character(focal_mean_type$group_size)
focal_mean_type$group_size <- factor(focal_mean_type$group_size,
                                     levels=c("1", "2", "6", "11"))

ggplot(focal_mean_type, aes(fill=treatment, y=mean, x=group_size))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
  facet_wrap(~twitch_type)



#nestmate
nestmate_mean_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = mean, data =  twitches_per_type[ which(twitches_per_type$ant_status=="NESTMATE"),  ] )
nestmate_error_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = standard_error, data =  twitches_per_type[ which(twitches_per_type$ant_status=="NESTMATE"),  ] )

nestmate_mean_type$se <- nestmate_error_type$frequency
colnames(nestmate_mean_type) <- c("group_size", "treatment", "twitch_type", "mean", "standard_error")
nestmate_mean_type$group_size <- as.character(nestmate_mean_type$group_size)
nestmate_mean_type$group_size <- factor(nestmate_mean_type$group_size,
                                        levels=c("1", "2", "6", "11"))

ggplot(nestmate_mean_type, aes(fill=treatment, y=mean, x=group_size))+
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
  facet_wrap(~twitch_type)



##new replicates barplot
##checking data from summer experiment
new_replicates <- subset(total_twitches, total_twitches$time == "SUMMER")

mean_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=mean, data= new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
sd_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=sd, data=new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=standard_error, data=new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])

mean_status_new$se <- error_status_new$frequency
colnames(mean_status_new) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status_new$Group_size <- factor(mean_status_new$Group_size,
                                     levels=c("1", "2", "6", "11"))

ggplot(mean_status_new, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
  #geom_text(aes(label = cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Number of body shakes per ant per minute")  +
  
  
  #compare with first replicates
  old_replicates <- subset(total_twitches, total_twitches$time== "WINTER")

mean_status_first <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=mean, data=old_replicates[ which(!names(old_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
sd_status_first <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=sd, data=old_replicates[ which(!names(old_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
error_status_first <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=standard_error, data=old_replicates[ which(!names(old_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])

mean_status_first$se <- error_status_first$frequency
colnames(mean_status_first) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
mean_status_first$Group_size <- factor(mean_status_first$Group_size,
                                       levels=c("1", "2", "6", "11"))


ggplot(mean_status_first, aes(fill=Treatment, y=Mean, x=Group_size)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
  #geom_text(aes(label = cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
  facet_wrap(~ant_status) +
  labs(x="Group size", y="Number of body shakes per ant per minute") 








###############
#####STATS#####
###############

## call this function form EXP1_FUNCTIONS

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

head(total_twitches)
hist(total_twitches$frequency)
100*sum(total_twitches$frequency == 0)/nrow(total_twitches) #85% of data is zeros
str(total_twitches)

##re-structuring data
total_twitches$interac_size_treat <- with (total_twitches, interac_size_treat <- paste(group_size,treatment,sep="_"))
total_twitches$Replicate <- factor(total_twitches$Replicate)
total_twitches$group_size <- factor(total_twitches$group_size)
total_twitches$treatment <- factor(total_twitches$treatment)
total_twitches$Petri_dish <- with (total_twitches, interac_size_treat <- paste(Replicate,group_size,treatment,sep="_"))

# #which replicate has the big focal twitchers in 11s
# ##replicate 19 and 21
# big_twitchers <- total_twitches [ which (total_twitches$interac_size_treat=="11_S"),]
# big_twitchers_ID <- aggregate( N ~ full_ant_ID, FUN=sum, data=big_twitchers)
# big_twitchers_ID[order(big_twitchers_ID$N),]

matrix_contrasts <- rbind( "1C-1S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 , 0),
                           "1C-1P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "1S-1P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 , 0 , 0 , 0 , 0),
                           "2C-2S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 ,-1,  0 , 0),
                           "2C-2P"   = c(0 , 0 , 0 , 0 ,-1,  0,-1,  0 , 0 , 0 , 0 , 0),
                           "2S-2P"   = c(0 , 0 , 0 , 0 ,-1,  1,-1,  0 , 0 , 1 , 0 , 0),
                           "6C-6S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 ,-1,  0),
                           "6C-6P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 ,-1,  0 , 0 , 0 , 0),
                           "6S-6P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 ,-1,  0 , 0 , 1 , 0),
                           "11C-11S" = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 ,-1),
                           "11C-11P" = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 ,-1,  0 , 0 , 0),
                           "11S-11P" = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 ,-1,  0 , 0 , 1),
                           "1C-2C"   = c(0 ,-1,  0 , 0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "1C-6C"   = c(0 , 0 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "1C-11C"  = c(0 , 0 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "2C-6C"   = c(0 , 1 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "2C-11C"  = c(0 , 1 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
                           "6C-11C"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  0,  0),
                           "1S-2S"   = c(0, -1,  0,  0,  0,  0, 0,  0,  0, -1,  0  ,0),
                           "1S-6S"   = c(0,  0, -1,  0,  0,  0, 0,  0,  0,  0, -1,  0),
                           "1S-11S"  = c(0,  0,  0, -1,  0,  0, 0,  0,  0,  0,  0, -1),
                           "2S-6S"   = c(0,  1, -1,  0,  0,  0, 0,  0,  0,  1, -1,  0),
                           "2S-11S"  = c(0,  1,  0, -1,  0,  0, 0,  0,  0,  1,  0, -1),
                           "6S-11S"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  1, -1),
                           "1P-2P"   = c(0, -1,  0,  0,  0,  0,-1,  0,  0,  0,  0,  0),
                           "1P-6P"   = c(0,  0, -1,  0,  0,  0, 0, -1,  0,  0,  0,  0),
                           "1P-11P"  = c(0,  0,  0, -1,  0,  0, 0,  0, -1,  0,  0,  0),
                           "2P-6P"   = c(0,  1, -1,  0,  0,  0, 1, -1,  0,  0,  0,  0),
                           "2P-11P"  = c(0,  1,  0, -1,  0,  0, 1,  0, -1,  0,  0,  0),
                           "6P-11P"  = c(0,  0,  1, -1,  0,  0, 0,  1, -1,  0,  0,  0))



matrix_contrasts_nestmates <- rbind("2C-2S"   = c(0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,-1 ,0 ,0 ),
                                    "2C-2P"   = c(0 ,0 ,0 ,0,-1 , 0 ,-1,0 ,0 ,0 ,0 ,0 ),
                                    "2S-2P"   = c(0 ,0 ,0 ,0,-1,1 ,-1,0 ,0 ,1 ,0 ,0 ),
                                    "6C-6S"   = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,-1,0 ),
                                    "6C-6P"   = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,-1,0 ,0 ,0 ,0 ),
                                    "6S-6P"   = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,-1,0 ,0 ,1 ,0 ),
                                    "11C-11S" = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,0 ,-1),
                                    "11C-11P" = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,-1,0 ,0 ,0 ),
                                    "11S-11P" = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,0 ,-1,0 ,0 ,1 ),
                                    "2C-6C"   = c(0 ,1,-1, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
                                    "2C-11C"  = c(0 ,1 ,0,-1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
                                    "6C-11C"  = c(0, 0, 1,-1,0,0,0,0,0,0,0,0),
                                    "2S-6S"   = c(0, 1,-1, 0,0,0,0,0,0,1,-1,0),
                                    "2S-11S"  = c(0, 1, 0,-1,0,0,0,0,0,1,0,-1),
                                    "6S-11S"  = c(0, 0, 1,-1,0,0,0,0,0,0,1,-1),
                                    "2P-6P"   = c(0, 1,-1, 0,0,0,1,-1,0,0,0,0),
                                    "2P-11P"  = c(0, 1, 0,-1,0,0,1,0,-1,0,0,0),
                                    "6P-11P"  = c(0, 0, 1,-1,0,0,0,1,-1,0,0,0))





############################
### FRREQUENCY OF SHAKES ###
############################

### FOCAL MODEL ###

###Poisson model
model_focal <- glmer(log(N+1/sqrt(2)) ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
                     family="gaussian", data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])


# qqnorm(residuals(model_focal))
# qqline(residuals(model_focal))
# hist(residuals(model_focal))

test_norm(residuals(model_focal))


Anova(model_focal, type="III") 


#interaction significant
model_focal_interac <- glmer(log(N+1/sqrt(2)) ~ interac_size_treat + offset(log(duration)) + (1|Replicate),
                             family="gaussian", data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])

# test_norm(residuals(model_focal_interac))

# hist(residuals(model_focal_interac))
# qqnorm(residuals(model_focal_interac))
# qqline(residuals(model_focal_interac))

post_hoc_focal_interac_tukey<- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# post_hoc_focal_interac <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))

print(cld(post_hoc_focal_interac_tukey))








### NESTMATE MODEL ###


model_nestmate_ni_log <- glmer(N^0.1 ~ group_size * treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                               family="gaussian", data=total_twitches[which(total_twitches$ant_status=="NESTMATE") ,])

test_norm(residuals(model_nestmate_ni_log))


Anova(model_nestmate_ni_log, type="III")


model_nestmate_ni_log <- glmer(N^0.1 ~ group_size + treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                               family="gaussian", data=total_twitches[which(total_twitches$ant_status=="NESTMATE") ,])


Anova(model_nestmate_ni_log)

post_hoc_nestmate_treat_log <- summary(glht(model_nestmate_ni_log, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
print(cld(post_hoc_nestmate_treat_log))
#  C   P   S 
# "a" "c" "b"










############################
### CATEGORIES OF SHAKES ###
############################

### FOCAL MODELS ###

### REPETITIVE ###

### TOTAL NUMBER OF REPETITVE SHAKES PERFORMED BY FOCAL ANTS

total_twitches_repetitive <- subset(total_twitches, twitch_repetition=="repetitive")


model_focal_rep <- glmer(log(N+1/sqrt(2)) ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
                         family="gaussian", data=total_twitches_repetitive[which(total_twitches_repetitive$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_rep))



qqnorm(residuals(model_focal_rep))
qqline(residuals(model_focal_rep))
hist(residuals(model_focal_rep))

Anova(model_focal_rep, type="III")



#significant interaction 

model_focal_rep_ni <- glmer(log(N+1/sqrt(2)) ~ interac_size_treat + offset(log(duration)) + (1|Replicate),
                            family="gaussian", data=total_twitches_repetitive[which(total_twitches_repetitive$ant_status=="FOCAL"),])


Anova(model_focal_rep_ni, type="II")

test_norm(residuals(model_focal_rep_ni))


h <- summary(glht(model_focal_rep_ni, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(h))

# 1_C  1_P  1_S 11_C 11_P 11_S  2_C  2_P  2_S  6_C  6_P  6_S 
# "a"  "a"  "a"  "a" "bc"  "b"  "a" "bc"  "c"  "a"  "d" "cd" 





### SINGLE ###

### TOTAL NUMBER OF SINGLE SHAKES PERFORMED BY FOCAL ANTS

total_twitches_single <- subset(total_twitches, twitch_repetition=="single")


model_focal_single <- glmer(sqrt(N) ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
                            family="gaussian", data=total_twitches_single[which(total_twitches_single$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_single))


qqnorm(residuals(model_focal_single))
qqline(residuals(model_focal_single))
hist(residuals(model_focal_single))

Anova(model_focal_single, type="III")



#interaction is significant
model_focal_single_interac <- glmer(sqrt(N) ~ interac_size_treat + offset(log(duration)) + (1|Replicate),
                                    family="gaussian", data=total_twitches_single[which(total_twitches_single$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_single_interac))


model_focal_single_interac_tukey<- summary(glht(model_focal_single_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(model_focal_single_interac_tukey)) 
# 1_C  1_P  1_S 11_C 11_P 11_S  2_C  2_P  2_S  6_C  6_P  6_S 
# "a"  "b" "ab"  "a" "de" "cd" "ab"  "c"  "c" "ab"  "e"  "e" 






##### SIZE #####



### SMALL ###

### TOTAL NUMBER OF SMALL SHAKES PERFORMED BY FOCAL ANT

total_twitches_small <- subset(total_twitches, twitch_size =="small")

model_focal_small <- glmer(sqrt(N) ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
                           family="gaussian", data=total_twitches_small[which(total_twitches_small$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_small))



qqnorm(residuals(model_focal_small))
qqline(residuals(model_focal_small))
hist(residuals(model_focal_small))

Anova(model_focal_small, type="III")

#interact sig
model_focal_small_ni <- glmer(sqrt(N) ~ interac_size_treat + offset(log(duration)) + (1|Replicate),
                           family="gaussian", data=total_twitches_small[which(total_twitches_small$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_small_ni))


Anova(model_focal_small_ni)


post_hoc <- summary(glht(model_focal_small_ni, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(post_hoc))
# 1_C  1_P  1_S 11_C 11_P 11_S  2_C  2_P  2_S  6_C  6_P  6_S 
# "a" "bc" "ab" "ab" "ef" "cd" "ab" "de"  "d" "ab"  "f"  "e" 






### LARGE ###


### TOTAL NUMBER OF LARGE SHAKES PERFORMED BY FOCAL ANT

total_twitches_large <- subset(total_twitches, twitch_size =="large")


model_focal_large <- glmer(log(N+1/sqrt(2)) ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
                           family="gaussian", data=total_twitches_large[which(total_twitches_large$ant_status=="FOCAL"),])

test_norm(residuals(model_focal_large))



qqnorm(residuals(model_focal_large))
qqline(residuals(model_focal_large))
hist(residuals(model_focal_large))

Anova(model_focal_large, type="III")

#interaction significant

model_focal_interac_large <- glmer(log(N+1/sqrt(2)) ~ interac_size_treat + offset(log(duration)) + (1|Replicate),
                                   family="gaussian", data=total_twitches_large[which(total_twitches_large$ant_status=="FOCAL"),])





model_focal_large_interac_tukey<- summary(glht(model_focal_interac_large, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
print(cld(model_focal_large_interac_tukey))

# 1_C  1_P  1_S 11_C 11_P 11_S  2_C  2_P  2_S  6_C  6_P  6_S 
# "a"  "a"  "a"  "a"  "b"  "b"  "a"  "b"  "b"  "a"  "c"  "c" 








#### NESTMATE MODELS ####



### REPETITIVE ###

### TOTAL NUMBER OF REPETITVE SHAKES PERFORMED BY NESTMATES
##issues with this model
## high kurtosis


model_nestmate_repetitive <- glmer(N^0.1 ~ group_size*treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                                   family="gaussian", data=total_twitches_repetitive[which(total_twitches_repetitive$ant_status=="NESTMATE"),])

test_norm(residuals(model_nestmate_repetitive))



qqnorm(residuals(model_nestmate_repetitive))
qqline(residuals(model_nestmate_repetitive))
hist(residuals(model_nestmate_repetitive))

Anova(model_nestmate_repetitive, type="III")

#interaction non significant
model_nestmate_repetitive_ni <- glmer(N^0.1 ~ group_size+treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                                      family="gaussian", data=total_twitches_repetitive[which(total_twitches_repetitive$ant_status=="NESTMATE"),])



Anova(model_nestmate_repetitive_ni, type="II")

summary(glht(model_nestmate_repetitive_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))








### SINGLE ###

### TOTAL NUMBER OF SINGLE SHAKES PERFORMED BY NESTMATES

model_nestmate_single <- glmer(N^0.1 ~ group_size*treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                               family="gaussian", data=total_twitches_single[which(total_twitches_single$ant_status=="NESTMATE"),])


test_norm(residuals(model_nestmate_single))
# no transformation: skew = 2.52, kurt = 11.96
# sqrt: skew=1.6, kurt=3.10
# log: skew=1.66 ,kurt=3.28


qqnorm(residuals(model_nestmate_single))
qqline(residuals(model_nestmate_single))
hist(residuals(model_nestmate_single))

Anova(model_nestmate_single, type="III")

#interaction non-sig
model__nestmate_single_ni <- glmer(N^0.1  ~ group_size+treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                                   family="gaussian", data=total_twitches_single[which(total_twitches_single$ant_status=="NESTMATE"),])


Anova(model__nestmate_single_ni, type="II")

post_hoc_nestmate <- summary(glht(model__nestmate_single_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))


print(cld(post_hoc_nestmate))







##### SIZE #####



### SMALL ###

### TOTAL NUMBER OF SMALL SHAKES PERFORMED BY NESTMATES


model_nestmate_small <- glmer(N^0.1 ~ group_size*treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                              family="gaussian", data=total_twitches_small[which(total_twitches_small$ant_status=="NESTMATE"),])

test_norm(residuals(model_nestmate_small))



qqnorm(residuals(model_nestmate_small))
qqline(residuals(model_nestmate_small))
hist(residuals(model_nestmate_small))

Anova(model_nestmate_small, type="III")

#interact non sig
model_nestmate_small_ni <- glmer(N^0.1 ~ group_size+treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                                 family="gaussian", data=total_twitches_small[which(total_twitches_small$ant_status=="NESTMATE"),])

test_norm(residuals(model_nestmate_small_ni))


Anova(model_nestmate_small_ni)


post_hocs <- summary(glht(model_nestmate_small_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))


print(cld(post_hocs))
#   C   P   S 
# "a" "c" "b" 




### LARGE ###


### TOTAL NUMBER OF LARGE SHAKES PERFORMED BY NESTMATES

model_nestmate_large <- glmer(N^0.1 ~ group_size*treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                              family="gaussian", data=total_twitches_large[which(total_twitches_large$ant_status=="NESTMATE"),])


test_norm(residuals(model_nestmate_large))



qqnorm(residuals(model_nestmate_large))
qqline(residuals(model_nestmate_large))
hist(residuals(model_nestmate_large))

Anova(model_nestmate_large, type="III")


#interaction non-significant
model_nestmate_large_ni <- glmer(N^0.1 ~ group_size+treatment + offset(log(duration)) + (1|Replicate) + (1|time) + (1|Petri_dish),
                                 family="gaussian", data=total_twitches_large[which(total_twitches_large$ant_status=="NESTMATE"),])


test_norm(residuals(model_nestmate_large_ni))


Anova(model_nestmate_large_ni, type="II")


post_hocs <- summary(glht(model_nestmate_large_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
 

print(cld(post_hocs))







# ### FULL MODEL ATTEMPT
# 
# model_focal_cat <- glmer(N ~ group_size*treatment*twitch_size*twitch_repetition + offset(log(duration)) + (1|Replicate) + (1|time),
#                      family="poisson", data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# 
# post_hoc_focal<- summary(glht(model_focal_cat, linfct=mcp(twitch_repetition="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal)) #makes a little bit of sense
# 
# qqnorm(residuals(model_focal_cat))
# qqline(residuals(model_focal_cat))
# hist(residuals(model_focal_cat))
# 
# Anova(model_focal_cat, type="III")
# 
# 
# model_focal_cat_sqrt <- glmer(sqrt(N) ~ group_size*treatment*twitch_size*twitch_repetition + offset(log(duration)) + (1|Replicate) + (1|time),
#                          family="poisson", data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# 
# 
# model_focal_cat_log <- glmer(log(N+1/sqrt(2)) ~ group_size*treatment*twitch_size*twitch_repetition + offset(log(duration)) + (1|Replicate) + (1|time),
#                               family="poisson", data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# 


# 
# ###OTHER MODELS TESTED###
# 
# ##negative binom model
# ##focal
# #model not zero inflated or overdispersed
# model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate), data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# check_zeroinflation(model_focal_NB)
# check_overdispersion(model_focal_NB)
# Anova(model_focal_NB)
# 
# #significant interaction
# #post hocs don't look too bad, check against graphs
# model_focal_interac <- glmer.nb(N ~ interac_size_treat + offset(log(duration)) + (1|rep_set:Replicate),
#                                 data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])
# post_hoc_focal_interac_nb <- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# post_hoc_focal_interac_nb_matrix <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))
# 
# 
# 
# ##nestmate 
# #not zero inflated or overdispersed
# ##output looks v wrong
# model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# check_zeroinflation(model_nestmate_NB)
# check_overdispersion(model_nestmate_NB)
# Anova(model_nestmate_NB)
# 
# 
# 
# 
# 
# 
# 

####DRAFTS####

# ##trying a model check I saw online, not sure if it is good or not
# #add more iterations 
# #output looks more reliable
# ss <- getME(model_nestmate_NB,c("theta", "fixef"))
# m2 <- update(model_nestmate_NB, start=ss, control=glmerControl(optCtrl = list(maxfun=2e4)))
# Anova(m2)
# 
# #interaction non sig
# #don't trust result
# model_nestmate_NB_ni <- glmer.nb(N ~ group_size+treatment + offset(log(duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# Anova(model_nestmate_NB_ni)
# 
# post_hoc_nestmate_treatment <- summary(glht(model_nestmate_NB_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# post_hoc_nestmate_group <- summary(glht(model_nestmate_NB_ni, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# 
# 
#
# dd <- getME(model_nestmate_NB_ni, c("theta", "fixef"))
# m3 <- update(model_nestmate_NB_ni, start=dd, control=glmerControl(optCtrl = list(maxfun=2e4)))
# 
# Anova(m3)
# post_hoc_nestmate_treatment <- summary(glht(m3, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))


# ##Hurdle model
# ###not working
# sum(total_twitches$N <1)
# 
# mod.hurdle <- hurdle(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate),
#                      data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# 
# hurdle_attempt_focal <- glmmTMB(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate),
#                                 data=total_twitches[which(total_twitches$ant_status=="FOCAL"),], zi=T, family=truncated_poisson)
##on new rep set only
##new replicates from summer experiment
# 
# 
# model_focal <- glmer(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate),
#                      family="poisson", data=[which($ant_status=="FOCAL"),])
# check_zeroinflation(model_focal)
# check_overdispersion(model_focal)
# Anova(model_focal)
# 
# model_focal_interac <- glmer(N ~ interac_size_treat + offset(log(duration)) + (1|rep_set:Replicate), family="poisson",
#                              data=[which($ant_status=="FOCAL") , ])
# post_hoc_focal_interac<- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# 
# 
# print(cld(post_hoc_focal_interac))

# model_nestmate_interac <- glmer(N ~ interac_size_treat + offset(log(duration)) + (1|Replicate) + (1|Petri_dish),
#                                 family="poisson", data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# Anova(model_nestmate_interac)
# post_hoc_nestmate_interac <- summary(glht(model_nestmate_interac, linfct=matrix_contrasts_nestmates), test=adjusted("BH"))

##negative binom model
# ##focal
# #model not zero inflated or overdispersed
# model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate), data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# check_zeroinflation(model_focal_NB)
# check_overdispersion(model_focal_NB)
# Anova(model_focal_NB)
# 
# #significant interaction
# #post hocs don't look too bad, check against graphs
# ###no sig differences between S and P, but sig diff between control and experimental
# model_focal_interac <- glmer.nb(N ~ interac_size_treat + offset(log(duration)) + (1|rep_set:Replicate),
#                              data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])
# post_hoc_focal_interac_nb <- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# post_hoc_focal_interac_nb_matrix <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))
# 
# 
# ##nestmate 
# #not zero inflated or overdispersed
# model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# check_zeroinflation(model_nestmate_NB)
# check_overdispersion(model_nestmate_NB)
# Anova(model_nestmate_NB)
# 
# #add more iterations 
# ss <- getME(model_nestmate_NB,c("theta", "fixef"))
# m2 <- update(model_nestmate_NB, start=ss, control=glmerControl(optCtrl = list(maxfun=2e4)))
# Anova(m2)
# 
# #now interaction non sig
# model_nestmate_NB_ni <- glmer.nb(N ~ group_size+treatment + offset(log(duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# Anova(model_nestmate_NB_ni)
# 
# post_hoc_nestmate_treatment <- summary(glht(model_nestmate_NB_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# post_hoc_nestmate_group <- summary(glht(model_nestmate_NB_ni, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# 
# # dd <- getME(model_nestmate_NB_ni, c("theta", "fixef"))
# # m3 <- update(model_nestmate_NB_ni, start=dd, control=glmerControl(optCtrl = list(maxfun=2e4)))


# #check singularity
# #not an issue
# tt <- getME(model_nestmate_NB,"theta")
# ll <- getME(model_nestmate_NB, "lower")
# min(tt[ll==0])
# 
# #double checking gradient calculations
# derivs1 <- model_nestmate_NB@optinfo$derivs
# sc_grad1 <- with(derivs1, solve(Hessian, gradient))
# max(abs(sc_grad1))
# max(pmin(abs(sc_grad1), abs(derivs1$gradient)))

# 
# #significant interaction
# model_nestmate_NB_interac <- glmer.nb(N ~ interac_size_treat + offset(log(duration)) + (1|rep_set:Replicate)+ 
#                                         (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# post_hoc_nestmate_interac_nb <- summary(glht(model_nestmate_NB_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_interac_nb))


# 
# ###interaction is not significant, proceed with simple post-hocs
# post_hoc_nestmate_treatment<- summary(glht(model_nestmate_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_treatment))
# 
# post_hoc_nestmate_size<- summary(glht(model_nestmate_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_size))



# ###interaction is not significant, proceed with simple post-hocs
# post_hoc_focal_treatment<- summary(glht(model_focal_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal_treatment))
# 
# post_hoc_focal_size<- summary(glht(model_focal_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal_size))

###very odd post-hoc outcomes where all group size come out as equal!!
##try an alternative
# zinbinom_focal<- glmmTMB(N~treatment*group_size+
#                            offset(log(duration))+(1|Replicate),
#                          data=total_twitches[which(total_twitches$ant_status=="FOCAL"),],
#                          zi=~treatment+group_size,
#                          family=nbinom2)
# Anova(zinbinom_focal)
# post_hoc_focal_treatment<- summary(glht(zinbinom_focal, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal_treatment))
# post_hoc_focal_size<- summary(glht(zinbinom_focal, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal_size))
###STRANGE POST_HOCS WITH NEGATIVE BINOMIAL GLMER!



###Poisson model is both overdispersed and zero-inflated
# model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# check_zeroinflation(model_nestmate_NB)
# check_overdispersion(model_nestmate_NB)
# ###still overdispersed
# Anova(model_nestmate_NB)
# 
# ###interaction is not significant, proceed with simple post-hocs
# post_hoc_nestmate_treatment<- summary(glht(model_nestmate_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_treatment))
# 
# post_hoc_nestmate_size<- summary(glht(model_nestmate_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_size))


##zero inflation models
# zipoisson_focal <- glmmTMB(N~treatment*group_size+
#                              offset(log(duration))+(1|Replicate),
#                            data=total_twitches[which(total_twitches$ant_status=="FOCAL"),],
#                            ziformula=~.,
#                            family=poisson)
# Anova(zipoisson_focal)


# zipoisson_focal <-glmmadmb(frequency~treatment*group_size+(1|Replicate),zeroInflation=TRUE,data=twitch_table[which(twitch_table$ant_status=="FOCAL"),],family = "poisson")
# Anova(zipoisson_focal)
# check_zeroinflation(zipoisson_focal) ##overfitting

# zipoisson_focal <- glmmadmb(N ~ group_size*treatment + offset(log(duration)) + (1|Replicate),
#                             zeroInflation = TRUE, data=twitch_table[which(twitch_table$ant_status=="FOCAL"),],family = "poisson")


# model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|Replicate), data=twitch_table[which(twitch_table$ant_status=="FOCAL"),])
# 
# 
# model_focal_interac_NB <- glmer.nb(N ~ interac_size_treat + offset(log(duration)) + (1|Replicate), family="poisson", data=twitch_table[which(twitch_table$ant_status=="FOCAL"),])
# Anova(model_focal_interac_NB)
# post_hoc_focal_interac_NB <- summary(glht(model_focal_interac_NB, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_focal_interac_NB))
#   
# matrix_contrasts <- rbind( "1C-1S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 , 0),
#                           "1C-1P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "1S-1P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "2C-2S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 ,-1,  0 , 0),
#                           "2C-2P"   = c(0 , 0 , 0 , 0 ,-1,  0,-1,  0 , 0 , 0 , 0 , 0),
#                           "2S-2P"   = c(0 , 0 , 0 , 0 ,-1,  1,-1,  0 , 0 , 1 , 0 , 0),
#                           "6C-6S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 ,-1,  0),
#                           "6C-6P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 ,-1,  0 , 0 , 0 , 0),
#                           "6S-6P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 ,-1,  0 , 0 , 1 , 0),
#                           "11C-11S" = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 ,-1),
#                           "11C-11P" = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 ,-1,  0 , 0 , 0),
#                           "11S-11P" = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 ,-1,  0 , 0 , 1),
#                           "1C-2C"   = c(0 ,-1,  0 , 0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "1C-6C"   = c(0 , 0 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "1C-11C"  = c(0 , 0 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "2C-6C"   = c(0 , 1 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "2C-11C"  = c(0 , 1 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                           "6C-11C"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  0,  0),
#                           "1S-2S"   = c(0, -1,  0,  0,  0,  0, 0,  0,  0, -1,  0  ,0),
#                           "1S-6S"   = c(0,  0, -1,  0,  0,  0, 0,  0,  0,  0, -1,  0),
#                           "1S-11S"  = c(0,  0,  0, -1,  0,  0, 0,  0,  0,  0,  0, -1),
#                           "2S-6S"   = c(0,  1, -1,  0,  0,  0, 0,  0,  0,  1, -1,  0),
#                           "2S-11S"  = c(0,  1,  0, -1,  0,  0, 0,  0,  0,  1,  0, -1),
#                           "6S-11S"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  1, -1),
#                           "1P-2P"   = c(0, -1,  0,  0,  0,  0,-1,  0,  0,  0,  0,  0),
#                           "1P-6P"   = c(0,  0, -1,  0,  0,  0, 0, -1,  0,  0,  0,  0),
#                           "1P-11P"  = c(0,  0,  0, -1,  0,  0, 0,  0, -1,  0,  0,  0),
#                           "2P-6P"   = c(0,  1, -1,  0,  0,  0, 1, -1,  0,  0,  0,  0),
#                           "2P-11P"  = c(0,  1,  0, -1,  0,  0, 1,  0, -1,  0,  0,  0),
#                           "6P-11P"  = c(0,  0,  1, -1,  0,  0, 0,  1, -1,  0,  0,  0))

# matrix_contrasts_nestmates <- rbind("2C-2S"   = c(0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,-1 ,0 ,0 ),
#                                     "2C-2P"   = c(0 ,0 ,0 ,0,-1 , 0 ,-1,0 ,0 ,0 ,0 ,0 ),
#                                     "2S-2P"   = c(0 ,0 ,0 ,0,-1,1 ,-1,0 ,0 ,1 ,0 ,0 ),
#                                     "6C-6S"   = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,-1,0 ),
#                                     "6C-6P"   = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,-1,0 ,0 ,0 ,0 ),
#                                     "6S-6P"   = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,-1,0 ,0 ,1 ,0 ),
#                                     "11C-11S" = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,0 ,-1),
#                                     "11C-11P" = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,-1,0 ,0 ,0 ),
#                                     "11S-11P" = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,0 ,-1,0 ,0 ,1 ),
#                                     "2C-6C"   = c(0 ,1,-1, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
#                                     "2C-11C"  = c(0 ,1 ,0,-1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
#                                     "6C-11C"  = c(0, 0, 1,-1,0,0,0,0,0,0,0,0),
#                                     "2S-6S"   = c(0, 1,-1, 0,0,0,0,0,0,1,-1,0),
#                                     "2S-11S"  = c(0, 1, 0,-1,0,0,0,0,0,1,0,-1),
#                                     "6S-11S"  = c(0, 0, 1,-1,0,0,0,0,0,0,1,-1),
#                                     "2P-6P"   = c(0, 1,-1, 0,0,0,1,-1,0,0,0,0),
#                                     "2P-11P"  = c(0, 1, 0,-1,0,0,1,0,-1,0,0,0),
#                                     "6P-11P"  = c(0, 0, 1,-1,0,0,0,1,-1,0,0,0))
# 
# summary(glht(model_focal_NB, linfct=matrix_contrasts), test=adjusted("BH"))
# 
# 
# 
# #negative binomial nestmate data
# 
# model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(duration)) + (1|Replicate), data=twitch_table[which(twitch_table$ant_status=="NESTMATE"),])



