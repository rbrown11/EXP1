rm(list=ls())

#setwd("/home/bzniks/Dropbox/SeniorLectureship_Bristol/Students_postdocs/PhD_students/2021_Rachael_Brown/experiment1/with_added_replicates")
#setwd("C:/Users/racha/OneDrive/Documents/PhD/Winter 2022/Experiment/DATA ANALYSIS")
#setwd("C:/Users/racha/Downloads")
setwd("/home/rb17990/Documents/Data analysis and exploration/R scripts EXP1")


events <- read.csv(file="events.csv", header=T, sep=",")

# install.packages("pscl")
# install.packages("performance")
# install.packages("ggbeeswarm")
# install.packages("R2admb")
# install.packages("glmmADMB",
#            	repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                    	getOption("repos")),
#            	type="source")
# install.packages("glmmTMB")
install.packages("carData")
install.packages("usethis")
install.packages("devtools")

library(curl)
library(usethis)
library(devtools)
library(performance)
library(stringr)
library(lme4)
library(carData)
library(car)
library(pscl)
#library(glmmADMB)
library(ggplot2)
library(multcomp)
library(glmmTMB)
library(ggbeeswarm)
library(numDeriv)
library(e1071)

behavior_list <- c("T","S","N","G","A","B")
names(behavior_list) <- c("twitch","alarm","escape_or_death","self_grooming","allo_grooming","brood_tending")


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
events$Trial <- as.character(events$Trial)
events$group_size <- 1+as.numeric(sapply(strsplit(events$Trial, ",_"), getElement, 1))
events$treatment <- sapply(strsplit(events$Trial, ",_"), getElement, 2)
events$ant_ID <- events$Subject
events <- within (events, Petri_dish <- paste(Replicate,group_size,treatment,sep="_"))

# replicate time
events$time_of_year <- "WINTER"
events[which(grepl("\\/08\\/",events$Observation_date)),"time_of_year"] <- "SUMMER"

no_events <- events[ which(   events$ant_ID  == "NO_EVENTS"  	)  ,   ]
discards  <- events[ which(   events$ant_ID  == "DISCARD"  	)  ,   ]
events	<- events[ which( !  events$ant_ID   %in% c("DISCARD","NO_EVENTS")  	)  ,   ]



###get non-discarded replicate list
replicate_list <- rbind(unique(events[c("Replicate","treatment","group_size","Petri_dish","time_of_year")]),unique(no_events[c("Replicate","treatment","group_size","Petri_dish","time_of_year")]))

#ant status
events$ant_status <- "NESTMATE"
events[ which (events$ant_ID=="ANT_1")  , "ant_status"	] <- "FOCAL"
# unique(events$ant_status)



#replace behaviour names with more explicit behaviour names
events$Behavior_explicit <- names(behavior_list)[match(events$Behavior,behavior_list)]


##adding twitch modifiers into separate columns
#twitch size, repetition, proximity to brood and nestmates
events$Behavior <- as.character(events$Behavior)
events$Modifiers <- as.character(events$Modifiers)
events[which(events$Behavior == "T"),"twitch_size"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 1)
events[which(events$Behavior == "T"),"twitch_size"] <- sapply(strsplit(events[which(events$Behavior == "T"),"twitch_size"], "_"), getElement, 1)
events [ which(events$twitch_size == "Small") , "twitch_size" ] <- "small"
unique(events$twitch_size)

events[which(events$Behavior == "T"),"twitch_repetition"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 2)
events[ which(events$twitch_repetition == "single_twitch")  , "twitch_repetition" ] <- "single"
events[ which(events$twitch_repetition == "part_of_a_repetitive_set")  , "twitch_repetition" ] <- "repetitive"
unique(events$twitch_repetition)

events[which(events$Behavior == "T"),"proximity_to_nestmates"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 3)
events[ which(events$proximity_to_nestmates == "isolated_twitch/alone") , "proximity_to_nestmates" ] <- "alone"
events[ which(events$proximity_to_nestmates == "in_contact_with_nestmate") , "proximity_to_nestmates" ] <- "contact"
events[ which(events$proximity_to_nestmates == "in_contact_with_multiple_nestmates") , "proximity_to_nestmates" ] <- "multiple_contact"
unique(events$proximity_to_nestmates)

events[which(events$Behavior == "T"),"proximity_to_brood"] <- sapply(strsplit(events[which(events$Behavior == "T"),"Modifiers"], ","), getElement, 4 )
unique(events$proximity_to_brood)

##adding grooming modifiers into grooming
events[which(events$Behavior == "A"),"Behavior_explicit"] <- events[which(events$Behavior == "A"),"Modifiers"]
unique(events$Behavior_explicit)

#adding column to count each event
events$N <- 1

###rename event duration column
names(events)[which(names(events)=="Duration_.s.")] <- "event_duration"

#each ant gets a bloc of all possible combinations of twitch modifiers
ant_bloc_twitches <- expand.grid (  
  twitch_size = na.omit(unique(events$twitch_size))
  ,
  twitch_repetition =  na.omit(unique (events$twitch_repetition))
  ,
  proximity_to_nestmates =  na.omit(unique (events$proximity_to_nestmates))
  ,
  proximity_to_brood =  na.omit(unique (events$proximity_to_brood))
  ,
  recording_duration = 20*60
  
)
ant_bloc_behaviour <- expand.grid (  
  behaviour = na.omit(unique(events$Behavior_explicit))
  ,
  recording_duration = 20*60
)

###create a twitch table counting the number of twitches of each type for each ant
###and a behaviour table counting the number of behaviours of each type for each ant
twitch_table	<- NULL
behaviour_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list[ 1:   replicate_list[repl,"group_size"] ])
  repl_twitches <- expand.grid.df(   data.frame( replicate_list[repl,] )
                                     ,
                                     ants
                                     ,
                                     ant_bloc_twitches
  )
  repl_behaviours <- expand.grid.df(   data.frame( replicate_list[repl,] )
                                       ,
                                       ants
                                       ,
                                       ant_bloc_behaviour
  )
  
  twitch_table	<- rbind(twitch_table , repl_twitches)
  behaviour_table <- rbind(behaviour_table , repl_behaviours)
}

###count number of twitches of each type and fill twitch_table
events_Nb <- events[which(events$Behavior=="T"),c(names(twitch_table)[which(names(twitch_table)%in%names(events))],"ant_status","Behavior", "N")]
events_Nb <- aggregate(N ~ . , FUN=sum, data= events_Nb)
twitch_table <- merge(twitch_table,events_Nb,all.x=T)
###fill in missing cases
twitch_table  [ is.na(twitch_table$N)  ,"N"] <- 0
twitch_table[ which(twitch_table$ant_ID=="ANT_1")   ,	"ant_status"] <- "FOCAL"
twitch_table[ which(twitch_table$ant_ID!="ANT_1")   ,	"ant_status"] <- "NESTMATE"
###add column for whether behaviour was performed or not
twitch_table$performs <- as.numeric(twitch_table$N>0)

###count number of behaviours of each type and fill behaviour_table
events_Nb <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"ant_status","Behavior_explicit", "N")]
names(events_Nb)[which(names(events_Nb)=="Behavior_explicit")] <- "behaviour"
events_Nb <- aggregate(na.rm=T,na.action="na.pass",N ~ . , FUN=sum, data= events_Nb)

###do the same for event duration
events_duration <- events[c(names(behaviour_table)[which(names(behaviour_table)%in%names(events))],"Behavior_explicit", "event_duration")]
names(events_duration)[which(names(events_duration)=="Behavior_explicit")] <- "behaviour"
events_duration$event_duration <- as.numeric(events_duration$event_duration)
events_duration <- aggregate(na.rm=T,na.action="na.pass",event_duration ~ . , FUN=sum, data= events_duration)
###merge into behaviour table
behaviour_table <- merge(behaviour_table,events_Nb,all.x=T)
behaviour_table <- merge(behaviour_table,events_duration,all.x=T)
###fill in missing cases
behaviour_table  [ which(is.na(behaviour_table$N) ) ,"N"] <- 0
behaviour_table  [ which(is.na(behaviour_table$event_duration)&behaviour_table$behaviour!="twitch" ) ,"event_duration"] <- 0
behaviour_table  [ which(behaviour_table$behaviour=="twitch" ) ,"event_duration"] <- NA
###complement ant status column
behaviour_table[ which(behaviour_table$ant_ID=="ANT_1")   ,	"ant_status"] <- "FOCAL"
behaviour_table[ which(behaviour_table$ant_ID!="ANT_1")   ,	"ant_status"] <- "NESTMATE"
###add column for whether behaviour was performed or not
behaviour_table$performs <- as.numeric(behaviour_table$N>0)

######add a full ant ID column
twitch_table$full_ant_ID	<- with (twitch_table,paste(Replicate,treatment,group_size,ant_ID,sep="_"))
behaviour_table$full_ant_ID <- with (behaviour_table,paste(Replicate,treatment,group_size,ant_ID,sep="_"))

##deaths
deaths <- subset(events, Comment_start == "DEAD")
dead_ants <- with (deaths, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
###remove dead ants from twitch_table and behaviour_table
twitch_table <- twitch_table [ which( ! twitch_table$full_ant_ID %in% dead_ants) ,   ]
behaviour_table <- behaviour_table [ which( ! behaviour_table$full_ant_ID %in% dead_ants) ,   ]
###correct the group size
actual_group_size <- aggregate (ant_ID ~ Replicate + treatment + group_size, function(x)length(unique(x)),data=twitch_table)
names(actual_group_size)[which(names(actual_group_size)=="ant_ID")] <- "actual_group_size"
twitch_table <- merge(twitch_table,actual_group_size,all.x=T)
behaviour_table <- merge(behaviour_table,actual_group_size,all.x=T)

###escapes
escapes <- subset(events, Comment_start == "ESCAPE")
escapes$full_ant_ID <-  with (escapes, paste(Replicate,treatment,group_size,ant_ID,sep="_"))
twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] <- twitch_table[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] - as.numeric(escapes[match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   twitch_table$full_ant_ID,escapes$full_ant_ID  )))],"event_duration"])
behaviour_table[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] <- behaviour_table[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  ))), "recording_duration"] - as.numeric(escapes[match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  )[which(!is.na(match   (   behaviour_table$full_ant_ID,escapes$full_ant_ID  )))],"event_duration"])

##################################
####ANALYSIS OF BEHAVIOURS #######
##################################

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
statistical_analysis <- function (form, dataset,fam,name){
  print(paste("ANALYSIS: ",name,sep=""))
  if (fam!="normal"){
    model <-  glmer ( form, family=fam, data=dataset)
  }else{
    model <-  lmer ( form, data=dataset)
  }
  
  ###testing residuals (do only if family is normal)
  #if (fam=="normal"){
  p <- test_norm(residuals(model))
  #}
  print(Anova(model))
  ###if interaction is significant: create dummy variable
  if (Anova(model)["treatment:group_size","Pr(>Chisq)"]<0.05){
    dataset$dummy <- interaction(dataset$treatment,dataset$group_size)
    form_dummy <- update.formula(form, ~ . -treatment:group_size - treatment - group_size +dummy )
    if (fam!="normal"){
      model_dummy <-  glmer ( form_dummy, family=fam, data=dataset)
    }else{
      model_dummy <-  lmer ( form_dummy, data=dataset)
    }
    post_hocs <- summary(glht(model_dummy,linfct = mcp(dummy="Tukey")),test=adjusted("BH"))
    significance <- cld(post_hocs)
    print(significance)
    
    if (Anova(model)["treatment","Pr(>Chisq)"]<0.05){
      post_hocs_treatment	<- summary(glht(model,linfct = mcp(treatment="Tukey")),test=adjusted("BH"))
      significance_treatment <- cld(post_hocs_treatment)
      print(post_hocs_treatment)
      print(significance_treatment)
    }
    
  }else{
    form_minus_interaction <- update.formula(form, ~ . -treatment:group_size  )
    if (fam!="normal"){
      model_minus_interaction <-  glmer ( form_minus_interaction, family=fam, data=dataset)
    }else{
      model_minus_interaction <-  lmer ( form_minus_interaction, data=dataset)
    }
    if (Anova(model_minus_interaction)["treatment","Pr(>Chisq)"]<0.05){
      post_hocs_treatment	<- summary(glht(model_minus_interaction,linfct = mcp(treatment="Tukey")),test=adjusted("BH"))
      significance_treatment <- cld(post_hocs_treatment)
      print(post_hocs_treatment)
      print(significance_treatment)
    }
    
    
  }
  
}



# for (behav in unique(behaviour_table$behaviour)[which( unique(behaviour_table$behaviour)!="escape_or_death")]){
for (behav in c("twitch","self_grooming","grooming_focal_ant", "brood_tending")){
  
  ###1. subset behaviour table to only have the behaviour of interest, and split between focals and nestmates. For the focals, only use winter replicates
  dataset_focals	<- behaviour_table[which(behaviour_table$behaviour==behav&behaviour_table$ant_status=="FOCAL"),]
  dataset_nestmates <- behaviour_table[which(behaviour_table$behaviour==behav&behaviour_table$ant_status=="NESTMATE"),]
  
  ###2. Analysis for focals
  if (behav !="grooming_focal_ant"){
    #######a. Probability of performing the behaviour
    # statistical_analysis(   form = as.formula("performs ~ treatment*group_size + (1|Replicate) + (1|time_of_year)")
    #                     	,
    #                     	dataset=dataset_focals
    #                     	,
    #                     	fam="binomial"
    #                     	,
    #                     	name=paste("Probability of focals performing ",behav)
    #                     	)
    #######b. Number of occurrences across all focal ants
    statistical_analysis(   form = as.formula("N ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|time_of_year)")
                            ,
                            dataset=dataset_focals
                            ,
                            fam="poisson"
                            ,
                            name=paste("Number of ",behav," per focal ant - all focals")
    )
    
    #######c. Number of occurrences across focal ants that do peform the behaviour
    # statistical_analysis(   form = as.formula("N ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|time_of_year)")
    #                     	,
    #                     	dataset=dataset_focals[which(dataset_focals$performs==1),]
    #                     	,
    #                     	fam="poisson"
    #                     	,
    #                     	name=paste("Number of ",behav," per focal ant - only those performing the behaviour")
    # )
    if (behav!="twitch"){
      # #######d. Total duration of behaviours across all focals
      statistical_analysis(   form = as.formula("event_duration ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|time_of_year)")
                              ,
                              dataset=dataset_focals
                              ,
                              fam="normal"
                              ,
                              name=paste("Duration of ",behav," per focal ant - all focals")
      )
      #######e. Total duration of behaviours focals that perform the behaviour
      statistical_analysis(   form = as.formula("event_duration ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|time_of_year)")
                              ,
                              dataset=dataset_focals[which(dataset_focals$performs==1),]
                              ,
                              fam="normal"
                              ,
                              name=paste("Duration of ",behav," per focal ant - only focals performing the behaviour")
      )
      
      
    }
    
  }
  
  
  ###2. Analysis for nestmates
  #######a. Probability of performing the behaviour
  # statistical_analysis(   form = as.formula("performs ~ treatment*group_size + (1|Replicate) + (1|Petri_dish) + (1|time_of_year)")
  #                     	,
  #                     	dataset=dataset_nestmates
  #                     	,
  #                     	fam="binomial"
  #                     	,
  #                     	name=paste("Probability of nestmates performing ",behav)
  # )
  #######b. Number of occurrences across all nestmate ants
  
  statistical_analysis(   form = as.formula("log10(N+1/sqrt(2)) ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) + (1|time_of_year)")
                          ,
                          dataset=dataset_nestmates
                          ,
                          fam="normal"
                          ,
                          name=paste("Number of ",behav," per nestmate ant - all nestmates")
  )
  
  # #######c. Number of occurrences across focal ants that do peform the behaviour
  # statistical_analysis(   form = as.formula("log10(N+1/sqrt(2))  ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) + (1|time_of_year)")
  #                     	,
  #                     	dataset=dataset_nestmates[which(dataset_nestmates$performs==1),]
  #                     	,
  #                     	fam="normal"
  #                     	,
  #                     	name=paste("Number of ",behav," per nestmate ant - only those performing the behaviour")
  # )
  if (behav!="twitch"){
    # #######d. Total duration of behaviours across all focals
    statistical_analysis(   form = as.formula("event_duration ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) + (1|time_of_year)")
                            ,
                            dataset=dataset_nestmates
                            ,
                            fam="normal"
                            ,
                            name=paste("Duration of ",behav," per nestmate ant - all nestmates")
    )
    # #######e. Total duration of behaviours across nestmates that perform the behaviour
    statistical_analysis(   form = as.formula("event_duration ~ treatment*group_size + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish) + (1|time_of_year)")
                            ,
                            dataset=dataset_nestmates[which(dataset_nestmates$performs==1),]
                            ,
                            fam="normal"
                            ,
                            name=paste("Duration of ",behav," per nestmate ant - only nestmates performing the behaviour")
    )
    
    
  }
  
}










#
#
#
#
# ##add a replicate set
# ##to separate replicates from spring experiment and summer experiment
# total_twitches$Replicate <- as.numeric(total_twitches$Replicate)
# total_twitches$rep_set <- "1"
# total_twitches [ which (total_twitches$Replicate > 16 ) , "rep_set" ] <- "2"
#
#
# ###total number of twitches
# # total_twitches <-aggregate (cbind(frequency,N) ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
# total_twitches <-aggregate (N ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
# total_twitches$frequency <- total_twitches$N/total_twitches$recording_duration
#
#
# ##medians
# median_twitches <- aggregate ( N ~ treatment + group_size + ant_status, FUN=median,
#                            	data = total_twitches[ which(!names(total_twitches)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
#
#
#
# #####GRAPHS#####
#
#
# ##focal vs nestmate mean twitches
# standard_error <- function(x) sd(x) / sqrt(length(x))
# mean_status <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=mean, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# sd_status <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=sd, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# error_status <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=standard_error, data=total_twitches[ which(!names(total_twitches)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
#
# mean_status$se <- error_status$frequency
# colnames(mean_status) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
# mean_status$Group_size <- factor(mean_status$Group_size,
#                              	levels=c("1", "2", "6", "11"))
#
# ggplot(mean_status, aes(fill=Treatment, y=Mean, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label = cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")
#
#
# ##boxplot of frequency of twitches
# boxplot1 <- total_twitches
#
# boxplot1$group_size <- factor(boxplot1$group_size,
#                           	levels=c("1", "2", "6", "11"))
#
# ggplot(boxplot1, aes(fill=treatment, y=N, x=group_size, )) +
#   geom_boxplot(alpha=0.3) +
#   labs ( x = "Group size", y = "Frequency of body shakes" , fill = "Treatment" ) +
#   scale_fill_discrete(labels = c("Control", "Pathogen", "Sham"))+
#   facet_wrap(~ant_status)
#
# #frequency polygon of twitch decay
# scatter <- events1
# scatter$N <- 1
# scatter1 <- cbind(as.numeric(scatter$N), as.numeric(scatter$Start_.s.), scatter$group_size, scatter$treatment, scatter$ant_status)
# scatter1 <- data.frame(scatter1)
# colnames(scatter1) <- c("N", "Time", "Group_size", "Treatment", "Ant_status")
#
# scatter1$N <- as.numeric(scatter1$N)
# scatter1$Time <- as.numeric(scatter1$Time)
#
#
# scatter1$Group_size <- factor(scatter1$Group_size,
#                           	levels = c("1", "2", "6", "11"))
#
# scatter1_focal <- subset(scatter1, Ant_status=="FOCAL")
# scatter1_nestmate <- subset(scatter1, Ant_status=="NESTMATE")
#
# ggplot(scatter1, aes(Time, colour=Treatment)) +
#   geom_freqpoly(binwidth=60, size=1) +
#   xlab("Time (s)") + ylab("Body shake count") +
#   xlim(0, 1200) + facet_wrap( ~ Group_size)
#
#
# ggplot(scatter1_focal, aes(Time, colour=Treatment))+
#   geom_freqpoly(binwidth=60, size=1)+
#   xlab("Time (s)") + ylab("Body shake count") + labs(subtitle = "Focal twitches") +
#   xlim(0,1200) + ylim(0,80) + facet_wrap ( ~ Group_size)
#
#
# ggplot(scatter1_nestmate, aes(Time, colour=Treatment))+
#   geom_freqpoly(binwidth=60, size=1)+
#   xlab("Time (s)") + ylab("Body shake count") + labs(subtitle="Nestmate twitches")+
#   xlim(0,1200) + ylim(0,80) + facet_wrap ( ~ Group_size)
#
#
# #barplot mean twitches
# #not separated into focal and nestmates
# mean_twitch_per_treatment_group <- aggregate ( frequency ~ group_size + treatment, FUN=mean,  data = total_twitches[ which(!names(total_twitches)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
# standard <- aggregate ( frequency ~ group_size + treatment, FUN=sd, data = total_twitches[ which(!names(total_twitches)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ] )
# mean_twitch_per_treatment_group$sd <- standard$frequency
# error <- aggregate(frequency ~ group_size + treatment, FUN=standard_error, data = total_twitches[ which(!names(total_twitches)%in%c("twitch_size","twitch_repetition","proximity_to_nestmates","proximity_to_brood"))   ])
# mean_twitch_per_treatment_group$se <- error$frequency
# colnames(mean_twitch_per_treatment_group) <- c("group_size", "treatment", "mean", "sd", "se")
#
# ##plotting mean twitches per group per treatment
# mean_twitch_per_treatment_group$group_size <- as.character(mean_twitch_per_treatment_group$group_size)
#
# data1 <- mean_twitch_per_treatment_group
# data1$group_size <- factor(data1$group_size,
#                        	levels=c("1", "2", "6", "11"))
# ggplot(data1, aes(fill=treatment, y=mean, x=group_size, )) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
#   labs ( x = "Group size", y = "Mean body shakes per ant per minute" , fill = "Treatment" ) +
#   scale_fill_discrete(labels = c("Control", "Pathogen", "Sham"))+
#   facet_wrap(~ant_status)
#
#
# ### number of twitches of each type
# twitch_table$twitch_type <- with ( twitch_table, paste (twitch_size, twitch_repetition, sep="_"))
# twitches_per_type <-aggregate (N ~ . , FUN=sum, data = twitch_table[ which(!names(twitch_table)%in%c("proximity_to_nestmates","proximity_to_brood"))   ])
# twitches_per_type$frequency <- twitches_per_type$N/twitches_per_type$recording_duration
#
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
#                            	levels=c("1", "2", "6", "11"))
#
# ggplot(mean_type, aes(fill=treatment, y=mean, x=group_size))+
#   geom_bar(position="dodge", stat="identity")+
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
#   facet_wrap(~twitch_type)
#
#
# ##just repetitive
# mean_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition, FUN = mean, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )
# standard_dev_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition, FUN = sd, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )
# error_rep <- aggregate ( frequency ~ group_size + treatment + twitch_repetition, FUN = standard_error, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size"))   ] )
#
# mean_rep$se <- error_rep$frequency
# colnames(mean_rep) <- c("group_size", "treatment", "repetition", "mean", "standard_error")
# mean_rep$group_size <- factor(mean_rep$group_size,
#                            	levels=c("1", "2", "6", "11"))
#
# ggplot(mean_rep, aes(fill=treatment, y=mean, x=group_size))+
#   geom_bar(position="dodge", stat="identity")+
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
#   facet_wrap(~repetition)
#
# ##just size
# mean_size <- aggregate (frequency ~ group_size + treatment + twitch_size, FUN = mean, data = twitches_per_type[ which(!names(twitches_per_type)%in%c("proximity_to_nestmates", "proximity_to_brood", "twitch_repetition") ) ] )
# error_size <- aggregate ( frequency ~ group_size + treatment + twitch_size, FUN=standard_error, data = twitches_per_type [ which (!names(twitches_per_type)%in%c("proximity_to_nestmates", "proximity_to_brood", "twitch_repetition"))] )
#
# mean_size$se <- error_size$frequency
# colnames(mean_size) <- c("group_size", "treatment", "size", "mean", "standard_error")
# mean_size$group_size <- factor(mean_size$group_size,
#                            	levels=c("1", "2", "6", "11"))
#
# ggplot(mean_size, aes(fill=treatment, y=mean, x=group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9)) +
#   facet_wrap(~size)
#
#
#
#
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
#
#
# ##focal vs nestmate size and repetition
# #focal
# focal_mean_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = mean, data = twitches_per_type[ which(twitches_per_type$ant_status=="FOCAL"),  ] )
# focal_error_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = standard_error, data =  twitches_per_type[ which(twitches_per_type$ant_status=="FOCAL"),  ]  )
#
# focal_mean_type$se <- focal_error_type$frequency
# colnames(focal_mean_type) <- c("group_size", "treatment", "twitch_type", "mean", "standard_error")
# focal_mean_type$group_size <- as.character(focal_mean_type$group_size)
# focal_mean_type$group_size <- factor(focal_mean_type$group_size,
#                            	levels=c("1", "2", "6", "11"))
#
# ggplot(focal_mean_type, aes(fill=treatment, y=mean, x=group_size))+
#   geom_bar(position="dodge", stat="identity")+
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
#   facet_wrap(~twitch_type)
#
#
#
# #nestmate
# nestmate_mean_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = mean, data =  twitches_per_type[ which(twitches_per_type$ant_status=="NESTMATE"),  ] )
# nestmate_error_type <- aggregate ( frequency ~ group_size + treatment + twitch_type, FUN = standard_error, data =  twitches_per_type[ which(twitches_per_type$ant_status=="NESTMATE"),  ] )
#
# nestmate_mean_type$se <- nestmate_error_type$frequency
# colnames(nestmate_mean_type) <- c("group_size", "treatment", "twitch_type", "mean", "standard_error")
# nestmate_mean_type$group_size <- as.character(nestmate_mean_type$group_size)
# nestmate_mean_type$group_size <- factor(nestmate_mean_type$group_size,
#                                  	levels=c("1", "2", "6", "11"))
#
# ggplot(nestmate_mean_type, aes(fill=treatment, y=mean, x=group_size))+
#   geom_bar(position="dodge", stat="identity")+
#   geom_errorbar(aes(ymin=mean-standard_error, ymax=mean+standard_error), width=.2, position=position_dodge(0.9))+
#   facet_wrap(~twitch_type)
#
#
#
# ##new replicates barplot
# ##checking data from summer experiment
# new_replicates <- subset(total_twitches, total_twitches$rep_set == "2")
#
# mean_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=mean, data=new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# sd_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=sd, data=new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
# error_status_new <- aggregate( frequency ~ group_size + treatment + ant_status, FUN=standard_error, data=new_replicates[ which(!names(new_replicates)%in%c("proximity_to_nestmates","proximity_to_brood", "twitch_size", "twitch_repetition"))   ])
#
# mean_status_new$se <- error_status_new$frequency
# colnames(mean_status_new) <- c("Group_size", "Treatment", "ant_status", "Mean", "standard_error")
# mean_status_new$Group_size <- factor(mean_status_new$Group_size,
#                              	levels=c("1", "2", "6", "11"))
#
# ggplot(mean_status_new, aes(fill=Treatment, y=Mean, x=Group_size)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_errorbar(aes(ymin=Mean-standard_error, ymax=Mean+standard_error), width=.2, position=position_dodge(0.9)) +
#   #geom_text(aes(label = cld, y = Mean + standard_error), vjust= -0.5, position = position_dodge(0.9)) +
#   facet_wrap(~ant_status) +
#   labs(x="Group size", y="Number of body shakes per ant per minute")
#
#
#
#
#
#
#
# #####STATS#####
#
#
# head(total_twitches)
# hist(total_twitches$frequency)
# 100*sum(total_twitches$frequency == 0)/nrow(total_twitches) #69% of data is zeros
# str(total_twitches)
#
# ##re-structuring data
# total_twitches$interac_size_treat <- with (total_twitches, interac_size_treat <- paste(group_size,treatment,sep="_"))
# total_twitches$Replicate <- factor(total_twitches$Replicate)
# total_twitches$group_size <- factor(total_twitches$group_size)
# total_twitches$treatment <- factor(total_twitches$treatment)
# total_twitches$Petri_dish <- with (total_twitches, interac_size_treat <- paste(Replicate,group_size,treatment,sep="_"))
#
# matrix_contrasts <- rbind( "1C-1S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 , 0),
#                        	"1C-1P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"1S-1P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"2C-2S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 ,-1,  0 , 0),
#                        	"2C-2P"   = c(0 , 0 , 0 , 0 ,-1,  0,-1,  0 , 0 , 0 , 0 , 0),
#                        	"2S-2P"   = c(0 , 0 , 0 , 0 ,-1,  1,-1,  0 , 0 , 1 , 0 , 0),
#                        	"6C-6S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 ,-1,  0),
#                        	"6C-6P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 ,-1,  0 , 0 , 0 , 0),
#                        	"6S-6P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 ,-1,  0 , 0 , 1 , 0),
#                        	"11C-11S" = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 ,-1),
#                        	"11C-11P" = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 ,-1,  0 , 0 , 0),
#                        	"11S-11P" = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 ,-1,  0 , 0 , 1),
#                        	"1C-2C"   = c(0 ,-1,  0 , 0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"1C-6C"   = c(0 , 0 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"1C-11C"  = c(0 , 0 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"2C-6C"   = c(0 , 1 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"2C-11C"  = c(0 , 1 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#                        	"6C-11C"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  0,  0),
#                        	"1S-2S"   = c(0, -1,  0,  0,  0,  0, 0,  0,  0, -1,  0  ,0),
#                        	"1S-6S"   = c(0,  0, -1,  0,  0,  0, 0,  0,  0,  0, -1,  0),
#                        	"1S-11S"  = c(0,  0,  0, -1,  0,  0, 0,  0,  0,  0,  0, -1),
#                        	"2S-6S"   = c(0,  1, -1,  0,  0,  0, 0,  0,  0,  1, -1,  0),
#                        	"2S-11S"  = c(0,  1,  0, -1,  0,  0, 0,  0,  0,  1,  0, -1),
#                        	"6S-11S"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  1, -1),
#                        	"1P-2P"   = c(0, -1,  0,  0,  0,  0,-1,  0,  0,  0,  0,  0),
#                        	"1P-6P"   = c(0,  0, -1,  0,  0,  0, 0, -1,  0,  0,  0,  0),
#                        	"1P-11P"  = c(0,  0,  0, -1,  0,  0, 0,  0, -1,  0,  0,  0),
#                        	"2P-6P"   = c(0,  1, -1,  0,  0,  0, 1, -1,  0,  0,  0,  0),
#                        	"2P-11P"  = c(0,  1,  0, -1,  0,  0, 1,  0, -1,  0,  0,  0),
#                        	"6P-11P"  = c(0,  0,  1, -1,  0,  0, 0,  1, -1,  0,  0,  0))
#
# ###Poisson model
# ###Poisson model is both overdispersed and zero-inflated
# model_focal <- glmer(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate),
#                  	family="poisson", data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# check_zeroinflation(model_focal)
# check_overdispersion(model_focal)
# Anova(model_focal)
#
# #interaction significant
# model_focal_interac <- glmer(N ~ interac_size_treat + offset(log(recording_duration)) + (1|rep_set:Replicate), family="poisson",
#                          	data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])
# post_hoc_focal_interac_tukey<- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# post_hoc_focal_interac <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))
#
# print(cld(post_hoc_focal_interac_tukey))
# #not matching graphs
#
# ###nestmate poisson
# model_nestmate <- glmer(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate) + (1|Petri_dish),
#                     	family="poisson", data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# check_zeroinflation(model_nestmate)
# check_overdispersion(model_nestmate)
# Anova(model_nestmate)
#
# #interaction non significant
# model_nestmate_ni <- glmer(N ~ group_size + treatment + offset(log(recording_duration)) + (1|rep_set:Replicate) + (1|Petri_dish),
#                        	family="poisson", data=total_twitches[which(total_twitches$ant_status=="NESTMATE") ,])
# Anova(model_nestmate_ni)
#
# post_hoc_nestmate_treat<- summary(glht(model_nestmate_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_treat))
#
# post_hoc_nestmate_group<- summary(glht(model_nestmate_ni, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# print(cld(post_hoc_nestmate_group))
#
#
#
#
# ##negative binom model
# ##focal
# #model not zero inflated or overdispersed
# model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate), data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# check_zeroinflation(model_focal_NB)
# check_overdispersion(model_focal_NB)
# Anova(model_focal_NB)
#
# #significant interaction
# #post hocs don't look too bad, check against graphs
# model_focal_interac <- glmer.nb(N ~ interac_size_treat + offset(log(recording_duration)) + (1|rep_set:Replicate),
#                             	data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])
# post_hoc_focal_interac_nb <- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# post_hoc_focal_interac_nb_matrix <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))
#
#
#
# ##nestmate
# #not zero inflated or overdispersed
# ##output looks v wrong
# model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# check_zeroinflation(model_nestmate_NB)
# check_overdispersion(model_nestmate_NB)
# Anova(model_nestmate_NB)
#
#
# ##trying a model check I saw online, not sure if it is good or not!
# #add more iterations
# #output looks more reliable
# ss <- getME(model_nestmate_NB,c("theta", "fixef"))
# m2 <- update(model_nestmate_NB, start=ss, control=glmerControl(optCtrl = list(maxfun=2e4)))
# Anova(m2)
#
# #interaction non sig
# #don't trust result
# model_nestmate_NB_ni <- glmer.nb(N ~ group_size+treatment + offset(log(recording_duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# Anova(model_nestmate_NB_ni)
#
# post_hoc_nestmate_treatment <- summary(glht(model_nestmate_NB_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# post_hoc_nestmate_group <- summary(glht(model_nestmate_NB_ni, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
#
#
#
# # dd <- getME(model_nestmate_NB_ni, c("theta", "fixef"))
# # m3 <- update(model_nestmate_NB_ni, start=dd, control=glmerControl(optCtrl = list(maxfun=2e4)))
# #
# # Anova(m3)
# # post_hoc_nestmate_treatment <- summary(glht(m3, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
#
#
#
#
#
#
# ####DRAFTS####
#
# # ##Hurdle model
# # ###not working
# # sum(total_twitches$N <1)
# #
# # mod.hurdle <- hurdle(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate),
# #                  	data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# #
# # hurdle_attempt_focal <- glmmTMB(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate),
# #                             	data=total_twitches[which(total_twitches$ant_status=="FOCAL"),], zi=T, family=truncated_poisson)
# ##on new rep set only
# ##new replicates from summer experiment
# #
# #
# # model_focal <- glmer(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate),
# #                  	family="poisson", data=new_replicates[which(new_replicates$ant_status=="FOCAL"),])
# # check_zeroinflation(model_focal)
# # check_overdispersion(model_focal)
# # Anova(model_focal)
# #
# # model_focal_interac <- glmer(N ~ interac_size_treat + offset(log(recording_duration)) + (1|rep_set:Replicate), family="poisson",
# #                          	data=new_replicates[which(new_replicates$ant_status=="FOCAL") , ])
# # post_hoc_focal_interac<- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# #
# #
# # print(cld(post_hoc_focal_interac))
#
# # model_nestmate_interac <- glmer(N ~ interac_size_treat + offset(log(recording_duration)) + (1|Replicate) + (1|Petri_dish),
# #                             	family="poisson", data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# # Anova(model_nestmate_interac)
# # post_hoc_nestmate_interac <- summary(glht(model_nestmate_interac, linfct=matrix_contrasts_nestmates), test=adjusted("BH"))
#
# ##negative binom model
# # ##focal
# # #model not zero inflated or overdispersed
# # model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate), data=total_twitches[which(total_twitches$ant_status=="FOCAL"),])
# # check_zeroinflation(model_focal_NB)
# # check_overdispersion(model_focal_NB)
# # Anova(model_focal_NB)
# #
# # #significant interaction
# # #post hocs don't look too bad, check against graphs
# # ###no sig differences between S and P, but sig diff between control and experimental
# # model_focal_interac <- glmer.nb(N ~ interac_size_treat + offset(log(recording_duration)) + (1|rep_set:Replicate),
# #                          	data=total_twitches[which(total_twitches$ant_status=="FOCAL") , ])
# # post_hoc_focal_interac_nb <- summary(glht(model_focal_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# # post_hoc_focal_interac_nb_matrix <- summary(glht(model_focal_interac, linfct=matrix_contrasts), test=adjusted("BH"))
# #
# #
# # ##nestmate
# # #not zero inflated or overdispersed
# # model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# # check_zeroinflation(model_nestmate_NB)
# # check_overdispersion(model_nestmate_NB)
# # Anova(model_nestmate_NB)
# #
# # #add more iterations
# # ss <- getME(model_nestmate_NB,c("theta", "fixef"))
# # m2 <- update(model_nestmate_NB, start=ss, control=glmerControl(optCtrl = list(maxfun=2e4)))
# # Anova(m2)
# #
# # #now interaction non sig
# # model_nestmate_NB_ni <- glmer.nb(N ~ group_size+treatment + offset(log(recording_duration)) + (1|rep_set:Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# # Anova(model_nestmate_NB_ni)
# #
# # post_hoc_nestmate_treatment <- summary(glht(model_nestmate_NB_ni, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# # post_hoc_nestmate_group <- summary(glht(model_nestmate_NB_ni, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# #
# # # dd <- getME(model_nestmate_NB_ni, c("theta", "fixef"))
# # # m3 <- update(model_nestmate_NB_ni, start=dd, control=glmerControl(optCtrl = list(maxfun=2e4)))
#
#
# # #check singularity
# # #not an issue
# # tt <- getME(model_nestmate_NB,"theta")
# # ll <- getME(model_nestmate_NB, "lower")
# # min(tt[ll==0])
# #
# # #double checking gradient calculations
# # derivs1 <- model_nestmate_NB@optinfo$derivs
# # sc_grad1 <- with(derivs1, solve(Hessian, gradient))
# # max(abs(sc_grad1))
# # max(pmin(abs(sc_grad1), abs(derivs1$gradient)))
#
# #
# # #significant interaction
# # model_nestmate_NB_interac <- glmer.nb(N ~ interac_size_treat + offset(log(recording_duration)) + (1|rep_set:Replicate)+
# #                                     	(1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# # post_hoc_nestmate_interac_nb <- summary(glht(model_nestmate_NB_interac, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_nestmate_interac_nb))
#
#
# #
# # ###interaction is not significant, proceed with simple post-hocs
# # post_hoc_nestmate_treatment<- summary(glht(model_nestmate_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_nestmate_treatment))
# #
# # post_hoc_nestmate_size<- summary(glht(model_nestmate_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_nestmate_size))
#
#
#
# # ###interaction is not significant, proceed with simple post-hocs
# # post_hoc_focal_treatment<- summary(glht(model_focal_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_focal_treatment))
# #
# # post_hoc_focal_size<- summary(glht(model_focal_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_focal_size))
#
# ###very odd post-hoc outcomes where all group size come out as equal!!
# ##try an alternative
# # zinbinom_focal<- glmmTMB(N~treatment*group_size+
# #                        	offset(log(recording_duration))+(1|Replicate),
# #                      	data=total_twitches[which(total_twitches$ant_status=="FOCAL"),],
# #                      	zi=~treatment+group_size,
# #                      	family=nbinom2)
# # Anova(zinbinom_focal)
# # post_hoc_focal_treatment<- summary(glht(zinbinom_focal, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_focal_treatment))
# # post_hoc_focal_size<- summary(glht(zinbinom_focal, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_focal_size))
# ###STRANGE POST_HOCS WITH NEGATIVE BINOMIAL GLMER!
#
#
#
# ###Poisson model is both overdispersed and zero-inflated
# # model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|Replicate)+ (1|Petri_dish), data=total_twitches[which(total_twitches$ant_status=="NESTMATE"),])
# # check_zeroinflation(model_nestmate_NB)
# # check_overdispersion(model_nestmate_NB)
# # ###still overdispersed
# # Anova(model_nestmate_NB)
# #
# # ###interaction is not significant, proceed with simple post-hocs
# # post_hoc_nestmate_treatment<- summary(glht(model_nestmate_NB, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_nestmate_treatment))
# #
# # post_hoc_nestmate_size<- summary(glht(model_nestmate_NB, linfct=mcp(group_size="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_nestmate_size))
#
#
# ##zero inflation models
# # zipoisson_focal <- glmmTMB(N~treatment*group_size+
# #                          	offset(log(recording_duration))+(1|Replicate),
# #                        	data=total_twitches[which(total_twitches$ant_status=="FOCAL"),],
# #                        	ziformula=~.,
# #                        	family=poisson)
# # Anova(zipoisson_focal)
#
#
# # zipoisson_focal <-glmmadmb(frequency~treatment*group_size+(1|Replicate),zeroInflation=TRUE,data=twitch_table[which(twitch_table$ant_status=="FOCAL"),],family = "poisson")
# # Anova(zipoisson_focal)
# # check_zeroinflation(zipoisson_focal) ##overfitting
#
# # zipoisson_focal <- glmmadmb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|Replicate),
# #                         	zeroInflation = TRUE, data=twitch_table[which(twitch_table$ant_status=="FOCAL"),],family = "poisson")
#
#
# # model_focal_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|Replicate), data=twitch_table[which(twitch_table$ant_status=="FOCAL"),])
# #
# #
# # model_focal_interac_NB <- glmer.nb(N ~ interac_size_treat + offset(log(recording_duration)) + (1|Replicate), family="poisson", data=twitch_table[which(twitch_table$ant_status=="FOCAL"),])
# # Anova(model_focal_interac_NB)
# # post_hoc_focal_interac_NB <- summary(glht(model_focal_interac_NB, linfct=mcp(interac_size_treat="Tukey")), test=adjusted("BH"))
# # print(cld(post_hoc_focal_interac_NB))
# #   
#  # matrix_contrasts <- rbind( "1C-1S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"1C-1P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"1S-1P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"2C-2S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 ,-1,  0 , 0),
#  #                       	"2C-2P"   = c(0 , 0 , 0 , 0 ,-1,  0,-1,  0 , 0 , 0 , 0 , 0),
#  #                       	"2S-2P"   = c(0 , 0 , 0 , 0 ,-1,  1,-1,  0 , 0 , 1 , 0 , 0),
#  #                       	"6C-6S"   = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 ,-1,  0),
#  #                       	"6C-6P"   = c(0 , 0 , 0 , 0 ,-1,  0 ,0 ,-1,  0 , 0 , 0 , 0),
#  #                       	"6S-6P"   = c(0 , 0 , 0 , 0 ,-1,  1 ,0 ,-1,  0 , 0 , 1 , 0),
#  #                       	"11C-11S" = c(0 , 0 , 0 , 0 , 0 ,-1, 0 , 0 , 0 , 0 , 0 ,-1),
#  #                       	"11C-11P" = c(0 , 0 , 0 , 0 ,-1,  0 ,0 , 0 ,-1,  0 , 0 , 0),
#  #                       	"11S-11P" = c(0 , 0 , 0 , 0 ,-1,  1 ,0 , 0 ,-1,  0 , 0 , 1),
#  #                       	"1C-2C"   = c(0 ,-1,  0 , 0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"1C-6C"   = c(0 , 0 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"1C-11C"  = c(0 , 0 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"2C-6C"   = c(0 , 1 ,-1,  0 , 0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"2C-11C"  = c(0 , 1 , 0 ,-1,  0 , 0 ,0 , 0 , 0 , 0 , 0 , 0),
#  #                       	"6C-11C"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  0,  0),
#  #                       	"1S-2S"   = c(0, -1,  0,  0,  0,  0, 0,  0,  0, -1,  0  ,0),
#  #                       	"1S-6S"   = c(0,  0, -1,  0,  0,  0, 0,  0,  0,  0, -1,  0),
#  #                       	"1S-11S"  = c(0,  0,  0, -1,  0,  0, 0,  0,  0,  0,  0, -1),
#  #                       	"2S-6S"   = c(0,  1, -1,  0,  0,  0, 0,  0,  0,  1, -1,  0),
#  #                       	"2S-11S"  = c(0,  1,  0, -1,  0,  0, 0,  0,  0,  1,  0, -1),
#  #                       	"6S-11S"  = c(0,  0,  1, -1,  0,  0, 0,  0,  0,  0,  1, -1),
#  #                       	"1P-2P"   = c(0, -1,  0,  0,  0,  0,-1,  0,  0,  0,  0,  0),
#  #                       	"1P-6P"   = c(0,  0, -1,  0,  0,  0, 0, -1,  0,  0,  0,  0),
#  #                       	"1P-11P"  = c(0,  0,  0, -1,  0,  0, 0,  0, -1,  0,  0,  0),
#  #                       	"2P-6P"   = c(0,  1, -1,  0,  0,  0, 1, -1,  0,  0,  0,  0),
#  #                       	"2P-11P"  = c(0,  1,  0, -1,  0,  0, 1,  0, -1,  0,  0,  0),
#  #                       	"6P-11P"  = c(0,  0,  1, -1,  0,  0, 0,  1, -1,  0,  0,  0))
#
# # matrix_contrasts_nestmates <- rbind("2C-2S"   = c(0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,-1 ,0 ,0 ),
# #                                 	"2C-2P"   = c(0 ,0 ,0 ,0,-1 , 0 ,-1,0 ,0 ,0 ,0 ,0 ),
# #                                 	"2S-2P"   = c(0 ,0 ,0 ,0,-1,1 ,-1,0 ,0 ,1 ,0 ,0 ),
# #                                 	"6C-6S"   = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,-1,0 ),
# #                                 	"6C-6P"   = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,-1,0 ,0 ,0 ,0 ),
# #                                 	"6S-6P"   = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,-1,0 ,0 ,1 ,0 ),
# #                                 	"11C-11S" = c(0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,0 ,-1),
# #                                 	"11C-11P" = c(0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,-1,0 ,0 ,0 ),
# #                                 	"11S-11P" = c(0 ,0 ,0 ,0 ,-1,1 ,0 ,0 ,-1,0 ,0 ,1 ),
# #                                 	"2C-6C"   = c(0 ,1,-1, 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
# #                                 	"2C-11C"  = c(0 ,1 ,0,-1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ),
# #                                 	"6C-11C"  = c(0, 0, 1,-1,0,0,0,0,0,0,0,0),
# #                                 	"2S-6S"   = c(0, 1,-1, 0,0,0,0,0,0,1,-1,0),
# #                                 	"2S-11S"  = c(0, 1, 0,-1,0,0,0,0,0,1,0,-1),
# #                                 	"6S-11S"  = c(0, 0, 1,-1,0,0,0,0,0,0,1,-1),
# #                                 	"2P-6P"   = c(0, 1,-1, 0,0,0,1,-1,0,0,0,0),
# #                                 	"2P-11P"  = c(0, 1, 0,-1,0,0,1,0,-1,0,0,0),
# #                                 	"6P-11P"  = c(0, 0, 1,-1,0,0,0,1,-1,0,0,0))
# #
# # summary(glht(model_focal_NB, linfct=matrix_contrasts), test=adjusted("BH"))
# #
# #
# #
# # #negative binomial nestmate data
# #
# # model_nestmate_NB <- glmer.nb(N ~ group_size*treatment + offset(log(recording_duration)) + (1|Replicate), data=twitch_table[which(twitch_table$ant_status=="NESTMATE"),])
#
#
#