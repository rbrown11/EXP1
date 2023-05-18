
########################################
##### NB SHAKES PER MINUTE PER ANT #####
########################################


###### LINE GRAPH #######

events1$Start_.s.<- as.numeric(events1$Start_.s.)
unique(events1$Start_.s.)

events1$scatter <- round_any(events1$Start_.s., 60)

scat <- aggregate(N ~ scatter + ant_status + treatment + group_size, FUN=sum, data=events1)

ggplot(scat, aes(x=scatter, y=N, color=treatment))+
  geom_point(aes(fill=treatment))+
  facet_grid(treatment~ant_status) +
  labs(x="Time (seconds)", y="Number of body shakes performed")+
  # geom_smooth(method=lm, se=F, aes(fill=treatment))+
  geom_line(size=1)+
  theme_bw()+
  # ylim(0, 100)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )







##### NB SHAKES IN 5 MINUTE TIME BINS #####



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



ggplot(shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),],
       aes(x=group_size, y=N, fill=treatment))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~smins)

ggplot(shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),],
       aes(x=group_size, y=N, fill=treatment))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~smins)

## just count data, have a look at means too
## statistics for this 
## time effect into model
## N ~ group_size*treatment*smins ???


###### GRAPHS FOR MEAN BODY SHAKES PER ANT ###### 
#DURATION?

mean_shake_mins <- aggregate(frequency ~ group_size + treatment + ant_status + smins, FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~ group_size + treatment + ant_status + smins, FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c("group_size", "Treatment", "ant_status", "time_bins", "mean_shakes", "se")

ggplot(mean_shake_mins, aes(fill=Treatment, y=mean_shakes, x=group_size))+
  geom_errorbar(aes(ymin=mean_shakes-se, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(ant_status ~ time_bins)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )




##without group size

mean_shake_mins <- aggregate(frequency ~  treatment + ant_status + smins, FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~  treatment + ant_status + smins, FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c("Treatment", "ant_status", "time_bins", "mean_shakes", "se")

ggplot(mean_shake_mins, aes(fill=Treatment, y=mean_shakes, x=Treatment))+
  geom_errorbar(aes(ymin=mean_shakes-se, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(ant_status ~ time_bins)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )


## without group size and treatment

mean_shake_mins <- aggregate(frequency ~  ant_status + smins, FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~ ant_status + smins, FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c( "ant_status", "time_bins", "mean_shakes", "se")



ggplot(shake_time_table1, aes(fill=ant_status, y=N, x=smins))+
  # geom_errorbar(aes(ymin=mean_shakes-se, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  # geom_bar(position="dodge", stat="identity") +
  # geom_violin(trim=T)+
  # geom_boxplot()+
  geom_point (aes(color=ant_status))+
  stat_summary(fun.y=mean, geom="bar", color="black", size=1)+
  # geom_line(size=1, color="red")+
  facet_wrap( ~ ant_status) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )

  # ggplot(sum_status_int, aes(fill=treatment, y=N, x=group_size, color=treatment))+
#   geom_beeswarm(cex=1, priority = "density") +
#   labs(x="group size", y="Number of body shakes of each ant") +
#   facet_wrap(~ant_status)




### STATISTICS ###

shake_time_table1$interac_size_treat <- with (shake_time_table1, interac_size_treat <- paste(group_size,treatment,sep="_"))
shake_time_table1$Petri_dish <- with (shake_time_table1, interac_size_treat <- paste(Replicate,group_size,treatment,sep="_"))

shake_time_table1$Replicate <- as.numeric(shake_time_table1$Replicate)
shake_time_table1$time <- "WINTER"
shake_time_table1 [ which (shake_time_table1$Replicate > 16 ) , "time" ] <- "SUMMER"



### FOCAL MODEL ###


model_focal <- glmer(N^0.1 ~ group_size*treatment+smins + offset(log(duration)) + (1|Replicate),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),])

qqnorm(residuals(model_focal))
qqline(residuals(model_focal))
hist(residuals(model_focal))

test_norm(residuals(model_focal))

Anova(model_focal, type="III")

post_hoc_smins <- summary(glht(model_focal, linfct=mcp(smins="Tukey")), test=adjusted("BH"))

### brain dump ###

## model was first ran with all interactions, but then only group_size:treatment was significant
## none of the smins interactions were significant so i removed them
## once the model was re ran smins was significant, with post hoc scores that look like:
#                             Estimate Std. Error z value Pr(>|z|)    
# 5-10mins - 0-5mins == 0    -0.07678    0.02235  -3.436  0.00118 ** 
# 10-15mins - 0-5mins == 0   -0.09893    0.02235  -4.427 2.87e-05 ***
# 15-20mins - 0-5mins == 0   -0.11759    0.02235  -5.262 8.57e-07 ***
# 10-15mins - 5-10mins == 0  -0.02215    0.02235  -0.991  0.38594    
# 15-20mins - 5-10mins == 0  -0.04080    0.02235  -1.826  0.10184    
# 15-20mins - 10-15mins == 0 -0.01865    0.02235  -0.835  0.40397    

## which suggests that 0-5 mins contains more shakes than the rest of the time for focal ants
## regardless of treatment?
## but the interaction between group size and treatment is still there so model_focal_interac in exp code (3) is still valid
# should try to find a nice way to show distribution of data 



### NESTMATES MODEL ###


model_nestmates <- glmer(N^0.1 ~ group_size+treatment+smins + offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),])


# qqnorm(residuals(model_nestmates))
# qqline(residuals(model_nestmates))
# hist(residuals(model_nestmates))

test_norm(residuals(model_nestmates))

Anova(model_nestmates, type="III")


post_smin <- summary(glht(model_nestmates, linfct=mcp(smins="Tukey")), test=adjusted("BH"))
#                             Estimate Std. Error z value Pr(>|z|)    
# 5-10mins - 0-5mins == 0    -0.05400    0.01336  -4.042 0.000106 ***
# 10-15mins - 0-5mins == 0   -0.07393    0.01336  -5.534 9.37e-08 ***
# 15-20mins - 0-5mins == 0   -0.09401    0.01336  -7.037 1.18e-11 ***
# 10-15mins - 5-10mins == 0  -0.01993    0.01336  -1.492 0.135672    
# 15-20mins - 5-10mins == 0  -0.04001    0.01336  -2.995 0.004117 ** 
# 15-20mins - 10-15mins == 0 -0.02008    0.01336  -1.503 0.135672    

print(cld(post_smin))
# 0-5mins  5-10mins 10-15mins 15-20mins 
# "c"       "b"      "ab"       "a" 


post_treat <- summary(glht(model_nestmates, linfct=mcp(treatment="Tukey")), test=adjusted("BH"))
#             Estimate Std. Error z value Pr(>|z|)   
# P - C == 0  0.13608    0.02149   6.332 7.28e-10 ***
# S - C == 0  0.07590    0.02033   3.733 0.000283 ***
# S - P == 0 -0.06018    0.02127  -2.830 0.004657 ** 

print(cld(post_treat))
#  C   P   S 
# "a" "c" "b" 



### brain dump ###

## model ran with all random factors included for nestmates 
## interaction terms non significant so removed, and treatment and time window significant
## post hocs ran, results above, and they match the graph


