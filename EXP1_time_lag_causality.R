
#####################################################
############# TIME SERIES RELATIONSHIPS #############
########### SHAKING AND OTHER BEHAVIOURS ############
####### BOTH BINNED TIME SERIES AND RAW DATA ########
#####################################################


### time between allogrooming (or self grooming?) and first contact with twitch
## need to dissect the in contact with column
## include the start times of the events and the ant IDs
### something like 'ant 1 twitches in contact with ant 4, then ant 4 grooms with this time lag'
### not just focal twitches though - what about nestmates twitching with nestmates then allogrooming of focal?
## then look at time lagged causality stuff for p vs s vs c
#### twitch decay over time
## do pathgoen trials keep consistently shaking for longer?

events$petri_dish <- with ( events, paste (Replicate, treatment, group_size, sep="_"))


########################################
##### NB SHAKES PER MINUTE PER ANT #####
########################################

##subset events to only include shake behaviours
events2 <- events[ which(events$Behavior == "T")   ,   ]

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

##sanity checks
shake_per_min <- aggregate(N ~ group_size + treatment + ant_status + time_bins, FUN=sum, data=events2)
shakes_min <- aggregate(N ~ group_size + treatment + ant_ID + time_bins + Replicate + rep_set + Observation_date,
                        FUN=sum, data=events2)

###create a data frame with all options to add zeros in
replicate_list <- rbind(unique(events[c("Replicate","treatment","group_size")]),unique(no_events[c("Replicate","treatment","group_size")]))

ant_list <- paste("ANT",c(1:11),sep="_")

shake_bloc <- expand.grid (
  time_bins = unique (events2$time_bins)
)


shake_time_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list [ 1 : replicate_list [repl, "group_size"] ] )
  repl_shake <- expand.grid.df( data.frame ( replicate_list[repl,])
                                , ants,
                                shake_bloc)
  shake_time_table <- rbind(shake_time_table, repl_shake)
}


shake_Nb <- aggregate(N ~ Replicate + treatment + group_size + ant_ID + time_bins + rep_set, data=events2, FUN=sum  )
shake_time_table1 <- merge(shake_time_table, shake_Nb, all.x=T)


##add in missing values
shake_time_table1 [is.na (shake_time_table1$N), "N" ] <- 0
shake_time_table1$ant_status <- "NESTMATE"
shake_time_table1 [ which ( shake_time_table1$ant_ID=="ANT_1"), "ant_status"] <- "FOCAL"
shake_time_table1$rep_set <- "1"
shake_time_table1 [ which (shake_time_table1$Replicate > 16 ) , "rep_set" ] <- "2"


##create a column to identify unique ants
shake_time_table1$full_ant_ID <-  with (shake_time_table1,paste(Replicate,treatment,group_size,ant_ID,sep="_"))




##############################################
##### NB ALLOGROOMING PER MINUTE PER ANT #####
##############################################

##subset so only allogrooming behaviours
allogrooming_table <- subset(events, Behavior=="A")

##add a rep set for winter and summer
allogrooming_table$rep_set <- "1"
allogrooming_table [ which (allogrooming_table$Replicate > 16 ) , "rep_set" ] <- "2"

##remove the summer focals
events2_focals <- subset(allogrooming_table, ant_status=="FOCAL")
events2_nestmates <- subset(allogrooming_table, ant_status=="NESTMATE")
events2_focals_first_rep <- subset(events2_focals, rep_set=="1")
allogrooming_table <- rbind(events2_nestmates, events2_focals_first_rep)

##add the time bins for events
allogrooming_table$Start_.s. <- as.numeric(allogrooming_table$Start_.s.)
allogrooming_table$time_bins <- round_any(allogrooming_table$Start_.s., 60)

##include a tally value
allogrooming_table$N <- 1
allo_per_min <- aggregate(N ~ group_size + treatment + ant_status + time_bins, FUN=sum, data=allogrooming_table)


###create a data frame with all options to add zeros in
allo_bloc <- expand.grid( time_bins = unique ( allogrooming_table$time_bins )
                          # ,
                          # modifiers = unique ( allogrooming_table$Modifiers )
)

allo_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list [ 1 : replicate_list [repl, "group_size"] ] )
  repl_groom <- expand.grid.df( data.frame ( replicate_list[repl,])
                                , ants,
                                allo_bloc)
  allo_table <- rbind(allo_table, repl_groom)
}

allo_Nb <- aggregate(N ~ Replicate + treatment + group_size + ant_ID + time_bins + rep_set, data=allogrooming_table, FUN=sum  )
allo_time_table1 <- merge(allo_table, allo_Nb, all.x=T)

##add in missing values
allo_time_table1 [is.na (allo_time_table1$N), "N" ] <- 0
allo_time_table1$ant_status <- "NESTMATE"
allo_time_table1 [ which ( allo_time_table1$ant_ID=="ANT_1"), "ant_status"] <- "FOCAL"
allo_time_table1$rep_set <- "1"
allo_time_table1 [ which (allo_time_table1$Replicate > 16 ) , "rep_set" ] <- "2"

##add full ant ID 
allo_time_table1$full_ant_ID <-  with (allo_time_table1,paste(Replicate,treatment,group_size,ant_ID,sep="_"))


##combine the tally with the shake table
shake_time_table2 <- shake_time_table1 
shake_time_table2 <- shake_time_table2 %>% 
  right_join(allo_time_table1, by = c("full_ant_ID", "time_bins"))

shake_time_table2 <- select(shake_time_table2, -c("Replicate.y", "treatment.y", "group_size.y",
                                                  "ant_ID.y", "rep_set.y", "ant_status.y" ))


colnames(shake_time_table2) <- c("Replicate", "treatment", "group_size", "ant_ID", "time_bins",
                                 "rep_set", "Nb_shakes", "ant_status", "full_ant_ID", "Nb_allogrooming")




##############################################
##### NB SELFGROOMING PER MINUTE PER ANT #####
##############################################

##subset so only selfgrooming behaviours
self_grooming_table <- subset(events, Behavior=="G")

##add rep set
self_grooming_table$rep_set <- "1"
self_grooming_table [ which (self_grooming_table$Replicate > 16 ) , "rep_set" ] <- "2"

##remove summer focals
events2_focals <- subset(self_grooming_table, ant_status=="FOCAL")
events2_nestmates <- subset(self_grooming_table, ant_status=="NESTMATE")
events2_focals_first_rep <- subset(events2_focals, rep_set=="1")
self_grooming_table <- rbind(events2_nestmates, events2_focals_first_rep)


##add time bins
self_grooming_table$Start_.s. <- as.numeric(self_grooming_table$Start_.s.)
self_grooming_table$time_bins <- round_any(self_grooming_table$Start_.s., 60)

##add tally / sanity checks
self_grooming_table$N <- 1
allo_per_min <- aggregate(N ~ group_size + treatment + ant_status + time_bins, FUN=sum, data=self_grooming_table)


###create a data frame with all options to add zeros in
self_bloc <- expand.grid( time_bins = unique ( self_grooming_table$time_bins ))

self_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list [ 1 : replicate_list [repl, "group_size"] ] )
  repl_groom <- expand.grid.df( data.frame ( replicate_list[repl,])
                                , ants,
                                self_bloc)
  self_table <- rbind(self_table, repl_groom)
}




self_Nb <- aggregate(N ~ Replicate + treatment + group_size + ant_ID + time_bins + rep_set, data=self_grooming_table, FUN=sum  )
self_time_table1 <- merge(self_table, self_Nb, all.x=T)

##add missing values in
self_time_table1 [is.na (self_time_table1$N), "N" ] <- 0
self_time_table1$ant_status <- "NESTMATE"
self_time_table1 [ which ( self_time_table1$ant_ID=="ANT_1"), "ant_status"] <- "FOCAL"
self_time_table1$rep_set <- "1"
self_time_table1 [ which (self_time_table1$Replicate > 16 ) , "rep_set" ] <- "2"

##add full ant ID
self_time_table1$full_ant_ID <-  with (self_time_table1,paste(Replicate,treatment,group_size,ant_ID,sep="_"))


##combine tally to existing table
shake_time_table2 <- shake_time_table2 %>% 
  right_join(self_time_table1, by = c("full_ant_ID", "time_bins"))

shake_time_table2 <- select(shake_time_table2, -c("Replicate.y", "treatment.y", "group_size.y",
                                                  "ant_ID.y", "rep_set.y", "ant_status.y" ))


colnames(shake_time_table2) <- c("Replicate", "treatment", "group_size", "ant_ID", "time_bins",
                                 "rep_set", "Nb_shakes", "ant_status", "full_ant_ID", "Nb_allogrooming",
                                 "Nb_selfgrooming")




###############################################
##### NB BROOD TENDING PER MINUTE PER ANT #####
###############################################

##subset so only brood tending behaviours
brood_time_table <- subset(events, Behavior=="B")

##add rep set 
brood_time_table$rep_set <- "1"
brood_time_table [ which (brood_time_table$Replicate > 16 ) , "rep_set" ] <- "2"

##remove summer focals
events2_focals <- subset(brood_time_table, ant_status=="FOCAL")
events2_nestmates <- subset(brood_time_table, ant_status=="NESTMATE")
events2_focals_first_rep <- subset(events2_focals, rep_set=="1")
brood_time_table <- rbind(events2_nestmates, events2_focals_first_rep)

##add time bins
brood_time_table$Start_.s. <- as.numeric(brood_time_table$Start_.s.)
brood_time_table$time_bins <- round_any(brood_time_table$Start_.s., 60)

##add tally / sanity checks
brood_time_table$N <- 1
allo_per_min <- aggregate(N ~ group_size + treatment + ant_status + time_bins, FUN=sum, data=brood_time_table)

###create a data frame with all options to add zeros in
brood_bloc <- expand.grid( time_bins = unique ( brood_time_table$time_bins ))

brood_table <- NULL
for (repl in 1:nrow(replicate_list)){
  ants <- data.frame(ant_ID=ant_list [ 1 : replicate_list [repl, "group_size"] ] )
  repl_groom <- expand.grid.df( data.frame ( replicate_list[repl,])
                                , ants,
                                brood_bloc)
  brood_table <- rbind(brood_table, repl_groom)
}


brood_Nb <- aggregate(N ~ Replicate + treatment + group_size + ant_ID + time_bins + rep_set, data=brood_time_table, FUN=sum  )
brood_table1 <- merge(brood_table, brood_Nb, all.x=T)

##add in missing values
brood_table1 [is.na (brood_table1$N), "N" ] <- 0
brood_table1$ant_status <- "NESTMATE"
brood_table1 [ which ( brood_table1$ant_ID=="ANT_1"), "ant_status"] <- "FOCAL"
brood_table1$rep_set <- "1"
brood_table1 [ which (brood_table1$Replicate > 16 ) , "rep_set" ] <- "2"

##include a full ant ID column
brood_table1$full_ant_ID <-  with (brood_table1,paste(Replicate,treatment,group_size,ant_ID,sep="_"))

##combine to existing table
shake_time_table2 <- shake_time_table2 %>% 
  right_join(brood_table1, by = c("full_ant_ID", "time_bins"))

shake_time_table2 <- select(shake_time_table2, -c("Replicate.y", "treatment.y", "group_size.y",
                                                  "ant_ID.y", "rep_set.y", "ant_status.y" ))


colnames(shake_time_table2) <- c("Replicate", "treatment", "group_size", "ant_ID", "time_bins",
                                 "rep_set", "Nb_shakes", "ant_status", "full_ant_ID", "Nb_allogrooming",
                                 "Nb_selfgrooming", "Nb_broodtending")

###### result should be a table with tally for nb of all behaviours per minute per ant ######




#################################
##### GRAPHS OF TIME SERIES #####
#################################



allo <- aggregate(Nb_allogrooming ~ time_bins + treatment + ant_status, FUN=sum, data=shake_time_table2)
shk <- aggregate(Nb_shakes ~ time_bins + treatment + ant_status, FUN=sum, data=shake_time_table2)
self <- aggregate(Nb_selfgrooming ~ time_bins + treatment + ant_status, FUN=sum, data=shake_time_table2)

##barplots of shakes 
ggplot(shk, aes(x=time_bins, y=Nb_shakes, fill=treatment))+ 
  geom_col(position="dodge") +
  facet_wrap(~ant_status)

##barplots of allogrooming events
ggplot(allo, aes(x=time_bins, y=Nb_allogrooming, fill=treatment))+ 
  geom_col(position="dodge") +
  facet_wrap(~ant_status)

##barplots of self grooming events
ggplot(self, aes(x=time_bins, y=Nb_selfgrooming, fill=treatment))+ 
  geom_col(position="dodge") +
  facet_wrap(~ant_status)

##put all tallys into one dataframe
all <- allo %>% right_join(shk, by= c("time_bins", "ant_status", "treatment"))
all <- all %>% right_join(self, by = c ( "time_bins", "ant_status", "treatment"))


##need to separate graphs into the treatments
ggplot(all, aes(x=time_bins, y=Nb_shakes, color=treatment))+
  geom_point()






##############################
##### BINNED TIME SERIES #####
### SHAKING VS ALLOGROOMING ##
##############################

# install.packages("lmtest")
# install.packages("tseries")

library(lmtest)
library(tseries)


# ###granger test of causality
# grangertest(Nb_shakes ~ Nb_selfgrooming, order=1, data=shake_time_table2)



###for each trial we want to work out what the best lag time is 
## loop through each trial and calculate the cross correlation using ccf, extract the ccf values to determine which is best
## store the best values per trial and which lags
## granger test with these lag values?
## interpret the lags biologically

shake_time_table2$petri_dish <-  with (shake_time_table2,paste(Replicate,treatment,group_size,sep="_"))


optimal_lag_per_trial_binned_SA <- data.frame(acf=numeric(0), lag=character(0), trial=character(0), distance=numeric(0))
cross_cor_values = NULL


for (PETRI in unique(shake_time_table2$petri_dish)) {
  
  SHAKE_PETRI <- shake_time_table2[which(shake_time_table2$petri_dish==PETRI),]
  
  SHAKE_PETRI_S <- sum(SHAKE_PETRI$Nb_shakes)
  SHAKE_PETRI_A <- sum(SHAKE_PETRI$Nb_allogrooming)
  
  if (SHAKE_PETRI_S != 0 & SHAKE_PETRI_A != 0 ){
    
    cross_cor_values <- ccf(SHAKE_PETRI$Nb_shakes, SHAKE_PETRI$Nb_allogrooming)
    
    cross_cor_values <- as.data.frame(do.call(cbind, cross_cor_values))
    
    cross_cor_values <- cross_cor_values %>% select ( -one_of( 'type', 'series', 'n.used'))
    
    #cross_cor_values$trial <- sapply(strsplit(cross_cor_values$snames, " "), getElement, 3)
    
    #cross_cor_values$trial <- gsub('["),]', '', cross_cor_values$trial)
    
    cross_cor_values <- cross_cor_values %>% select(-snames)
    
    cross_cor_values$trial <- unique(SHAKE_PETRI$petri_dish)
    
    cross_cor_values$acf <- as.numeric(cross_cor_values$acf)
    
    cross_cor_values$distance <- abs ( 0 - cross_cor_values$acf )
    
    cross_cor_values <- cross_cor_values %>% arrange(desc(distance))
    
    optimal_lag <- cross_cor_values[1,]
  }
  
  else if (SHAKE_PETRI_A == 0 | SHAKE_PETRI_S == 0){
    next()
  }
  
  optimal_lag_per_trial_binned_SA <- rbind(optimal_lag_per_trial_binned_SA, optimal_lag)
  
}






################################
##### RAW DATA TIME SERIES #####
#### SHAKING VS ALLOGROOMING ###
################################



### same code but without the time bins ###

shake_time_series <- events[which(events$Behavior=="T"),] %>% 
  select(-one_of('Observation_id', 'Observation_date', 'Media_file', 'Total_length', 'FPS',
                 'Modifiers', 'Behavior_type', 'Stop_.s.', 'Duration_.s.'))

shake_time_series$petri_dish <- with ( shake_time_series, paste (Replicate, treatment, group_size, sep="_"))


allo_time_series <- events[which(events$Behavior=="A"),] %>%
  select(-one_of('Observation_id', 'Observation_date', 'Media_file', 'Total_length', 'FPS',
                 'Modifiers', 'Behavior_type', 'Stop_.s.', 'Duration_.s.'))

allo_time_series$petri_dish <- with ( allo_time_series, paste (Replicate, treatment, group_size, sep="_"))

shake_time_series$Start_.s. <- as.numeric(shake_time_series$Start_.s.)
allo_time_series$Start_.s. <- as.numeric(allo_time_series$Start_.s.)




time_series_data <- rbind(shake_time_series, allo_time_series)



optimal_lag_per_trial_raw_SA <- data.frame(acf=numeric(0), lag=character(0), trial=character(0), distance=numeric(0))
cross_cor_values = NULL

for (PETRI in unique(time_series_data$petri_dish)) {
  
  SHAKE_PETRI <- time_series_data[which(time_series_data$petri_dish==PETRI),]
  
  
  if (nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]) > 1 & nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="A"),]) > 1 ) {
    
    cross_cor_values <- ccf(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]$Start_.s., SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="A"),]$Start_.s. )
    
    cross_cor_values <- as.data.frame(do.call(cbind, cross_cor_values))
    
    cross_cor_values <- cross_cor_values %>% select ( -one_of( 'type', 'series', 'n.used', 'snames'))
    
    cross_cor_values$trial <- unique(SHAKE_PETRI$petri_dish)
    
    cross_cor_values$acf <- as.numeric(cross_cor_values$acf)
    
    cross_cor_values$distance <- abs ( 0 - cross_cor_values$acf )
    
    cross_cor_values <- cross_cor_values %>% arrange(desc(distance))
    
    optimal_lag <- cross_cor_values[1,]
    
  }
  
  
  else if (nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]) <= 1 | nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="A"),]) <= 1 ){
    next()
    
    
  }
  
  
  optimal_lag_per_trial_raw_SA <- rbind(optimal_lag_per_trial_raw_SA, optimal_lag)
  
}


### optimal_lag_per_trial contains 57 trials where both allogrooming and shaking happened





##############################
##### BINNED TIME SERIES #####
### SHAKING VS SELFGROOMING ##
##############################


optimal_lag_per_trial_binned_SS <- data.frame(acf=numeric(0), lag=character(0), trial=character(0), distance=numeric(0))
cross_cor_values = NULL

for (PETRI in unique(shake_time_table2$petri_dish)) {
  
  SHAKE_PETRI <- shake_time_table2[which(shake_time_table2$petri_dish==PETRI),]
  
  SHAKE_PETRI_S <- sum(SHAKE_PETRI$Nb_shakes)
  SHAKE_PETRI_G <- sum(SHAKE_PETRI$Nb_selfgrooming)
  
  if (SHAKE_PETRI_S != 0 & SHAKE_PETRI_G != 0 ){
    
    cross_cor_values <- ccf(SHAKE_PETRI$Nb_shakes, SHAKE_PETRI$Nb_selfgrooming)
    
    cross_cor_values <- as.data.frame(do.call(cbind, cross_cor_values))
    
    cross_cor_values <- cross_cor_values %>% select ( -one_of( 'type', 'series', 'n.used'))
    
    #cross_cor_values$trial <- sapply(strsplit(cross_cor_values$snames, " "), getElement, 3)
    
    #cross_cor_values$trial <- gsub('["),]', '', cross_cor_values$trial)
    
    cross_cor_values <- cross_cor_values %>% select(-snames)
    
    cross_cor_values$trial <- unique(SHAKE_PETRI$petri_dish)
    
    cross_cor_values$acf <- as.numeric(cross_cor_values$acf)
    
    cross_cor_values$distance <- abs ( 0 - cross_cor_values$acf )
    
    cross_cor_values <- cross_cor_values %>% arrange(desc(distance))
    
    optimal_lag <- cross_cor_values[1,]
  }
  
  else if (SHAKE_PETRI_S == 0 | SHAKE_PETRI_G == 0){
    next()
  }
  
  optimal_lag_per_trial_binned_SS <- rbind(optimal_lag_per_trial_binned_SS, optimal_lag)
  
}






################################
##### RAW DATA TIME SERIES #####
#### SHAKING VS SELFGROOMING ###
################################



# shake_time_series <- events[which(events$Behavior=="T"),] %>% 
#   select(-one_of('Observation_id', 'Observation_date', 'Media_file', 'Total_length', 'FPS',
#                  'Modifiers', 'Behavior_type', 'Stop_.s.', 'Duration_.s.'))
# 
# shake_time_series$petri_dish <- with ( shake_time_series, paste (Replicate, treatment, group_size, sep="_"))
# 
# shake_time_series$Start_.s. <- as.numeric(shake_time_series$Start_.s.)




self_time_series <- events[which(events$Behavior=="G"),] %>%
  select(-one_of('Observation_id', 'Observation_date', 'Media_file', 'Total_length', 'FPS',
                 'Modifiers', 'Behavior_type', 'Stop_.s.', 'Duration_.s.'))

self_time_series$petri_dish <- with ( self_time_series, paste (Replicate, treatment, group_size, sep="_"))


self_time_series$Start_.s. <- as.numeric(self_time_series$Start_.s.)



time_series_data <- rbind(shake_time_series, self_time_series)



par(mfrow=c(6,12))


## TO DO : EXCLUDE GROUP SIZE=1 from this analysis!

optimal_lag_per_trial_raw_SS <- data.frame(acf=numeric(0), lag=character(0), trial=character(0), distance=numeric(0))
cross_cor_values = NULL
cross_cor_values_ALL <- NULL

for (PETRI in unique(time_series_data$petri_dish)) {
  
  SHAKE_PETRI <- time_series_data[which(time_series_data$petri_dish==PETRI),]
  
  
  if (nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]) > 1 & nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="G"),]) > 3 ) {
    
    cross_cor_values <- ccf(x=SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]$Start_.s., 
                            y=SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="G"),]$Start_.s., lag.max=100, plot=F)
    
    cross_cor_values <- as.data.frame(do.call(cbind, cross_cor_values))
    plot(acf ~ lag, cross_cor_values); abline(v=0); abline(h=0)
    
    cross_cor_values_ALL <- rbind(cross_cor_values_ALL, 
                                  data.frame(PETRI, cross_cor_values))
    
    cross_cor_values <- cross_cor_values %>% select ( -one_of( 'type', 'series', 'n.used', 'snames'))
    
    cross_cor_values$trial <- unique(SHAKE_PETRI$petri_dish)
    
    cross_cor_values$acf <- as.numeric(cross_cor_values$acf)
    
    cross_cor_values$distance <- abs ( 0 - cross_cor_values$acf )
    
    cross_cor_values <- cross_cor_values %>% arrange(desc(distance))
    
    optimal_lag <- cross_cor_values[1,]
    
  }
  
  
  else if (nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="T"),]) <= 1 | nrow(SHAKE_PETRI[which(SHAKE_PETRI$Behavior=="G"),]) <= 1 ){
    
    next()
    
    
  }
  
  
  optimal_lag_per_trial_raw_SS <- rbind(optimal_lag_per_trial_raw_SS, optimal_lag)
  
}

cross_cor_values_ALL$acf <- as.numeric(cross_cor_values_ALL$acf)
cross_cor_values_ALL$snames <- NULL

## wrangling
cross_cor_values_ALL$REPLICATE <- sapply(strsplit(cross_cor_values_ALL$PETRI,"_"),"[[",1)
cross_cor_values_ALL$TREATMENT <- sapply(strsplit(cross_cor_values_ALL$PETRI,"_"),"[[",2)
cross_cor_values_ALL$GROUPSIZE <- sapply(strsplit(cross_cor_values_ALL$PETRI,"_"),"[[",3)

## averaging across combinations
cross_cor_values_MEANS <- aggregate( acf ~ lag + TREATMENT + GROUPSIZE, FUN=mean, na.rm=T,  cross_cor_values_ALL, na.action=NULL)
cross_cor_values_SD    <- aggregate( acf ~ lag + TREATMENT + GROUPSIZE, FUN=standard_error,     cross_cor_values_ALL, na.action=NULL)

# kludge
cross_cor_values_SD$acf [ which(is.na(cross_cor_values_SD$acf))] <- 0

##
cross_cor_values_MEANS$lag <- as.numeric(cross_cor_values_MEANS$lag)
cross_cor_values_SD$lag <- as.numeric(cross_cor_values_SD$lag)


## exclude group size 1
cross_cor_values_MEANS <- cross_cor_values_MEANS[which(cross_cor_values_MEANS$GROUPSIZE!="1"),]
cross_cor_values_SD <- cross_cor_values_SD[which(cross_cor_values_SD$GROUPSIZE!="1"),]

## plot eaach combination
par(mfrow=c(3,3))
for (TREAT in unique(cross_cor_values_MEANS$TREATMENT))
{
  for (GRP in unique(cross_cor_values_MEANS$GROUPSIZE))
  {
    ## show the averagec correlation
    Ylimmy <- c(-0.5,0.5)
    plot(acf ~ lag, cross_cor_values_MEANS[which(cross_cor_values_MEANS$TREATMENT==TREAT & cross_cor_values_MEANS$GROUPSIZE==GRP),], main=paste(TREAT,GRP), ylim=Ylimmy)
    ## show the sd around those averaves
    segments(x0 = cross_cor_values_SD$lag   [which(cross_cor_values_SD$TREATMENT==TREAT & cross_cor_values_SD$GROUPSIZE==GRP)],
             y0 = cross_cor_values_MEANS$acf[which(cross_cor_values_MEANS$TREATMENT==TREAT & cross_cor_values_MEANS$GROUPSIZE==GRP)] - cross_cor_values_SD$acf [which(cross_cor_values_SD$TREATMENT==TREAT & cross_cor_values_SD$GROUPSIZE==GRP)],
             x1 = cross_cor_values_SD$lag   [which(cross_cor_values_SD$TREATMENT==TREAT & cross_cor_values_SD$GROUPSIZE==GRP)],
             y1 = cross_cor_values_MEANS$acf[which(cross_cor_values_MEANS$TREATMENT==TREAT & cross_cor_values_MEANS$GROUPSIZE==GRP)] + cross_cor_values_SD$acf [which(cross_cor_values_SD$TREATMENT==TREAT & cross_cor_values_SD$GROUPSIZE==GRP)])
    abline(v=0, lty=2); abline(h=0, lty=2)
  }
}

