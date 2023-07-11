# 
# install.packages("cowplot")



library(gridExtra)
library(cowplot)
library(igraph)


####function
get_posthoc_groups <- function(model,contrast_matrix,which_levels,dataset,interaction_variables=NULL){
  levels <- get(which_levels)
  if (is.null(names(levels))){
    names(levels) <- paste("Delta_",levels,sep="")
  }
  
  post_hoc <- summary(glht(model,contrast_matrix),test=adjusted("BH"))
  # print("z value");print(post_hoc$test$tstat);print("Pr>|z|");print(post_hoc$test$pvalues);
  p_values <- as.numeric(post_hoc$test$pvalues);names(p_values) <- names(post_hoc$test$coefficients)
  
  post_hoc_levels <- names(levels)
  post_hoc_mat <- matrix(NA,nrow=length(post_hoc_levels)-1,ncol=length(post_hoc_levels)-1)
  rownames(post_hoc_mat) <- post_hoc_levels[2:length(post_hoc_levels)]
  colnames(post_hoc_mat) <- post_hoc_levels[1:(length(post_hoc_levels)-1)]
  for (i in 1:nrow(post_hoc_mat)){
    for (j in 1:i){
      post_hoc_mat[i,j] <- as.logical(as.numeric(p_values[paste(colnames(post_hoc_mat)[j]," - ",rownames(post_hoc_mat)[i],sep="")])>0.05)
    }
  }
  g <- post_hoc_mat
  g <- cbind(rbind(NA, g), NA)
  g <- replace(g, is.na(g), FALSE)
  g <- g + t(g)
  diag(g) <- 1
  n <- length(post_hoc_levels)
  rownames(g) <- 1:n
  colnames(g) <- 1:n
  #g
  same <- which(g==1)
  topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
  topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
  g3 <- simplify(graph.data.frame(topology,directed = FALSE))
  #get.data.frame(g3)
  #plot(g3)
  res <- maximal.cliques(g3)
  
  # Reorder given the smallest level
  clique_value <- NULL
  if (which_levels=="groups"){
    form <- as.formula(paste("variable ~ ", paste(interaction_variables, collapse= "+")))
  }else{
    form <- as.formula(paste("variable ~ ", paste(gsub("_order","",which_levels), collapse= "+")))
  }
  
  means <- aggregate(form,FUN=mean,data=dataset)
  means$predictor <- names(levels)[match(paste(means[,interaction_variables[1]],means[,interaction_variables[2]],sep="_"),levels)]
  for (i in 1:length(res)){
    clique_value <- c(clique_value,mean(means[as.numeric(unlist(res[[i]])),"variable"]))
  }
  res <- res[order(clique_value)]
  
  # Get group letters
  lab.txt <- vector(mode="list", n)
  lab <- letters[seq(res)]
  for(i in seq(res)){
    for(j in res[[i]]){
      lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
    }
  }
  post_hoc_groups <- unlist(lab.txt); names(post_hoc_groups) <- levels[post_hoc_levels]
  return(post_hoc_groups)
}



########################################
##### NB SHAKES PER MINUTE PER ANT #####
########################################



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


##sanity checks##

ggplot(shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),],
       aes(x=treatment, y=N, fill=treatment))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~smins)

ggplot(shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),],
       aes(x=group_size, y=N, fill=treatment))+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~smins)




###### GRAPHS FOR MEAN BODY SHAKES PER ANT ###### 


### FOCAL ANTS ###



mean_shake_mins <- aggregate(frequency  ~ treatment + ant_status + smins, FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~   treatment + ant_status + smins, FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c("Treatment", "ant_status", "time_bins", "mean_shakes", "se")


mean_shake_mins[which(mean_shake_mins$Treatment=="C"), "Treatment"] <- "Control"
mean_shake_mins[which(mean_shake_mins$Treatment=="S"), "Treatment"] <- "Sham"
mean_shake_mins[which(mean_shake_mins$Treatment=="P"), "Treatment"] <- "Pathogen"

mean_shake_mins <- subset(mean_shake_mins, ant_status=="FOCAL")
mean_shake_mins1 <- subset(mean_shake_mins, time_bins=="0-5mins")

a <- ggplot(mean_shake_mins1, aes(fill=Treatment, y=mean_shakes, x=Treatment))+
  geom_errorbar(aes(ymin=mean_shakes-0.001, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  # facet_wrap(~ time_bins)+
  labs(x="Treatment", y="Mean number of shakes per ant", title="0-5 minutes")+
  theme_bw()+
  guides(fill="none")+
  geom_signif(comparisons = list(c("Control", "Pathogen")), annotations="***",
              y_position = 0.11, tip_length=0.025, vjust=0.4,  size=1)+
  geom_signif(comparisons = list(c("Control", "Sham")), annotations="***",
              y_position = 0.13, tip_length=0.025, vjust=0.4,  size=1)+
  ylim(c(0, 0.15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_blank(),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_text(size=18), legend.position = "bottom" ,  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(size=18,  hjust=0.5))

leg <- get_legend(a)

mean_shake_mins2 <- subset(mean_shake_mins, time_bins=="5-10mins")

b <- ggplot(mean_shake_mins2, aes(fill=Treatment, y=mean_shakes, x=Treatment))+
  geom_errorbar(aes(ymin=mean_shakes-0.001, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  # facet_wrap(~ time_bins)+
  labs(x="Treatment", y="Mean number of shakes per ant", title="5-10 minutes")+
  theme_bw()+
  guides(fill="none")+
  geom_signif(comparisons = list(c("Control", "Pathogen")), annotations="**",
              y_position = 0.06, tip_length=0.05, vjust=0.4,  size=1)+
  geom_signif(comparisons = list(c("Control", "Sham")), annotations="*",
              y_position = 0.08, tip_length=0.05, vjust=0.4,  size=1)+
  ylim(c(0, 0.15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_blank(),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(size=18,  hjust=0.5))

mean_shake_mins3 <- subset(mean_shake_mins, time_bins=="10-15mins")


c <- ggplot(mean_shake_mins3, aes(fill=Treatment, y=mean_shakes, x=Treatment))+
  geom_errorbar(aes(ymin=mean_shakes-0.001, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  # facet_wrap(~ time_bins)+
  labs(x="Treatment", y="Mean number of shakes per ant", title="10-15 minutes")+
  theme_bw()+
  guides(fill="none")+
  geom_signif(comparisons = list(c("Control", "Pathogen")), annotations="*",
              y_position = 0.04, tip_length=0.1, vjust=0.4,  size=1)+
  ylim(c(0, 0.15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_blank(),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(size=18,  hjust=0.5))


mean_shake_mins4 <- subset(mean_shake_mins, time_bins=="15-20mins")




d <- ggplot(mean_shake_mins4, aes(fill=Treatment, y=mean_shakes, x=Treatment))+
  geom_errorbar(aes(ymin=mean_shakes-0.001, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  # facet_wrap(~ time_bins)+
  labs(x="Treatment", y="Mean number of shakes per ant", title="15-20 minutes")+
  theme_bw()+
  guides(fill="none")+
  geom_signif(comparisons = list(c("Control", "Pathogen")), annotations="*",
              y_position = 0.03, tip_length=0.1, vjust=0.4,  size=1)+
  ylim(c(0, 0.15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_blank(),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.title = element_text(size=18, hjust=0.5))



hlay <- rbind(c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(1,1,2,2,3,3,4,4),
              c(NA,NA,NA,5,5,NA,NA,NA))


tg <- text_grob("Treatment", just='centre', family="Roboto", size=18)
                
allplots <- align_plots(a,b,c,d, align="h")

#grid.arrange(allplots[[1]], allplots[[2]], allplots[[3]], allplots[[4]], leg, layout_matrix=hlay)

plot_grid(allplots[[1]], allplots[[2]], allplots[[3]], allplots[[4]], leg,
          ncol=4, nrow=2, rel_widths = c(0.28,0.24,0.24,0.24),
          rel_heights = c(0.9,0.1)) 


plot1 <- plot_grid(allplots[[1]], allplots[[2]], allplots[[3]], allplots[[4]],
          ncol=4, rel_widths = c(0.28,0.24,0.24,0.24))


plot_grid(plot1, leg, ncol=1, rel_heights = c(0.9, 0.1))


### NESTMATES ###


##time bins only

mean_shake_mins <- aggregate(frequency  ~ smins + ant_status , FUN=mean, data=shake_time_table1)
se_shake_mins <- aggregate(frequency ~   smins + ant_status , FUN=standard_error, data=shake_time_table1)

mean_shake_mins$se <- se_shake_mins$frequency
colnames(mean_shake_mins) <- c("time_bins", "ant_status", "mean_shakes", "se")

mean_shake_mins <- subset(mean_shake_mins, ant_status=="NESTMATE")
mean_shake_mins$time_bins <- as.character(mean_shake_mins$time_bins)

mean_shake_mins[which(mean_shake_mins$time_bins=="0-5mins"), "time_bins"] <- "0 - 5"
mean_shake_mins[which(mean_shake_mins$time_bins=="5-10mins"), "time_bins"] <- "5 - 10"
mean_shake_mins[which(mean_shake_mins$time_bins=="10-15mins"), "time_bins"] <- "10 - 15"
mean_shake_mins[which(mean_shake_mins$time_bins=="15-20mins"), "time_bins"] <- "15 - 20"

mean_shake_mins$time_bins <- factor(mean_shake_mins$time_bins,
                                       levels=c("0 - 5", "5 - 10", "10 - 15",  "15 - 20"))


cld <- c("c", "b", "ab", "a")
ggplot(mean_shake_mins, aes(fill=time_bins, y=mean_shakes, x=time_bins))+
  geom_errorbar(aes(ymin=mean_shakes-se, ymax=mean_shakes+se), width=.6, position=position_dodge(0.9), size=1) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label= cld, y = mean_shakes + se), vjust= -1, position = position_dodge(0.9), size=6, family="Roboto") +
  # facet_wrap(~ time_bins)+
  labs(x="Time bins (minutes)", y="Mean number of shakes per ant", title="Mean number of shakes by nestmates")+
  theme_bw()+
  guides(fill="none")+
  scale_fill_viridis_d(option = "plasma")+
  # geom_signif(comparisons = list(c("Control", "Pathogen")), annotations="*",
  #             y_position = 0.03, tip_length=0.1, vjust=0.4,  size=1)+
  ylim(c(0, 0.04))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=18),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15, colour="black"), strip.text=element_text(size=15),
        axis.title.y=element_text(size=18), plot.title = element_text(size=18, hjust=0.5) )




### STATISTICS ###

shake_time_table1$interac_size_treat <- with (shake_time_table1, interac_size_treat <- paste(group_size,treatment,sep="_"))
shake_time_table1$Petri_dish <- with (shake_time_table1, interac_size_treat <- paste(Replicate,group_size,treatment,sep="_"))

shake_time_table1$Replicate <- as.numeric(shake_time_table1$Replicate)
shake_time_table1$time <- "WINTER"
shake_time_table1 [ which (shake_time_table1$Replicate > 16 ) , "time" ] <- "SUMMER"



### FOCAL MODEL ###
model_focal <- glmer(N^0.1 ~ group_size*treatment*smins + offset(log(duration)) + (1|Replicate),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),])

#remove non sig 3 way interaction
model_focal <- glmer(N^0.1 ~ group_size*treatment+group_size*smins+ treatment*smins+ offset(log(duration)) + (1|Replicate),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),])

model_focal <- glmer(N^0.1 ~ group_size*treatment+ treatment*smins+ offset(log(duration)) + (1|Replicate),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),])


qqnorm(residuals(model_focal))
qqline(residuals(model_focal))
hist(residuals(model_focal))

test_norm(residuals(model_focal))

Anova(model_focal, type="III")

# post_hoc_smins <- summary(glht(model_focal, linfct=mcp(smins="Tukey")), test=adjusted("BH"))
contrast_matrix <- rbind(
  "1 - 2"=c(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "1 - 3"=c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "1 - 4"=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "1 - 5"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0)
  ,
  "1 - 6"=c(0,0,0,0,0,-1,-1,0,0,0,0,0,-0.25,-0.25,-0.25,0,-1,0,0,0,0)
  ,
  "1 - 7"=c(0,0,0,0,0,-1,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0)
  ,
  "1 - 8"=c(0,0,0,0,0,-1,0,0,-1,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1)
  ,
  "1 - 9"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,0,0)
  ,
  "1 - 10"=c(0,0,0,0,-1,0,-1,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0,0,0,0)
  ,
  "1 - 11"=c(0,0,0,0,-1,0,0,-1,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1,0,0,0)
  ,
  "1 - 12"=c(0,0,0,0,-1,0,0,0,-1,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,-1,0)
  ,
  "2 - 3"=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "2 - 4"=c(0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "2 - 5"=c(0,0,0,0,0,-1,1,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0)
  ,
  "2 - 6"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,-1,0,0,0,0)
  ,
  "2 - 7"=c(0,0,0,0,0,-1,1,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0)
  ,
  "2 - 8"=c(0,0,0,0,0,-1,1,0,-1,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1)
  ,
  "2 - 9"=c(0,0,0,0,-1,0,1,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,0,0)
  ,
  "2 - 10"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0,0,0,0)
  ,
  "2 - 11"=c(0,0,0,0,-1,0,1,-1,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1,0,0,0)
  ,
  "2 - 12"=c(0,0,0,0,-1,0,1,0,-1,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,-1,0)
  ,
  "3 - 4"=c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0)
  ,
  "3 - 5"=c(0,0,0,0,0,-1,0,1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0)
  ,
  "3 - 6"=c(0,0,0,0,0,-1,-1,1,0,0,0,0,-0.25,-0.25,-0.25,0,-1,0,0,0,0)
  ,
  "3 - 7"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0)
  ,
  "3 - 8"=c(0,0,0,0,0,-1,0,1,-1,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1)
  ,
  "3 - 9"=c(0,0,0,0,-1,0,0,1,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,0,0)
  ,
  "3 - 10"=c(0,0,0,0,-1,0,-1,1,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0,0,0,0)
  ,
  "3 - 11"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1,0,0,0)
  ,
  "3 - 12"=c(0,0,0,0,-1,0,0,1,-1,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,-1,0)
  ,
  "4 - 5"=c(0,0,0,0,0,-1,0,0,1,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0)
  ,
  "4 - 6"=c(0,0,0,0,0,-1,-1,0,1,0,0,0,-0.25,-0.25,-0.25,0,-1,0,0,0,0)
  ,
  "4 - 7"=c(0,0,0,0,0,-1,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0)
  ,
  "4 - 8"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1)
  ,
  "4 - 9"=c(0,0,0,0,-1,0,0,0,1,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,0,0)
  ,
  "4 - 10"=c(0,0,0,0,-1,0,-1,0,1,-0.25,-0.25,-0.25,0,0,0,-1,0,0,0,0,0)
  ,
  "4 - 11"=c(0,0,0,0,-1,0,0,-1,1,-0.25,-0.25,-0.25,0,0,0,0,0,-1,0,0,0)
  ,
  "4 - 12"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,-1,0)
  ,
  "5 - 6"=c(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0)
  ,
  "5 - 7"=c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0)
  ,
  "5 - 8"=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1)
  ,
  "5 - 9"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,0,0)
  ,
  "5 - 10"=c(0,0,0,0,-1,1,-1,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,-1,0,0,0,0,0)
  ,
  "5 - 11"=c(0,0,0,0,-1,1,0,-1,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,-1,0,0,0)
  ,
  "5 - 12"=c(0,0,0,0,-1,1,0,0,-1,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,-1,0)
  ,
  "6 - 7"=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,1,0,-1,0,0)
  ,
  "6 - 8"=c(0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,0,1,0,0,0,-1)
  ,
  "6 - 9"=c(0,0,0,0,-1,1,1,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,1,0,0,0,0)
  ,
  "6 - 10"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,-1,1,0,0,0,0)
  ,
  "6 - 11"=c(0,0,0,0,-1,1,1,-1,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,1,-1,0,0,0)
  ,
  "6 - 12"=c(0,0,0,0,-1,1,1,0,-1,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,1,0,0,-1,0)
  ,
  "7 - 8"=c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,1,0,-1)
  ,
  "7 - 9"=c(0,0,0,0,-1,1,0,1,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,1,0,0)
  ,
  "7 - 10"=c(0,0,0,0,-1,1,-1,1,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,-1,0,0,1,0,0)
  ,
  "7 - 11"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,-1,1,0,0)
  ,
  "7 - 12"=c(0,0,0,0,-1,1,0,1,-1,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,1,-1,0)
  ,
  "8 - 9"=c(0,0,0,0,-1,1,0,0,1,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,0,1)
  ,
  "8 - 10"=c(0,0,0,0,-1,1,-1,0,1,-0.25,-0.25,-0.25,0.25,0.25,0.25,-1,0,0,0,0,1)
  ,
  "8 - 11"=c(0,0,0,0,-1,1,0,-1,1,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,-1,0,0,1)
  ,
  "8 - 12"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,-1,1)
  ,
  "9 - 10"=c(0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0)
  ,
  "9 - 11"=c(0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0)
  ,
  "9 - 12"=c(0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0)
  ,
  "10 - 11"=c(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,1,0,-1,0,0,0)
  ,
  "10 - 12"=c(0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,1,0,0,0,-1,0)
  ,
  "11 - 12"=c(0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,1,0,-1,0)
)

 post_hoc_interac_smins_treatment <- summary(glht(model_focal, linfct=contrast_matrix), test=adjusted("BH"))
 
 groups <- paste(rep(c("C","S","P"),each=4),rep(c("0-5mins","5-10mins","10-15mins","15-20mins"),3),sep="_")
 names(groups) <- as.character(1:12)
 dataset = shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),]
 dataset$variable <- dataset$N^0.1
 get_posthoc_groups(model=model_focal,contrast_matrix=contrast_matrix,which_levels="groups",dataset=dataset,interaction_variables=c("treatment","smins"))
 
 
 
 
 contrast_matrix_simplified <- rbind(
   "1 - 5"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0)
   ,
   "1 - 9"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,0,0)
   ,
   "5 - 9"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,0,0)
   ,
   "2 - 6"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,-1,0,0,0,0)
   ,
   "2 - 10"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0,0,0,0)
   ,
   "6 - 10"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,-1,1,0,0,0,0)
   ,
   "3 - 7"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,-1,0,0)
   ,
   "3 - 11"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1,0,0,0)
   ,
   "7 - 11"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,-1,1,0,0)
   ,
   "4 - 8"=c(0,0,0,0,0,-1,0,0,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,-1)
   ,
   "4 - 12"=c(0,0,0,0,-1,0,0,0,0,-0.25,-0.25,-0.25,0,0,0,0,0,0,0,-1,0)
,
  "8 - 12"=c(0,0,0,0,-1,1,0,0,0,-0.25,-0.25,-0.25,0.25,0.25,0.25,0,0,0,0,-1,1)
 ) 
 post_hoc_interac_smins_treatment_simplified <- summary(glht(model_focal, linfct=contrast_matrix_simplified), test=adjusted("BH"))
 
 
 
 groups <- paste(rep(c("C","S","P"),each=4),rep(c("0-5mins","5-10mins","10-15mins","15-20mins"),3),sep="_")
 names(groups) <- as.character(1:12)
 dataset = shake_time_table1[which(shake_time_table1$ant_status=="FOCAL"),]
 dataset$variable <- dataset$N^0.1
 get_posthoc_groups(model=model_focal,contrast_matrix=contrast_matrix_simplified,which_levels="groups",dataset=dataset,interaction_variables=c("treatment","smins"))
 
 
 
 
 
 
 

### NESTMATES MODEL ###

model_nestmates <- glmer(N^0.1 ~ group_size*treatment*smins + offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
                     family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),])

#remove non sig 3 way interaction
model_nestmates <- glmer(N^0.1 ~ group_size*treatment+group_size*smins+ treatment*smins+ offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
                           family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),])


#remove next largest p value interaction
model_nestmates <- glmer(N^0.1 ~ group_size*smins+ treatment*smins+ offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
                         family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),])


#remove next largest p value interaction 
model_nestmates <- glmer(N^0.1 ~ group_size + smins + treatment*smins + offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
                         family="gaussian", data=shake_time_table1[which(shake_time_table1$ant_status=="NESTMATE"),])

#remove last interaction
model_nestmates <- glmer(N^0.1 ~ group_size + smins + treatment + offset(log(duration)) + (1|Replicate) + (1|Petri_dish) + (1|time),
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



