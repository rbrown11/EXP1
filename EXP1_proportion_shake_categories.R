
##############################################################
############# PROPORTIONS OF SHAKES OF EACH TYPE #############
##############################################################





####### REPETITIVE SHAKES GRAPH ####### 


sum_rep_shakes <- aggregate(N ~  twitch_repetition + ant_status + treatment + group_size, FUN=sum, data=twitch_table)
sum_shakes <- aggregate(N ~ ant_status + treatment + group_size, FUN=sum, data=twitch_table)


rep_shakes <- merge(sum_rep_shakes, sum_shakes, by=c("ant_status", "treatment", "group_size"), all=F)
colnames(rep_shakes) <- c("ant_status", "Treatment", "group_size", "Repetition", "Nb_shakes", "Total_shakes")


rep_shakes$proportion <- rep_shakes$Nb_shakes / rep_shakes$Total_shakes


rep_shakes$group_size <- factor(rep_shakes$group_size,
                                   levels=c("1", "2", "6", "11"))

rep_shakes[which(rep_shakes$Treatment=="C"), "Treatment"] <- "Control"
rep_shakes[which(rep_shakes$Treatment=="S"), "Treatment"] <- "Sham"
rep_shakes[which(rep_shakes$Treatment=="P"), "Treatment"] <- "Pathogen"

rep_shakes$Repetition <- as.character(rep_shakes$Repetition)


rep_shakes[which(rep_shakes$Repetition=="single"), "Repetition"] <- "Single"
rep_shakes[which(rep_shakes$Repetition=="repetitive"), "Repetition"] <- "Repetitive"

rep_shakes$Repetition <- factor(rep_shakes$Repetition,
                                levels=c("Single", "Repetitive"))

ggplot(rep_shakes)+
  geom_bar(aes(x=group_size, y=proportion,
               fill=Treatment, color=Treatment, alpha=Repetition), 
           stat="identity")+
  facet_grid(Treatment~ant_status) +
  scale_alpha_manual(values = c(0.5, 1)) +
  # scale_alpha_discrete(aes(values=c(0, 1)))+
  labs(x="Group size", y="Proportion of total body shakes", title="REPETITION")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )







####### SIZE SHAKES GRAPH ####### 


sum_size_shakes <- aggregate(N ~  twitch_size + ant_status + treatment + group_size, FUN=sum, data=twitch_table)
sum_shakes <- aggregate(N ~ ant_status + treatment + group_size, FUN=sum, data=twitch_table)


size_shakes <- merge(sum_size_shakes, sum_shakes, by=c("ant_status", "treatment", "group_size"), all=F)
colnames(size_shakes) <- c("ant_status", "Treatment", "group_size", "Repetition", "Nb_shakes", "Total_shakes")


size_shakes$proportion <- size_shakes$Nb_shakes / size_shakes$Total_shakes


size_shakes$group_size <- factor(size_shakes$group_size,
                                levels=c("1", "2", "6", "11"))

size_shakes[which(size_shakes$Treatment=="C"), "Treatment"] <- "CONTROL"
size_shakes[which(size_shakes$Treatment=="S"), "Treatment"] <- "SHAM"
size_shakes[which(size_shakes$Treatment=="P"), "Treatment"] <- "PATHOGEN"

size_shakes$Repetition <- as.character(size_shakes$Repetition)


size_shakes[which(size_shakes$Repetition=="large"), "Repetition"] <- "Large"
size_shakes[which(size_shakes$Repetition=="small"), "Repetition"] <- "Small"

size_shakes$Size <- factor(size_shakes$Repetition,
                                levels=c("Small", "Large"))

ggplot(size_shakes)+
  geom_bar(aes(x=group_size, y=proportion,
               fill=Treatment, color=Treatment, alpha=Size), 
           stat="identity")+
  facet_grid(Treatment~ant_status) +
  scale_alpha_manual(values = c(0.5, 1)) +
  # scale_alpha_discrete(aes(values=c(0, 1)))+
  labs(x="Group size", y="Proportion of total body shakes", title="SIZE")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(family="Roboto"), axis.title = element_text(size=20),
        legend.title=element_text(size=20), legend.text=element_text(size=18), axis.text=element_text(size=15), strip.text=element_text(size=15),
        axis.title.y=element_text(size=16) )









############################ 
######## STATISTICS ######## 
############################ 


##duration?


total_rep_shakes <- aggregate(N ~ full_ant_ID +  twitch_repetition + ant_status + treatment + group_size + Replicate
                              + time + duration, FUN=sum, data=total_twitches)
total_shakes <- aggregate(N ~ full_ant_ID + ant_status + treatment + group_size + Replicate
                          + time + duration, FUN=sum, data=total_twitches)



total_rep_shakes1 <- merge(total_rep_shakes, total_shakes, by=c("ant_status", "treatment", "group_size", "Replicate", "full_ant_ID",
                                                                "time", "duration"), all=F)
colnames(total_rep_shakes1) <- c("ant_status", "Treatment", "group_size", "Replicate", "full_ant_ID", 
"time", "duration", "Repetition", "Nb_shakes", "Total_shakes")


total_rep_shakes1$proportion <- total_rep_shakes1$Nb_shakes / total_rep_shakes1$Total_shakes

total_rep_shakes1 [ which (is.nan(total_rep_shakes1$proportion)), "proportion" ] <- 0

total_rep_shakes1_focal <- subset(total_rep_shakes1, ant_status=="FOCAL")

focal_model <- glmer(proportion ~ Treatment*group_size*Repetition + (1|Replicate),
                     weights=Total_shakes, family='binomial', data=total_rep_shakes1_focal)



Anova(focal_model, type="III")

#no significant interactions

focal_model_ni <- glmer(proportion ~ Treatment+group_size+Repetition + (1|Replicate),
                     weights=Total_shakes, family='binomial', data=total_rep_shakes1_focal)

Anova(focal_model_ni)



