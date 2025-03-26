library(readr)
library(dplyr)
library(lme4)
library(emmeans)
library(effects)
library(car)
setwd("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut")
#setwd("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut")
mean_velocity <- read_csv("mean_velocity 2.csv", 
                          col_types = cols(ID = col_number()))

veldf<-aggregate(x = mean_velocity, by = list(mean_velocity$ID, mean_velocity$condition, mean_velocity$body_part), FUN = "mean")

colnames(veldf)<-list("ID", "condition", "bodypart", "skip", "skip", "skip", "mean_v")

veldf<- select(veldf, ID, condition, bodypart, mean_v)

veldf$ID<-as.factor(veldf$ID)



library(tidyr)
veldf <- separate(veldf, bodypart, into = c("bodypart", "dimension"), sep = "_")

veldf.music <- veldf[veldf$condition != 'si', ]                          # Remove row based on condition
#veldf.music$age <- as.factor(veldf.music$age)
veldf.music$dimension <- as.factor(veldf.music$dimension)

setwd("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R")
musicom_id_age <- read_delim("musicom_id_age.csv", 
                             delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                             trim_ws = TRUE)
musicom_id_age$ID<-as.factor(musicom_id_age$ID)
musicom_id_age$age<-as.factor(musicom_id_age$age)

# remove duplicates in age dataframe
musicom_id_age <- musicom_id_age %>%
  distinct(ID, .keep_all = TRUE)

veldf.f <- merge(veldf.music, musicom_id_age, by='ID', all.x=T)


# Define scaling function
scale_value <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Scale value within each age group
veldf.f <- veldf.f %>%
  group_by(age) %>%
  group_by(bodypart) %>%
  group_by(dimension) %>%
  mutate(mean.v.sc = scale_value(mean_v))

####
# create roi for body parts
# summarize head(eyes, mouth)
# summarize hands(wrist, hand)
# summarize feet(ankle, foot))
# knees and elbows separately

# quantity of movement in one plot, age groups in different labels and conditions in different colors

veldf.f <- veldf.f[complete.cases(veldf.f), ]

#####
library(plyr)

veldf.f$condition.n<-revalue(veldf.f$condition, c("bs"="Music", "ct"="Control","hv"="High Pitch","lb"="Low Pitch"))
veldf.f$age.n<-revalue(veldf.f$age, c("3"="3m", "6"="6m","12"="12m"))

library(DescTools)
veldf.f$vel.w<-Winsorize(veldf.f$mean.v.sc,probs = c(0.025, 0.925))
# veldf.f$vel.wo <- remove_outliers(veldf.f$mean.v.sc, sd_threshold = 3)

# Calculate mean and standard deviation per age group
grouped_stats <- aggregate(mean_v ~ age.n, veldf.f, function(x) c(mean = mean(x), sd = sd(x)))

# Merge grouped_stats with original data
merged_data <- merge(veldf.f, grouped_stats, by = "age.n", suffixes = c("", "_stats"))

# Set outliers to mean + 3 * sd within each age group
merged_data$mean.v.sc <- with(merged_data, ifelse(abs(mean.v.sc - mean.v.sc_stats[,"mean"]) > 3 * mean.v.sc_stats[,"sd"], mean.v.sc_stats[,"mean"] + 3 * mean.v.sc_stats[,"sd"], mean.v.sc))
merged_data$mean.v.sc <- with(merged_data, ifelse(abs(mean.v.sc - mean.v.sc_stats[,"mean"]) < 3 * mean.v.sc_stats[,"sd"], mean.v.sc_stats[,"mean"] - 3 * mean.v.sc_stats[,"sd"], mean.v.sc))

# Remove extra columns
veldf.f <- merged_data[, !(names(merged_data) %in% c("mean", "sd"))]


# Filter the DataFrame based on conditions using dplyr
veldf.f.mvsc <- veldf.f %>%
  filter(condition.n %in% c("Music", "Control")) 
# Filter the DataFrame based on conditions using dplyr
veldf.f.hvsl <- veldf.f %>%
  filter(condition.n %in% c("High Pitch", "Low Pitch")) 

m1 <- lmer(mean_v ~ condition.n*age.n*bodypart+dimension  + (1|ID), 
           data=veldf.f.mvsc)
Anova(m1)

plot(allEffects(m1))
plot(effect("condition.n", m1), ylab="Velocity (scaled)", xlab="Conditions", main='Main Effect: Condition')
plot(effect("condition.n:age.n:dimension", m1),ylab="Velocity (scaled)", xlab="Conditions", main='Interaction Effect: Condition x Age' )

plot(effect("age:bodypart", m1))
plot(effect("age:dimension", m1))

emmeans(m1, pairwise~condition.n, adjust='tukey')
emmeans(m1, pairwise~condition.n|age.n, adjust='tukey')
emmeans(m1, pairwise~condition|bodypart, adjust='tukey')


#veldf.plot <- subset(veldf.roi, bodypart == "head" |  bodypart == "torso")

# Specify the order of conditions
condition_order <- c("Music", "Control") 

# Reorder the factor levels
veldf.f.mvsc$condition.n <- factor(veldf.f.mvsc$condition.n, levels = condition_order)
condtition_colors <- c('#FF9966','#999966')



library(ggplot2)
library(ggpubr)

ggplot(veldf.f.mvsc, aes(x = condition.n, y = mean.v.sc))+#, fill = condition.n)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(color = condition.n), alpha = 0.3)+#, position = position_jitterdodge(dodge.width = 0.75)) +
  geom_violin(aes(fill = condition.n),)+
  geom_boxplot(aes(fill = condition.n), alpha = 0.3)+
  scale_color_manual(values = condtition_colors) +
  scale_fill_manual(values = condtition_colors) +
  facet_wrap(~ age.n  , nrow= 3) +
  labs(x = "Condition", y = "Velocity") +
  #scale_fill_manual(values = c("blue", "red")) +
  #theme(legend='')+
  theme_bw()+
  theme(text = element_text(size = 16), legend.position = "none" )+
  ylim(0.00,0.30)


condition_order <- c("High Pitch", "Low Pitch")
veldf.f.hvsl$condition.n <- factor(veldf.f.hvsl$condition.n, levels = condition_order)
condtition_colors <- c('#4DBEEE','#993399')

ggplot(veldf.f.hvsl, aes(x = condition.n, y = mean.v.sc))+#, fill = condition.n)) +
  #geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(aes(color = condition.n), alpha = 0.3)+#, position = position_jitterdodge(dodge.width = 0.75)) +
  geom_violin(aes(fill = condition.n))+
  geom_boxplot(aes(fill = condition.n))+
  scale_color_manual(values = condtition_colors) +
  scale_fill_manual(values = condtition_colors) +
  facet_wrap(  ~ age.n, nrow = 3) +
  labs(x = "Condition", y = "Velocity") +
  #scale_fill_manual(values = c("blue", "red")) +
  #theme(legend='')+
  theme_bw()+
  theme(text = element_text(size = 16), legend.position = "none" )+
  ylim(0.00,0.30)



###### create ROIs
# head
veldf.head <- subset(veldf.f, bodypart == "lefteye" | bodypart == "righteye" | bodypart == "mouth")
veldf.head<-aggregate(x = veldf.head, by = list(veldf.head$ID, veldf.head$condition, veldf.head$age, veldf.head$dimension), FUN = "mean")
colnames(veldf.head)<-list("ID", "condition", "age","dimension" ,"skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.head<- select(veldf.head, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.head$bodypart <- "head"

# torso
veldf.torso <- subset(veldf.f, bodypart == "leftshoulder" | bodypart == "rightshoulder" | bodypart == "chest")
veldf.torso<-aggregate(x = veldf.torso, by = list(veldf.torso$ID, veldf.torso$condition, veldf.torso$age, veldf.torso$dimension), FUN = "mean")
colnames(veldf.torso)<-list("ID", "condition", "age","dimension" , "skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.torso<- select(veldf.torso, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.torso$bodypart <- "torso"

# left hand
veldf.lhand <- subset(veldf.f, bodypart == "rightwrist" | bodypart == "righthand")
veldf.lhand<-aggregate(x = veldf.lhand, by = list(veldf.lhand$ID, veldf.lhand$condition, veldf.lhand$age, veldf.lhand$dimension), FUN = "mean")
colnames(veldf.lhand)<-list("ID", "condition", "age","dimension" , "skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.lhand<- select(veldf.lhand, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.lhand$bodypart <- "lefthand"

# right hand
veldf.rhand <- subset(veldf.f, bodypart == "leftwrist" | bodypart == "lefthand")
veldf.rhand<-aggregate(x = veldf.rhand, by = list(veldf.rhand$ID, veldf.rhand$condition, veldf.rhand$age, veldf.rhand$dimension), FUN = "mean")
colnames(veldf.rhand)<-list("ID", "condition", "age","dimension" , "skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.rhand<- select(veldf.rhand, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.rhand$bodypart <- "righthand"

# left foot
veldf.lfoot <- subset(veldf.f, bodypart == "rightankle" | bodypart == "rightfoot")
veldf.lfoot<-aggregate(x = veldf.lfoot, by = list(veldf.lfoot$ID, veldf.lfoot$condition, veldf.lfoot$age, veldf.lfoot$dimension), FUN = "mean")
colnames(veldf.lfoot)<-list("ID", "condition", "age","dimension" , "skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.lfoot<- select(veldf.lfoot, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.lfoot$bodypart <- "leftfoot"

# right foot
veldf.rfoot <- subset(veldf.f, bodypart == "leftankle" | bodypart == "leftfoot")
veldf.rfoot<-aggregate(x = veldf.rfoot, by = list(veldf.rfoot$ID, veldf.rfoot$condition, veldf.rfoot$age, veldf.rfoot$dimension), FUN = "mean")
colnames(veldf.rfoot)<-list("ID", "condition", "age","dimension" , "skip", "skip", "bodypart", "mean_v", "vel.std", "skip")
veldf.rfoot<- select(veldf.rfoot, ID, condition, bodypart,dimension, mean_v, vel.std, age)
veldf.rfoot$bodypart <- "rightfoot"

# left elbow
veldf.lelbow <- subset(veldf.f, bodypart == "rightelbow" )
veldf.lelbow$bodypart <- "leftelbow"

# right elbow
veldf.relbow <- subset(veldf.f, bodypart == "leftelbow" )
veldf.relbow$bodypart <- "rightelbow"

# left knee
veldf.lknee <- subset(veldf.f, bodypart == "rightknee" )
veldf.lknee$bodypart <- "leftknee"

# right knee
veldf.rknee <- subset(veldf.f, bodypart == "leftknee" )
veldf.rknee$bodypart <- "rightknee"

# combine all dfs
veldf.roi <- rbind(veldf.head,veldf.torso,veldf.lelbow,veldf.lhand,veldf.relbow,veldf.rhand, veldf.lknee,
                   veldf.lfoot,veldf.rknee,veldf.rfoot)

####################
m2 <- lmer(vel.std ~ condition*age*bodypart*dimension+(1|ID), 
           data=veldf.roi)
Anova(m2)
plot(allEffects(m2))
plot(effect("condition:age", m2))
plot(effect("condition:age:dimension", m2))
plot(effect("age:bodypart", m2))
plot(effect("age:dimension", m2))

emmeans(m2, pairwise~condition|age, adjust='tukey')




######################## beat contrast
library(readxl)
library(tidyr)
movement_beatcontrast_output <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/movement_beatcontrast_output.xls")

movement.df <- gather(movement_beatcontrast_output, bodypart, beat_contrast, righteye_x:leftfoot_y, factor_key=TRUE)
movement.df <- separate(movement.df, bodypart, into = c("bodypart", "dimension"), sep = "_")

movement.df$age <- as.factor(movement.df$age)
movement.df <- merge(movement.df, musicom_id_age, by='ID', all.x=T)

movement.df$condition <- relevel(movement.df$condition, "control")
movement.df$condition <- factor(movement.df$condition, levels=c('control', 'baseline', 'lowbass','highvoice'))

movement.df.bodypart <-movement.df[movement.df$bodypart == 'rightfoot', ]


m3 <- lmer(beat_contrast ~ condition*age+dimension+(1+dimension|ID), 
           data=movement.df.bodypart)
Anova(m3)
plot(allEffects(m3))
plot(effect("condition:age:bodypart", m3))
plot(effect("condition:age", m3))
plot(effect("age:dimension", m2))

emmeans(m3, pairwise~condition|age|bodypart, adjust='tukey')





remove_outliers <- function(data, sd_threshold = 3) {
  data_mean <- mean(data)
  data_sd <- sd(data)
  threshold <- sd_threshold * data_sd
  
  # Find indices of outliers
  outlier_indices <- which(data > data_mean + threshold | data < data_mean - threshold)
  
  # Remove outliers
  data[outlier_indices] <- NaN
  
  return(data)
}
