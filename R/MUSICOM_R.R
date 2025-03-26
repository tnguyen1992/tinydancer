library(readr)
library(tidyr)
library(lme4)
library(effects)
library(car)
library(ggplot2)
library(readxl)
library(emmeans)
library(dplyr)
 
setwd("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R")
ERP_output <- read_excel("ERP_all_peak.xls")
data_long_pc <- gather(ERP_output, electrode, amplitude, FP2:PO9, factor_key=TRUE)
data_long_pc

data_long_pc$age <- factor(data_long_pc$age,levels = c("3", "6" , "12"),
                              labels = c("3m", "6m", "12m"))

df_pc<- data_long_pc[data_long_pc$electrode %in% c('FCz', 'FC3', 'FC4', 'FC7','FC8',
                                                   'Fz', 'F3', 'F4','F7','F8',
                                                   'Cz', 'C3', 'C4'), ]

df_pc <- df_pc %>%
  group_by(electrode) %>%
  mutate(amplitude.z = scale(amplitude), center=F)

df_pc_mc <- df_pc[df_pc$condition=="baseline" | df_pc$condition=="control",]
df_pc_fr <- df_pc[df_pc$condition=="highvoice" | df_pc$condition=="lowbass",]

hist(df_pc_mc$amplitude)
library(datawizard)
df_pc_mc$amplitude.w <- winsorize(df_pc_mc$amplitude, method="zscore", threshold=3)
df_pc_fr$amplitude.w <- winsorize(df_pc_fr$amplitude, method="zscore", threshold=3)



m1 <- lmer(amplitude ~ condition * age + electrode + (1+condition|ID), 
           data=df_pc_mc, REML = T, na.action=na.omit)
Anova(m1,2)
plot(effect('condition:age',m1))
plot(effect('age',m1))

emmeans(m1, pairwise~condition|age)
emmeans(m1, pairwise~age)
# emmeans(m1, pairwise~age|condition)

# emmeans(m1, pairwise~condition)





df <- aggregate(amplitude.z ~ age + ID + condition, data=df_pc_mc, FUN=sum)
m2 <- lmer(V1 ~ condition * age  + (1|ID), 
           data=df, REML = F, na.action=na.omit)
Anova(m2,3)
plot(effect('condition:age:electrode',m1))
plot(effect('condition',m1))

emmeans(m1, pairwise~condition|age|electrode)



df.baseline <- df[df$condition=="baseline",]
df.control <- df[df$condition=="control",]

df.mvsc <- merge(df.baseline,df.control,by=c("ID", "age"))
df.mvsc$amplitude <- df.mvsc$amplitude.x - df.mvsc$amplitude.y

df.mvsc <- df.mvsc[,c(1,2,7)]

####################################
setwd("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM")
ds_final <- read_delim("KDE_enh_8beat_pm_gaus_bwf.csv", delim = ";",
                       escape_double = FALSE, col_types = cols(...1 = col_skip()),
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)




ds_final$Child <- as.factor(ds_final$Child)
ds_final$cond <- as.factor(ds_final$cond)
ds_final$PM <- as.factor(ds_final$PM)
ds_final$age <- as.factor(ds_final$age)

ds_final$lag <- as.numeric(round(ds_final$X_round*1000))

ds_final$Child <- factor(ds_final$Child,levels = c("1","2","3","4","5","6","7",
                                                   "8","9","10","11","12","13",
                                                   "14","15","16","17","18","19",
                                                   "20","21","22","23","24",
                                                   "101","102","103","104","105","106","107",
                                                   "108","109","110","111","112","113",
                                                   "114","115","116","117","118","119",
                                                   "120","121","122","123","124","125",
                                                   "201","202","203","204","205","206","207",
                                                   "208","209","210","211","212","213",
                                                   "214","215","216","217","218","219",
                                                   "220","221","222","223","224"),
                         labels = c("13", "19","23","30","31","36","39",
                                    "43","46","48","54","60","61",
                                    "62","63","64","65","66","68",
                                    "72","74","75","76","77",
                                    "4","9","11","34","38","40","47",
                                    "49","52","78","79","81","83",
                                    "84","85","86","87","90","91",
                                    "92","93","94","95","96","97",
                                    "1","6","7","10","12","18","22",
                                    "24","26","27","28","33","35",
                                    "41","20","55","57","58","59",
                                    "67","69","71","73","80"))



pm_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)


library(permutes)

df_music <- ds_final[ds_final$cond == "Music", ]
df_control <- ds_final[ds_final$cond == "Control", ]

df_mvsc <- merge(df_music,df_control, by=c("Child","age", "X_round", "PM"))
df_mvsc$Y <- df_mvsc$Y.x - df_mvsc$Y.y
df_mvsc <-  df_mvsc[,c(1,2,3,4,11)]
lme_df <- df_mvsc[df_mvsc$X_round>=0.400 & df_mvsc$X_round<=0.488,]

lme_df <- aggregate(Y ~ age + Child + PM, data=lme_df, FUN=mean)
lme_df$Child <- as.character(lme_df$Child)

# xdata <- merge(lme_df,df.mvsc, by.x=c("Child"), by.y=c("ID"))

## testing on power 
power_music$power.c <- power_music$power - power_control$power 

xdata <- merge(lme_df,power_music, by.x=c("Child"), by.y=c("ID"), all.x=T)







# Winsorize the data
# xdata$amplitude.w <- pmin(pmax(xdata$amplitude,
#                              quantile(xdata$amplitude, lower_percentile)),
#                         quantile(xdata$amplitude, upper_percentile))
xdata$Y.w <- pmin(pmax(xdata$Y,
                               quantile(xdata$Y, lower_percentile)),
                          quantile(xdata$Y, upper_percentile))




# # plot 
# scatterplot(formula=xdata$amplitude ~ xdata$Y | xdata$age.x , xlab="EEG amplitude", ylab="Movement Periodicity",smooth=FALSE)
# scatterplot(formula=xdata$amplitude.w ~ xdata$Y.w | xdata$age.x, xlab="EEG latency", ylab="Movement Periodicity",smooth=T)
library(ggplot2)
library(ggpubr)

xdata$Y.w.z <- scale(xdata$Y.w,center=F)
xdata$power.c.z <- scale(xdata$power.c,center=F)
xdata$age <- xdata$age.x
xdata$Y.z <- scale(xdata$Y,center=F)

# Create a scatterplot
g1 <- ggplot(xdata[xdata$age.x=="12m",], aes(power.c.z, Y.w.z)) +
  geom_point(col="#d31f11") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot 12m", x = "SSEP", y = "Periodicity")+
  theme_minimal()+
  facet_grid(age.x~PM)+
  ylim(-2,3)
g2 <-ggplot(xdata[xdata$age.x=="3m",], aes(power.c.z, Y.w.z)) +
  geom_point(col="#50ad9f") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot 3m", x = "SSEP", y = "Periodicity")+
  theme_minimal()+
  facet_grid(~age.x)+
  ylim(-2,3)
g3 <-ggplot(xdata[xdata$age.x=="6m",], aes(power.c.z, Y.w.z)) +
  geom_point(col="#8cc5e3") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot 6m", x = "SSEP", y = "Periodicity")+
  theme_minimal()+
  facet_grid(~age.x)+
  ylim(-2,3)
ggarrange(g2,g3,g1, ncol=1)


# median_value <- median(xdata$Y)
# xdata$Y_group <- ifelse(xdata$Y > median_value, "Above Median", "Below Median")
# 
# median_value <- median(xdata$amplitude)
# xdata$amplitude_group <- ifelse(xdata$amplitude > median(xdata$amplitude), "Above Median", "Below Median")

# library(dplyr)
# xdata <- xdata %>%
#   group_by(Child) %>%
#   mutate(amplitude_group = ifelse(amplitude > median(amplitude), "Above Median", "Below Median"),
#          Y_group = ifelse(Y > median(Y), "Above Median", "Below Median"))
# 




xdata.12m <- xdata[xdata$age.x=="12m",]


m1 <- lmer(Y.w.z ~ power.c.z * PM + (1|Child), #+ I(amplitude.w^2)
         data=xdata.12m, na.action=na.omit)
summary(m1)
Anova(m1,2)

plot(effect('power.c.z',m1))
emtrends(m1, pairwise~age.y,var="power.c.z" )

contingency_table <- table(xdata$amplitude_group, xdata$Y_group)

# Chi-squared test
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)


####################################
setwd("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM")
ds_final <- read_delim("KDEfreq_enh_8beat_pm_gaus_bwf.csv", delim = ";",
                       escape_double = FALSE, col_types = cols(...1 = col_skip()),
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)




ds_final$Child <- as.factor(ds_final$Child)
ds_final$cond <- as.factor(ds_final$cond)
ds_final$PM <- as.factor(ds_final$PM)
ds_final$age <- as.factor(ds_final$age)

ds_final$lag <- as.numeric(round(ds_final$X_round*1000))

ds_final$Child <- factor(ds_final$Child,levels = c("1","2","3","4","5","6","7",
                                                   "8","9","10","11","12","13",
                                                   "14","15","16","17","18","19",
                                                   "20","21","22","23","24",
                                                   "101","102","103","104","105","106","107",
                                                   "108","109","110","111","112","113",
                                                   "114","115","116","117","118","119",
                                                   "120","121","122","123","124","125",
                                                   "201","202","203","204","205","206","207",
                                                   "208","209","210","211","212","213",
                                                   "214","215","216","217","218","219",
                                                   "220","221","222","223","224"),
                         labels = c("13", "19","23","30","31","36","39",
                                    "43","46","48","54","60","61",
                                    "62","63","64","65","66","68",
                                    "72","74","75","76","77",
                                    "4","9","11","34","38","40","47",
                                    "49","52","78","79","81","83",
                                    "84","85","86","87","90","91",
                                    "92","93","94","95","96","97",
                                    "1","6","7","10","12","18","22",
                                    "24","26","27","28","33","35",
                                    "41","20","55","57","58","59",
                                    "67","69","71","73","80"))



pm_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

library(permutes)

df_hv <- ds_final[ds_final$cond == "HighVoice", ]
df_lb <- ds_final[ds_final$cond == "LowBass", ]


hv_df <- df_hv[df_hv$X_round>=0.400 & df_hv$X_round<=0.488,]
lb_df <- df_lb[df_hv$X_round>=0.400 & df_lb$X_round<=0.488,]

hv_df <- aggregate(Y ~ age + Child + PM, data=hv_df, FUN=mean)
lb_df <- aggregate(Y ~ age + Child + PM, data=lb_df, FUN=mean)


# for music vs. control
# lme_df <- lme_df[lme_df$PM!=10 & lme_df$PM!=3 & lme_df$PM!=6 ,]
# hv_df <- aggregate(Y ~ age + Child, data=hv_df, FUN=mean)
# lb_df <- aggregate(Y ~ age + Child, data=lb_df, FUN=mean)

xdata <- merge(hv_df,df[df$condition=="highvoice",], by.x=c("Child"), by.y=c("ID"))
xdata <- merge(lb_df,df[df$condition=="lowbass",], by.x=c("Child"), by.y=c("ID"))


# Create a scatterplot
g1 <- ggplot(xdata[xdata$age.x=="12m",], aes(amplitude, Y)) +
  geom_point(col="#d31f11") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot LP 12m", x = "EEG auc", y = "Periodicity")+
  theme_minimal()+
  facet_grid(age.x~PM)
g2 <-ggplot(xdata[xdata$age.x=="3m",], aes(amplitude, Y)) +
  geom_point(col="#50ad9f") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot LP 3m", x = "EEG auc", y = "Periodicity")+
  theme_minimal()+
  facet_grid(age.x~PM)
g3 <-ggplot(xdata[xdata$age.x=="6m",], aes(amplitude, Y)) +
  geom_point(col="#8cc5e3") +
  geom_smooth(method = "lm")+
  labs(title = "Scatterplot LP 6m", x = "EEG auc", y = "Periodicity")+
  theme_minimal()+
  facet_grid(age.x~PM)
ggprubr::ggarrange(g2,g3,g1, ncol=1)




median_value <- median(xdata$Y)
xdata$Y_group <- ifelse(xdata$Y > median_value, "Above Median", "Below Median")

median_value <- median(xdata$amplitude)
xdata$amplitude_group <- ifelse(xdata$amplitude > median_value, "Above Median", "Below Median")


contingency_table <- table(xdata$amplitude_group, xdata$Y_group)

# Chi-squared test
chi_squared_test <- chisq.test(contingency_table)
print(chi_squared_test)

m1 <- lm(Y ~ amplitude * PM , 
         data=xdata[xdata$age.x=="12m",], na.action=na.omit)
confint(m1)
Anova(m1)
plot(effect('amplitude',m1))
emtrends(m1, var="amplitude" )

xdata <- merge(lmedf,df_pc[df_pc$condition=="highvoice",], by.x=c("Child"), by.y=c("ID"))

m1 <- lmer(amplitude ~ Y * age.x + electrode + (1|Child), 
           data=xdata, na.action=na.omit)
Anova(m1)
plot(effect('Y',m1))
emtrends(m1, var="Y" )



#######################
library(readxl)
power_music <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/power_stats.xlsx", 
                          col_types = c("text", "text", "text", 
                                        "numeric"), 
                          sheet = "music")
power_control <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/power_stats.xlsx", 
                          col_types = c("text", "text", "text", 
                                        "numeric"), 
                          sheet = "control")
power_lowbass <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/power_stats.xlsx", 
                          col_types = c("text", "text", "text", 
                                        "numeric"), 
                          sheet = "lowbass")
power_highvoice <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/power_stats.xlsx", 
                            col_types = c("text", "text", "text", 
                                          "numeric"), 
                            sheet = "highvoice")
 
power_all <- rbind(power_music,power_control,power_lowbass,power_highvoice)
power_all$power.z <- scale(power_all$power,center = F)[,1]

power_mus <- power_all[power_all$condition=="music" | power_all$condition=="control",]
power_mus$age <- factor(power_mus$age,levels=c("3","6","12","ad"), labels=c("3m","6m","12m","ad"))
power_mus$condition <- factor(power_mus$condition)

m1 <- lmer(power ~ condition*age + (1|ID)  , 
         data=power_mus, na.action=na.omit)
Anova(m1,3)
emmeans::emmeans(m1, pairwise ~  condition|age )
power_mus_6m <- power_mus[power_mus$age=="6m",]
m6 <- lm(power.z ~ condition  , 
           data=power_mus_6m, na.action=na.omit)
Anova(m6)

power_mus_12m <- power_mus[power_mus$age=="12m",]
m12 <- lm(power.z ~ condition  , 
         data=power_mus_12m, na.action=na.omit)
Anova(m12)

power_mus_3m <- power_mus[power_mus$age=="3m",]
m3 <- lm(power.z ~ condition  , 
          data=power_mus_3m, na.action=na.omit)
Anova(m3)

power_mus_ad <- power_mus[power_mus$age=="ad",]
ad <- lm(power.z ~ condition  , 
         data=power_mus_ad, na.action=na.omit)
Anova(ad)

p.adjust(c(0.00732,0.00106,0.09158), method="fdr")

ggplot(data=power_mus,aes(condition, power,color=age))+
  # geom_bar(stat="identity", position = position_dodge())
  geom_jitter()+
  geom_violin(aes(fill=age, alpha=0.3))+
  geom_boxplot(aes(col=age, alpha=0.5))+
  theme_minimal()+
  facet_grid(~age)+
  theme(legend.position = "none")
  

power_music <- power_all[power_all$condition=="music",]
m.mus <- lm(power.z ~1 , 
           data=power_music, na.action=na.omit)
summary(m.mus,3)

power_control <- power_all[power_all$condition=="control",]
m.con <- lm(power.z ~1 , 
            data=power_control, na.action=na.omit)
summary(m.con,3)


power_pitch <- power_all[power_all$condition=="lowbass" | power_all$condition=="highvoice",]
power_pitch$age <- factor(power_pitch$age,levels=c("3","6","12",'ad'), labels=c("3m","6m","12m",'ad'))
power_pitch$condition <- factor(power_pitch$condition)


power_pitch_6m <- power_pitch[power_pitch$age=="6m",]
m6 <- lm(power.z ~ condition  , 
         data=power_pitch_6m, na.action=na.omit)
summary(m6)

power_pitch_12m <- power_pitch[power_pitch$age=="12m",]
m12 <- lm(power.z ~ condition  , 
          data=power_pitch_12m, na.action=na.omit)
summary(m12)

power_pitch_3m <- power_pitch[power_pitch$age=="3m",]
m3 <- lm(power.z ~ condition  , 
         data=power_pitch_3m, na.action=na.omit)
summary(m3)

power_pitch_ad <- power_pitch[power_pitch$age=="ad",]
ad <- lm(power.z ~ condition  , 
         data=power_pitch_ad, na.action=na.omit)
summary(ad)

ggplot(data=power_pitch,aes(condition, power,color=age))+
  # geom_bar(stat="identity", position = position_dodge())
  geom_jitter()+
  geom_violin(aes(fill=age, alpha=0.3))+
  geom_boxplot(aes(col=age, alpha=0.5))+
  theme_minimal()+
  facet_grid(~age)+
  theme(legend.position = "none")


