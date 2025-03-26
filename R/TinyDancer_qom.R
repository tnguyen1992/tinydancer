# QoM script

library(readr)
library(dplyr)
library(pracma)

# Replace "your_folder_path" with the actual path to your folder
ages<- c('3m', "6m", "12m")
# conditions <- c("HighVoice", "LowBass")
conditions <- c("Baseline", "Control","HighVoice", "LowBass")

pm_frame_qom_all <- list()  # Assuming the third dimension is for x and y values
# pm_frame_y_all <- list()



for (age in ages) {
  if (age=="6m") {
    NB<-25
  } else {
    NB=24
  }
  # Set condition
  # Set the range of baby numbers
  baby_numbers <- sprintf("%02d", 1:NB)  # Formatting to ensure two-digit representation (e.g., 01, 02, ..., 24)
  
  # folder_path <- paste0("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/tinyDancers_commonAges/tinyDancers_commonAges/csv_for_trinh/age",age)
  folder_path <- paste0("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/tinyDancers_commonAges/tinyDancers_commonAges/csv_for_trinh/age",age)
  
  # Create an empty data frame to store results
  pm_frame_qom_list <- list()  # Assuming the third dimension is for x and y values
  # pm_frame_y_list <- list()
  
  for (condition in conditions) {
    pm_frame_qom <- array(NA, dim = c(10, NB, 1))  # Assuming the third dimension is for x and y values
    # pm_frame_y <- array(NA, dim = c(10, NB, 91))  # Assuming the third dimension is for x and y values
    
    # Loop through each baby
    for (baby_number in baby_numbers) {
      # Baseline Condition
      pattern<-paste0(condition,".*\\.csv$")
      file_names <- list.files(paste0(folder_path,"/babe",baby_number), pattern = pattern,recursive=T, full.names = TRUE)
      
      # Read, transpose, and stack files dynamically
      list_of_dfs <- lapply(file_names, function(file) {
        df <- read.csv(file)
        transposed_df <- t(df)
        transposed_df <- transposed_df[-1, ]  # Exclude the first row
        transposed_ar <- diff(transposed_df)
        # result_matrix <- apply(transposed_df, 2, function(column) ifelse(column > mean(column, na.rm = TRUE) + 2 * sd(column, na.rm = TRUE), column, 0))
      })
      velocity_data <- do.call(rbind, list_of_dfs)
      
      # Set the range of PM values
      pm_values <- 1:10
      
      # Loop over each PM value
      for (PM in pm_values) {
        velocity_PM<-velocity_data[,PM]
        qom <- mean(abs(velocity_PM), na.rm=T)
        qom.sd <- sd(abs(velocity_PM), na.rm=T)
        qom.max <- max(abs(velocity_PM), na.rm=T)
        # SD movement variance, max movement 
        # Granger causality: Audio , RQA 
        pm_frame_qom[PM,as.numeric(baby_number),] <- qom.sd
        
        # Plot the time series
        # time <- 1:length(velocity_PM)
        # plot(time, velocity_PM, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Peaks", main = "Area Under the Curve")
        
        # Highlight the area under the curve
        # polygon(c(0, time, max(time)), c(0, velocity_PM, 0), col = "lightgray", border = NA)
        
      }  
    }
    pm_frame_qom_list[[condition]] <- pm_frame_qom
  }
  pm_frame_qom_all[[age]] <- pm_frame_qom_list
}

###############
# Create an empty data frame to store the long-format data
m12_df <- data.frame()
m6_df <- data.frame()
m3_df <- data.frame()

# Loop through each condition in pm_frame_x_list
for (condition in names(pm_frame_qom_list)) {
  # Extract the data from pm_frame_x_list for the current condition
  pm_frame_qom_condition <- pm_frame_qom_all[["12m"]][[condition]]
  
  # Combine the arrays into a long-format data frame
  long_df <- data.frame(
    Infant = rep(1:24, times = 1 , each= 10),  # Assuming 24 infants, 10 PMs, and 512 values per PM
    PM = rep(1:10, times = 24),  # Assuming 10 PMs
    Value_qom = c(pm_frame_qom_condition[, , ])  # Flatten the pm_frame_x array
  )
  
  # Add a Condition column to identify the condition
  long_df$Condition <- condition
  
  # Append the long-format data for the current condition to the overall data frame
  m12_df <- rbind(m12_df, long_df)
}
for (condition in names(pm_frame_qom_list)) {
  # Extract the data from pm_frame_x_list for the current condition
  pm_frame_qom_condition <- pm_frame_qom_all[["6m"]][[condition]]
  
  # Combine the arrays into a long-format data frame
  long_df <- data.frame(
    Infant = rep(1:25, times = 1 , each= 10),  # Assuming 24 infants, 10 PMs, and 512 values per PM
    PM = rep(1:10, times = 25 * 1),  # Assuming 10 PMs
    Value_qom = c(pm_frame_qom_condition[, , ])  # Flatten the pm_frame_x array
  )
  
  # Add a Condition column to identify the condition
  long_df$Condition <- condition
  
  # Append the long-format data for the current condition to the overall data frame
  m6_df <- rbind(m6_df, long_df)
}
for (condition in names(pm_frame_qom_list)) {
  # Extract the data from pm_frame_x_list for the current condition
  pm_frame_qom_condition <- pm_frame_qom_all[["3m"]][[condition]]
  
  # Combine the arrays into a long-format data frame
  long_df <- data.frame(
    Infant = rep(1:24, times = 1 , each= 10),  # Assuming 24 infants, 10 PMs, and 512 values per PM
    PM = rep(1:10, times = 24 * 1),  # Assuming 10 PMs
    Value_qom = c(pm_frame_qom_condition[, , ])  # Flatten the pm_frame_x array
  )
  
  # Add a Condition column to identify the condition
  long_df$Condition <- condition
  
  # Append the long-format data for the current condition to the overall data frame
  m3_df <- rbind(m3_df, long_df)
}

m12_df$age <- "12m"
m6_df$age <- "6m"

m6_df$Infant <- m6_df$Infant+100
m12_df$Infant <- m12_df$Infant+200

m3_df$age <- "3m"

long_format_df <- rbind(m12_df,m6_df,m3_df)

library(lme4)
library(car)
library(emmeans)
library(permutes)
library(effects)
long_format_df$Infant <- as.factor(long_format_df$Infant)
long_format_df$Condition <- factor(long_format_df$Condition, levels=c("Control","Baseline", "HighVoice", "LowBass"),
                                   labels=c("Control","Music","High Pitch", "Low Pitch"))
long_format_df$PM <- as.factor(long_format_df$PM)
long_format_df$age <- as.factor(long_format_df$age)

# sum over PMs
result <- aggregate(Value_qom ~ Infant+PM+age, data = long_format_df, FUN = mean)

xdata<-long_format_df[long_format_df$Condition=="Music" | long_format_df$Condition=="Control" ,]
# xdata<-long_format_df[long_format_df$Condition=="High Pitch" | long_format_df$Condition=="Low Pitch" ,]

hist(xdata$Value_qom, breaks = 200)
max(xdata$Value_qom)
mean(xdata$Value_qom)
sd(xdata$Value_qom)

outlier <- mean(xdata$Value_qom) + 5*sd(xdata$Value_qom)
xdata$Value_qom[xdata$Value_qom>outlier]<-NA
hist(xdata$Value_qom, breaks = 200)


ds_final_check <- xdata[!xdata$Infant==20,]
ds_final_check <- ds_final_check[!ds_final_check$Infant==103,]
ds_final_check <- ds_final_check[!ds_final_check$Infant==211,]
xdata <- ds_final_check[!ds_final_check$Infant==214,]



m1 <- lmer(Value_qom ~ Condition*age*PM+ (1+Condition|Infant), 
           data=xdata)
Anova(m1,2)
summary(m1)

plot(allEffects(m1))
plot(effect("Condition:age", m1))
emmeans(m1, pairwise ~ Condition|age|PM, adjust="none")
emmeans(m1, pairwise ~ age|PM, adjust="fdr")
emmeans(m1, pairwise ~ Condition|age, adjust="fdr")


# Define the desired order of factor levels
desired_order <- c("3m", "6m", "12m")

# Reorder the factor levels
xdata$age.n <- ordered(xdata$age, levels = desired_order)
library(plyr)
xdata_summary <- ddply(xdata,.(PM,Condition,age.n),summarize,avg_value=mean(Value_qom,na.rm=TRUE),se_value=sd(Value_qom,na.rm=TRUE)/length(Value_qom))
xdata_summary_superPM <- ddply(xdata,.(Condition,age.n),summarize,avg_value=mean(Value_qom,na.rm=TRUE),se_value=sd(Value_qom,na.rm=TRUE)/length(Value_qom))

library(ggplot2)
ggplot(data=xdata_summary, aes(x=age.n, y=avg_value,fill=Condition,group=Condition))+#, group=interaction(PM,Condition)))+
  geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
  geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
                position="dodge",width=0.5, linewidth=0.5)+
  # geom_dotplot()+
  facet_grid(PM~"")+
  theme_minimal()+
  theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
  ylab("")+xlab("")+ylim(0,420)

ggplot(data=xdata_summary_superPM, aes(x=age.n, y=avg_value,fill=Condition,group=Condition))+#, group=interaction(PM,Condition)))+
  geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
  geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
                position="dodge",width=0.5, linewidth=0.5)+
  # geom_dotplot()+
  # facet_grid(PM~"")+
  theme_minimal()+
  theme(legend.position = "none",panel.grid = element_blank(),strip.text=element_blank())+
  ylab("")+xlab("")+ylim(0,400)




#### load ERP data
ERP_all_auc <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/ERP_all_auc.xls")
library(dplyr)

# Rename condition variable
ERP_all_auc<- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "baseline", "Music", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "control", "Control", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "lowbass", "Low Pitch", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "highvoice", "High Pitch", condition))

# Rename age variable
ERP_all_auc <- ERP_all_auc %>%
  mutate(age = case_when(
    age == 3 ~ "3m",
    age == 6 ~ "6m",
    age == 12 ~ "12m"
  ))



# Renaming other variables, adjust as needed
names(ERP_all_auc)[names(ERP_all_auc) == "ID"] <- "Child"
names(ERP_all_auc)[names(ERP_all_auc) == "condition"] <- "cond"

library(tidyr)

data_long_pc <- gather(ERP_all_auc, electrode, amplitude, FP2:PO9, factor_key=TRUE)
df_pc<- data_long_pc[data_long_pc$electrode %in% c('FCz', 'FC3', 'FC4', 'FC7','FC8',
                                                   'Fz', 'F3', 'F4','F7','F8',
                                                   'Cz', 'C3', 'C4'), ]

df_pc <- df_pc %>%
  group_by(electrode) %>%
  mutate(amplitude.z = scale(amplitude)[,1])

df <- aggregate(amplitude ~ age + Child + cond, data=df_pc, FUN=mean)


df$cond <- as.factor(df$cond)
df$age <- as.factor(df$age)
df$Child <- as.factor(df$Child)


long_format_df$cond <- as.factor(long_format_df$Condition)
long_format_df$age <- as.factor(long_format_df$age)
long_format_df$Child <- as.factor(long_format_df$Infant)


long_format_df$Child <- factor(long_format_df$Child,levels = c("1","2","3","4","5","6","7",
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


em.data <- merge(long_format_df,df, by=c("Child","cond", "age"), all=T)
qe.data<-em.data[em.data$Condition=="Music" | em.data$Condition=="Control" ,]
qe.data<-em.data[em.data$Condition=="High Pitch" | em.data$Condition=="Low Pitch" ,]

qe.m1 <- lmer(Value_qom ~ Condition*age*amplitude+ (1|Child), 
           data=qe.data)
Anova(qe.m1,3)
summary(qe.m1)
confint(qe.m1)

plot(allEffects(qe.m1))
plot(effect("amplitude", qe.m1))
emmeans(m1, pairwise ~ Condition|age|PM, adjust="none")
emmeans(m1, pairwise ~ age|PM, adjust="fdr")
emmeans(m1, pairwise ~ Condition, adjust="fdr")


## load assr for assr qom analyses

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

power_mus <- power_all
power_mus$age <- factor(power_mus$age,levels=c("3","6","12","ad"), labels=c("3m","6m","12m","ad"))
power_mus$condition <- factor(power_mus$condition,levels=c("music","control", "highvoice", "lowbass"), labels=c("Music","Control", "High Pitch", "Low Pitch"))

power_mus$cond <- as.factor(power_mus$condition)
power_mus$age <- as.factor(power_mus$age)
power_mus$Child <- as.factor(power_mus$ID)

power_mus.data <- merge(long_format_df,power_mus, by=c("Child","cond", "age"), all=T)
pe.data<-power_mus.data[power_mus.data$Condition=="Music" | power_mus.data$Condition=="Control" ,]
pe.data<-power_mus.data[power_mus.data$Condition=="High Pitch" | power_mus.data$Condition=="Low Pitch" ,]

pe.m1 <- lmer(Value_qom ~ cond*age*power.z+PM+ (1|Child), 
              data=pe.data)
Anova(pe.m1,2)
summary(pe.m1)
confint(pe.m1)

plot(allEffects(pe.m1))
plot(effect("cond:power.z", pe.m1))
emtrends(pe.m1, pairwise ~ cond, var='power.z',adjust="fdr")


pe.data<-em.data[em.data$Condition=="High Pitch" | em.data$Condition=="Low Pitch" ,]

