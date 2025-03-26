# Auto-correlation script

library(readr)
library(dplyr)
library(pracma)
library(lmtest)
library(signal)

env_hopp_lp <- read_csv("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_hopp_lp.csv", 
                        col_names = FALSE)
env_hopp_hp <- read_csv("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_hopp_hp.csv",
                        col_names = FALSE)
env_lola_lp <- read_csv("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_lola_lp.csv", 
                        col_names = FALSE)
env_lola_hp <- read_csv("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_lola_hp.csv",
                        col_names = FALSE)



# Replace "your_folder_path" with the actual path to your folder
ages<- c('3m', "6m", "12m")
conditions <- c("HLowBass", "HHighVoice","LLowBass", "LHighVoice")
lags <- c(1,seq(11,165,by=11),seq(6,165,by=11))
pm_frame_x_all <- list()  # Assuming the third dimension is for x and y values

for (age in ages) {
  if (age == "6m") {
    NB <- 25
  } else {
    NB = 24
  }
  # Set condition
  # Set the range of baby numbers
  baby_numbers <- sprintf("%02d", 1:NB)  # Formatting to ensure two-digit representation (e.g., 01, 02, ..., 24)
  folder_path <- paste0("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/tinyDancers_commonAges/tinyDancers_commonAges/csv_for_trinh/age", age)
  
  # Create an empty data frame to store results
  pm_frame_x_list <- list()  # Assuming the third dimension is for x and y values
  
  for (condition in conditions) {
    granger_F_values <- array(NA, dim = c(10, NB,length(lags)))
    
    if (condition=="HLowBass") {
      stim=env_hopp_lp
    } 
    if (condition=="HHighVoice"){
      stim=env_hopp_hp
    } 
    if (condition=="LLowBass"){
      stim=env_lola_lp
    } 
    if (condition=="LHighVoice") {
      stim=env_lola_hp
    }
    
    
    # Loop through each baby
    for (baby_number in baby_numbers) {
      pattern <- paste0(condition, ".*\\.csv$")
      file_names <- list.files(paste0(folder_path, "/babe", baby_number), pattern = pattern, recursive = T, full.names = TRUE)
      collectF <- array(NA, dim = c(length(file_names),length(lags)))
      
      # Set the range of PM values
      pm_values <- 1:10
      # Loop over each PM value
      for (PM in pm_values) {
        
        for (i in seq_along(file_names)) {
          df <- read.csv(file_names[i])
          transposed_df <- t(df)
          transposed_df <- transposed_df[-1, ]  # Exclude the first row
          transposed_vel <- diff(transposed_df)
          
          velocity_PM <- transposed_vel[, PM]
          
          for (l in seq_along(lags)) {
            # Perform Granger causality test
            # gctest <- suppressWarnings(grangertest(stim, velocity_PM, order = lags[l], na.action = na.omit))
            gctest <- suppressWarnings(grangertest(velocity_PM, stim, order = lags[l], na.action = na.omit))
            
            collectF[i,l] <- gctest$F[2]
          }
        }
        granger_F_values[PM,as.numeric(baby_number),] <- colMeans(collectF,na.rm=T)
        
      }
    }
    pm_frame_x_list[[condition]] <- granger_F_values
    
  }
  pm_frame_x_all[[age]] <- pm_frame_x_list
  
}


pm_frame_x_hcon_3m <-pm_frame_x_all[["3m"]][["HHighVoice"]]
pm_frame_x_lcon_3m <-pm_frame_x_all[["3m"]][["LHighVoice"]]

# Create a new 10 x 24 array by calculating the mean of corresponding values
pm_frame_x_con_3m <- array(NA, dim=c(10,24,length(lags)))

for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_con_3m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_3m[i, j,], pm_frame_x_lcon_3m[i, j,]))
  }
}


pm_frame_x_hmus_3m <-pm_frame_x_all[["3m"]][["HLowBass"]]
pm_frame_x_lmus_3m <-pm_frame_x_all[["3m"]][["LLowBass"]]

pm_frame_x_mus_3m <- array(NA, dim=c(10,24,length(lags)))

for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_mus_3m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_3m[i, j,], pm_frame_x_lmus_3m[i, j,]))
  }
}

pm_frame_x_hcon_6m <-pm_frame_x_all[["6m"]][["HHighVoice"]]
pm_frame_x_lcon_6m <-pm_frame_x_all[["6m"]][["LHighVoice"]]
pm_frame_x_con_6m <- array(NA, dim=c(10,25,length(lags)))
for (i in 1:10) {
  for (j in 1:25) {
    pm_frame_x_con_6m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_6m[i, j,], pm_frame_x_lcon_6m[i, j,]))
  }
}

pm_frame_x_hmus_6m <-pm_frame_x_all[["6m"]][["HLowBass"]]
pm_frame_x_lmus_6m <-pm_frame_x_all[["6m"]][["LLowBass"]]
pm_frame_x_mus_6m <- array(NA, dim=c(10,25,length(lags)))
for (i in 1:10) {
  for (j in 1:25) {
    pm_frame_x_mus_6m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_6m[i, j,], pm_frame_x_lmus_6m[i, j,]))
  }
}

pm_frame_x_hcon_12m <-pm_frame_x_all[["12m"]][["HHighVoice"]]
pm_frame_x_lcon_12m <-pm_frame_x_all[["12m"]][["LHighVoice"]]
pm_frame_x_con_12m <- array(NA, dim=c(10,24,length(lags)))
for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_con_12m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_12m[i, j,], pm_frame_x_lcon_12m[i, j,]))
  }
}

pm_frame_x_hmus_12m <-pm_frame_x_all[["12m"]][["HLowBass"]]
pm_frame_x_lmus_12m <-pm_frame_x_all[["12m"]][["LLowBass"]]
pm_frame_x_mus_12m <- array(NA, dim=c(10,24,length(lags)))
for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_mus_12m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_12m[i, j,], pm_frame_x_lmus_12m[i, j,]))
  }
}


##################################
# extracting 3m old info
# extracting 3m old info
combinations <- expand.grid(
  PM = 1:10,
  Child = 1:24,
  Lags = 1:length(lags)
)

# Create a data frame to store the results
result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),
  X = numeric(0)
)

# Loop through combinations and fill the data frame
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_con_3m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
# Print the resulting data frame
con_3m_df <- result_df

result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),
  X = numeric(0)
)
# Loop through combinations and fill the data frame
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_mus_3m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
mus_3m_df <- result_df

con_3m_df$cond <- 'HighPitch'
mus_3m_df$cond <- 'LowPitch'

cond_3m_df <- rbind(mus_3m_df,con_3m_df)


# extracting 6m old info
combinations <- expand.grid(
  PM = 1:10,
  Child = 1:25,
  Lags = 1:length(lags)
)

# Create a data frame to store the results
result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),
  X = numeric(0)
)

# Loop through combinations and fill the data frame
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_con_6m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
# Print the resulting data frame
con_6m_df <- result_df

result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),  
  X = numeric(0)
)
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_mus_6m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
mus_6m_df <- result_df

con_6m_df$cond <- 'HighPitch'
mus_6m_df$cond <- 'LowPitch'

cond_6m_df <- rbind(mus_6m_df,con_6m_df)


# now 12 m olds
combinations <- expand.grid(
  PM = 1:10,
  Child = 1:24,
  Lags = 1:length(lags)
)

# Create a data frame to store the results
result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),
  X = numeric(0)
)


# Loop through combinations and fill the data frame
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_con_12m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
# Print the resulting data frame
con_12m_df <- result_df

result_df <- data.frame(
  PM = numeric(0),
  Child = character(0),
  Lags = numeric(0),
  X = numeric(0)
)
# Loop through combinations and fill the data frame
for (i in 1:nrow(combinations)) {
  pm_idx <- combinations$PM[i]
  child_idx <- combinations$Child[i]
  Lags_idx <- combinations$Lags[i]
  x_val <- pm_frame_x_mus_12m[pm_idx,child_idx,Lags_idx]
  result_df <- rbind(result_df, data.frame(
    PM = pm_idx,
    Child = child_idx,
    Lags = lags[Lags_idx]/25,
    X = x_val
  ))
}
mus_12m_df <- result_df

con_12m_df$cond <- 'HighPitch'
mus_12m_df$cond <- 'LowPitch'

cond_12m_df <- rbind(mus_12m_df,con_12m_df)

cond_3m_df$age <- '3m'
cond_6m_df$age <- '6m'
cond_12m_df$age <- '12m'

cond_6m_df$Child <- cond_6m_df$Child+100
cond_12m_df$Child <- cond_12m_df$Child+200

cond_all_df <- rbind(cond_3m_df,cond_6m_df,cond_12m_df)


write.csv2(cond_all_df,"granger_freq_F_range_minuslags_upsample.csv")
# write.csv2(cond_all_df,"granger_freq_F_range_upsample.csv")


##########################################

pm_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
age_list <- c("3m","6m","12m")
library(readr )

library(ggplot2)
library(gridExtra )
library(glmmTMB )
library(car )
library(effects)
library(emmeans)
library(DescTools)
library(dplyr)

# ds_final <- cond_all_df

granger_freq_F_range <- 
  read_delim("granger_freq_F_range_upsample_110324.csv", 
                                  delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                  locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                  trim_ws = TRUE)
granger_freq_F_range$Lags <- as.numeric(granger_freq_F_range$Lags)
granger_freq_F_range$direction <- "music2move"

granger_freq_F_range_minuslags <- read_delim("granger_freq_F_minusrange_upsample_110324.csv", 
                                            delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                            locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                            trim_ws = TRUE)

granger_freq_F_range_minuslags$Lags <- as.numeric(granger_freq_F_range_minuslags$Lags)
granger_freq_F_range_minuslags$direction <- "move2music"

ds_final <- rbind(granger_freq_F_range_minuslags,granger_freq_F_range)


ds_final$PM <- as.factor(ds_final$PM)
# ds_final$X.w <- Winsorize(ds_final$X,na.rm=T)
ds_final$age <- factor(ds_final$age, levels=c('3m', '6m', '12m'))
ds_final$Lags.f <- as.factor(ds_final$Lags)#, levels = c('0.04','0.44','0.88',
                                                  # '1.32','1.76','2.20',
                                                  # '2.64','3.08','3.52',
                                                  # '3.96','4.40','4.84',
                                                  # '5.28','5.72','6.16',
                                                  # '6.60'))

# ########### ttest
# # Assuming your data frame is named ds_final
# 
# # Define the age groups, PMs, and conditions
# age_groups <- c("3m", "6m", "12m")
# PMs <- 1:10
# 
# # Create an empty list to store t-test results
# t_test_results <- list()
# 
# library(MKinfer)
# 
# # Loop through each combination of age group, PM, and condition
# for (age in age_groups) {
#   for (PM in PMs) {
#     # Subset the data for the current combination
#     
#     data1 <- ds_final[ds_final$age == age & ds_final$PM == PM & ds_final$cond == "Music", ]
#     data2 <- ds_final[ds_final$age == age & ds_final$PM == PM & ds_final$cond == "Control", ]
#     
#     # Perform paired t-test if both conditions have at least 2 observations
#     t_test_result <- boot.t.test(data1$X, data2$X, alternative = "greater", paired=T)
#     
#     # Store the result
#     t_test_results[[paste(age, "PM", PM)]] <- t_test_result$p.value
#     
#   }
#   
# }
# 
# # Accessing the results
# # For example, to access the result for age 3m, PM 1, Music condition
# print(t_test_results)
# 
# values_below_005 <- t_test_results[t_test_results < 0.05]


#### LMM
library(lme4)

ds_final$Category <- rep(NA, nrow(ds_final))

# Identify Upper body, Lower body, and Hybrid PMs
upper_body_pms <- c(1, 2, 5, 6)
lower_body_pms <- c(4, 7, 9, 10)
hybrid_pms <- c(3, 8)

# Assign categories based on PM numbers
ds_final$Category[ds_final$PM %in% upper_body_pms] <- "Upper body"
ds_final$Category[ds_final$PM %in% lower_body_pms] <- "Lower body"
ds_final$Category[ds_final$PM %in% hybrid_pms] <- "Hybrid"

# Convert Category to a factor
ds_final$Category <- factor(ds_final$Category, levels = c("Upper body", "Lower body", "Hybrid"))


library(MKinfer)
library(dplyr)
result <- ds_final %>%
  group_by(Lags.f, age,direction) %>%
  summarize(t_test_p_value = boot.t.test(X ~ cond, alternative="two.sided")$p.value)

# Display the result
print(result)


result <- result %>%
  mutate(sig_line = ifelse(t_test_p_value < 0.005, 0.5, NA))
plotdf <- aggregate(X ~ cond * age * Lags.f  * direction,  data=ds_final, FUN="mean", na.rm=T)

plotdf.se <- aggregate(X ~ cond * age * Lags.f  * direction,  data=ds_final, FUN="SD", na.rm=T)
colnames(plotdf.se) <- list("cond", "age", "Lags.f","direction","SD")

plotdf.f <- merge(plotdf,plotdf.se, by=c("age", "Lags.f","direction", "cond"),all.X=T)
plotdf.f <- merge(plotdf.f,result, by=c("age",  "Lags.f", "direction"))


# ggplot(data = plotdf.f[plotdf.f$cond=="HighPitch",], aes(x = Lags.f, y = X, group=interaction(direction,age), col=direction)) +
#   # geom_point()+
#   geom_line() +
#   geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
#   geom_vline(xintercept=0,linetype="dotted")+
#   # geom_vline(xintercept="0.44",linetype="dotted")+
#   # geom_vline(xintercept="0.88",linetype="dotted")+
#   # geom_vline(xintercept="1.76",linetype="dotted")+
#   # geom_vline(xintercept="3.52",linetype="dotted")+
#   # geom_line(aes(y = sig_line),
#   #           color = "red", size = 1) +
#   labs(title = "", x = "Lags in seconds", y = "Granger Causality F-value")+
#   theme_minimal()+
#   theme(legend.position = "top")+
#   facet_grid(age~"")

plotdf.f$Lags <- as.numeric(plotdf.f$Lags.f)
ggplot(data = plotdf.f[plotdf.f$direction=="music2move" & plotdf.f$Lags<15 ,], aes(x = Lags, y = X, group=interaction(cond,age), col=cond)) +
  # geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
  # geom_vline(xintercept=0,linetype="dotted")+
  # geom_vline(xintercept="0.44",linetype="dotted")+
  # geom_vline(xintercept="0.88",linetype="dotted")+
  # geom_vline(xintercept="1.76",linetype="dotted")+
  # geom_vline(xintercept="3.52",linetype="dotted")+
  geom_line(aes(y = sig_line),
            color = "black", size = 1) +
  labs(title = "", x = "Lags in frames (sr: 25 Hz)", y = "Granger Causality F-value")+
  theme_bw()+
  theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
  facet_grid(~age)









m1 <- lmer(X ~  cond * age * Category  + (1+cond|Child), 
           data=ds_final[ds_final$Lags.f==1.76,])
m1 <- lmer(X ~  cond * age * Category  + (1+cond|Child), 
           data=ds_final[ds_final$Lags.f==0.24,])
Anova(m1,2)
summary(m1)
emmeans(m1, pairwise ~ cond|age|Category)

# take ds.lag from music vs control as localizer
# library(dplyr)
# ds_lag <- ds_final[complete.cases(ds_final), ]
# ds_lag <- ds_lag[ds_lag$direction == "music2move",]
# ds.lag <- ds_lag %>%
#   group_by(Child, PM) %>%
#   slice_max(order_by = X, n = 1) %>%
#   ungroup() %>%
#   select(Child, PM, Max_X = X, Max_Lag_Value = Lags)

merged_df <- inner_join(ds_final, ds.lag, by = c("Child", "PM"))
filtered_df <- merged_df %>% filter(Lags == Max_Lag_Value)
filtered_df <- select(filtered_df, -Lags)

Desc(filtered_df$Max_Lag_Value)
Desc(ds.lag$Max_Lag_Value)

m1 <- lmer(X ~  direction  + (1+cond|Child), 
           data=filtered_df)
Anova(m1)
summary(m1)
emmeans(m1, pairwise ~ direction, adjust="sidak")

ggplot(data = filtered_df, aes(x = direction, y = X,  col=direction)) +
  geom_boxplot() + 
  labs(title = "", x = "Lags", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "none")



m1 <- lmer(X ~  cond * age * Category  + (1+cond|Child), 
           data=filtered_df[filtered_df$direction == "music2move",])
Anova(m1,2)
summary(m1)
emmeans(m1, pairwise ~ Category|age, adjust="sidak")
emmeans(m1, pairwise ~ cond, adjust="sidak")




ggplot(data = filtered_df, aes(x = cond, y = X, group=interaction(cond,age), col=cond)) +
  geom_boxplot() + 
  labs(title = "", x = "Lags", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "none")+
  facet_grid(Category~age)

ggplot(data = filtered_df, aes(x = Category, y = X, group=interaction(cond,Category), col=cond)) +
  geom_boxplot() + 
  labs(title = "", x = "Category", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "none")+
  facet_grid(~age)





ggplot() +
  geom_boxplot(data = ds_final, aes(x = cond, y = X, col=cond)) +
  labs(title = "", x = "Condition", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "right")+
  facet_grid(Category~age)
  # ylim(0.5,2.5)


plotdf <- aggregate(X ~ cond * age * Lags * Category,  data=ds_final, FUN="mean", na.rm=T)
ggplot(data = plotdf, aes(x = Lags, y = X, group=interaction(cond,age), col=interaction(cond))) +
  geom_point()+
  geom_line() +
  geom_vline(xintercept=0,linetype="dotted")+
  geom_vline(xintercept=0.444,linetype="dotted")+
  geom_vline(xintercept=0.888,linetype="dotted")+
  geom_vline(xintercept=1.776,linetype="dotted")+
  geom_vline(xintercept=3.552,linetype="dotted")+
  geom_vline(xintercept=-0.444,linetype="dotted")+
  geom_vline(xintercept=-0.888,linetype="dotted")+
  geom_vline(xintercept=-1.776,linetype="dotted")+
  geom_vline(xintercept=-3.552,linetype="dotted")+
  labs(title = "", x = "Lags in seconds", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "right")+
  facet_grid(age~Category)+
  xlim(-1,4)




# combine them
granger_mus_F_range <- read_csv("granger_mus_F_range.csv", 
                                col_types = cols(...1 = col_skip()))


xdata <- rbind(granger_mus_F_range,ds_final)

xdata$Lags <- as.numeric(xdata$Lags)
xdata <- xdata %>%
  arrange(Child, PM, cond, Lags) %>%  # Assuming there's a Time column for ordering
  group_by(Child, PM, cond) %>%
  mutate(Diff_X = X - lag(X)) %>%
  ungroup()

xdata.plot <- aggregate(Diff_X ~ cond * age * Lags * Category ,  data=xdata, FUN="mean")
ggplot(data = xdata.plot, aes(x = Lags, y = Diff_X, group=interaction(age), col=age)) +
  geom_point()+
  geom_line() +
  labs(title = "", x = "Lags", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "right")+
  facet_grid(Category~cond)


# create super music condition
xdata.music <- xdata[xdata$cond == "Music",]
xdata.lp <- xdata[xdata$cond == "LowPitch",]
xdata.hp <- xdata[xdata$cond == "HighPitch",]
xdata.control <- xdata[xdata$cond == "Control",]

xdata.super <- merge(xdata.music, xdata.lp[,c("Child","PM","Lags","X")], by=c("Child","PM","Lags"),all.x=T)
xdata.super <- merge(xdata.super, xdata.hp[,c("Child","PM","Lags","X")], by=c("Child","PM","Lags"),all.x=T)
xdata.super$X_new <- (xdata.super$X.x+xdata.super$X.y+xdata.super$X)/3 
xdata.super$X <- xdata.super$X_new
xdata.super <- xdata.super[,c("Child","age","PM","Category","Lags","cond","X")]

xdata.control <- xdata.control[,c("Child","age","PM","Category","Lags","cond","X")]


xdata.super <- merge(xdata.music, xdata.control[,c("Child","PM","Lags","X")], by=c("Child","PM","Lags"),all.x=T)
xdata.super$X <- xdata.super$X.x-xdata.super$X.y
xdata.super.plot <- aggregate(X ~ cond * age * Lags * Category ,  data=xdata.super, FUN="mean")
ggplot(data = xdata.super.plot, aes(x = Lags, y = X, group=interaction(cond), col=cond)) +
  geom_point()+
  geom_line() +
  labs(title = "", x = "Condition", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "right")+
  facet_grid(Category~age)
