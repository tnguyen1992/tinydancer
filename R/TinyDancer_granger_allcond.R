# Auto-correlation script

library(readr)
library(dplyr)
library(pracma)
library(lmtest)
library(signal)

env_hopp_music <- read_csv("C:/Users/User/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_hopp_music.csv", 
                           col_names = FALSE)
env_hopp_control <- read_csv("C:/Users/User/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_hopp_control.csv", 
                             col_names = FALSE)
env_lola_music <- read_csv("C:/Users/User/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_lola_music.csv", 
                           col_names = FALSE)
env_lola_control <- read_csv("C:/Users/User/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/AcousticStim/env_lola_control.csv", 
                             col_names = FALSE)


# plot(env_hopp_music$X1, type='l', col="red")
# lines(env_hopp_lp$X1, col="orange")
# lines(env_hopp_hp$X1, col="green")
# lines(stim$X1*0.01, col="blue")

# sd(env_hopp_music$X1)
# sd(env_hopp_lp$X1)
# sd(env_hopp_hp$X1)

# Replace "your_folder_path" with the actual path to your folder
ages<- c('3m', "6m", "12m")
conditions <- c("HBaseline", "HControl","LBaseline", "LControl")
# lags <- c(1,seq(11,165,by=11),seq(6,165,by=11))
lags <- seq(1,90)
pm_frame_x_all <- list()  # Assuming the third dimension is for x and y values

for (age in ages) {
  if (age == "6m") {
    NB <- 25
  } else {
    NB <- 24
  }
  # Set condition
  # Set the range of baby numbers
  baby_numbers <- sprintf("%02d", 1:NB)  # Formatting to ensure two-digit representation (e.g., 01, 02, ..., 24)
  folder_path <- paste0("C:/Users/tnguyen/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/Deeplabcut/TinyDancer_PM/tinyDancers_commonAges/tinyDancers_commonAges/csv_for_trinh/age", age)
  
  # Create an empty data frame to store results
  pm_frame_x_list <- list()  # Assuming the third dimension is for x and y values
  
  for (condition in conditions) {
    granger_F_values <- array(NA, dim = c(10, NB,length(lags)))
    
    if (condition=="HBaseline") {
      stim=env_hopp_music
    } 
    if (condition=="HControl"){
      stim=env_hopp_control
    } 
    if (condition=="LControl"){
      stim=env_lola_control
    } 
    if (condition=="LBaseline") {
      stim=env_lola_music
    }
    # stim_ts <- ts(scale(stim[1:524,],center=T, scale=T) )       
    stim_ts <- ts(stim[1:524,])         
    
    # acf(stim_ts,plot = T, demean=T,lag.max=50)
    
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
          transposed_vel <- diff(transposed_df) # get to velocity
          
          velocity_PM <- transposed_vel[, PM]
          velocity_PM_sample <- sample(velocity_PM)
          
          for (l in seq_along(lags)) {
            # Perform Granger causality test
            tsDat <- ts.union(stim_ts, velocity_PM) 
            tsVAR <- vars::VAR(tsDat, p = l)
            gctest <- vars::causality(tsVAR, cause = "stim_ts")$Granger
            # gctest <- vars::causality(tsVAR, cause = "velocity_PM")$Granger
            
            # gctest <- suppressWarnings(grangertest(stim, velocity_PM, order = lags[l], na.action = na.omit))
            # gctest <- suppressWarnings(grangertest(velocity_PM, stim, order = lags[l], na.action = na.omit))
            
            collectF[i,l] <- gctest$statistic
          }
        }
        granger_F_values[PM,as.numeric(baby_number),] <- colMeans(collectF,na.rm=T)
        
      }
    }
    pm_frame_x_list[[condition]] <- granger_F_values
    
  }
  pm_frame_x_all[[age]] <- pm_frame_x_list
  
}






pm_frame_x_hcon_3m <-pm_frame_x_all[["3m"]][["HControl"]]
pm_frame_x_lcon_3m <-pm_frame_x_all[["3m"]][["LControl"]]

# Create a new 10 x 24 array by calculating the mean of corresponding values
pm_frame_x_con_3m <- array(NA, dim=c(10,24,length(lags)))

for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_con_3m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_3m[i, j,], pm_frame_x_lcon_3m[i, j,]))
  }
}


pm_frame_x_hmus_3m <-pm_frame_x_all[["3m"]][["HBaseline"]]
pm_frame_x_lmus_3m <-pm_frame_x_all[["3m"]][["LBaseline"]]

pm_frame_x_mus_3m <- array(NA, dim=c(10,24,length(lags)))

for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_mus_3m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_3m[i, j,], pm_frame_x_lmus_3m[i, j,]))
  }
}

pm_frame_x_hcon_6m <-pm_frame_x_all[["6m"]][["HControl"]]
pm_frame_x_lcon_6m <-pm_frame_x_all[["6m"]][["LControl"]]
pm_frame_x_con_6m <- array(NA, dim=c(10,25,length(lags)))
for (i in 1:10) {
  for (j in 1:25) {
    pm_frame_x_con_6m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_6m[i, j,], pm_frame_x_lcon_6m[i, j,]))
  }
}

pm_frame_x_hmus_6m <-pm_frame_x_all[["6m"]][["HBaseline"]]
pm_frame_x_lmus_6m <-pm_frame_x_all[["6m"]][["LBaseline"]]
pm_frame_x_mus_6m <- array(NA, dim=c(10,25,length(lags)))
for (i in 1:10) {
  for (j in 1:25) {
    pm_frame_x_mus_6m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_6m[i, j,], pm_frame_x_lmus_6m[i, j,]))
  }
}

pm_frame_x_hcon_12m <-pm_frame_x_all[["12m"]][["HControl"]]
pm_frame_x_lcon_12m <-pm_frame_x_all[["12m"]][["LControl"]]
pm_frame_x_con_12m <- array(NA, dim=c(10,24,length(lags)))
for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_con_12m[i, j,] <- colMeans(rbind(pm_frame_x_hcon_12m[i, j,], pm_frame_x_lcon_12m[i, j,]))
  }
}

pm_frame_x_hmus_12m <-pm_frame_x_all[["12m"]][["HBaseline"]]
pm_frame_x_lmus_12m <-pm_frame_x_all[["12m"]][["LBaseline"]]
pm_frame_x_mus_12m <- array(NA, dim=c(10,24,length(lags)))
for (i in 1:10) {
  for (j in 1:24) {
    pm_frame_x_mus_12m[i, j,] <- colMeans(rbind(pm_frame_x_hmus_12m[i, j,], pm_frame_x_lmus_12m[i, j,]))
  }
}





##################################
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

con_3m_df$cond <- 'Control'
mus_3m_df$cond <- 'Music'

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

con_6m_df$cond <- 'Control'
mus_6m_df$cond <- 'Music'

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

con_12m_df$cond <- 'Control'
mus_12m_df$cond <- 'Music'

cond_12m_df <- rbind(mus_12m_df,con_12m_df)

cond_3m_df$age <- '3m'
cond_6m_df$age <- '6m'
cond_12m_df$age <- '12m'

cond_6m_df$Child <- cond_6m_df$Child+100
cond_12m_df$Child <- cond_12m_df$Child+200

cond_all_df <- rbind(cond_3m_df,cond_6m_df,cond_12m_df)


# write.csv2(cond_all_df,"granger_mus_F_range_minuslags_upsample.csv")
# write.csv2(cond_all_df,"granger_mus_F_range_upsample_110324.csv")

# write.csv2(cond_all_df,"granger_mus_F_minusrange_upsample_110324.csv")
write.csv2(cond_all_df,"granger_mus_F_samplerange_010724.csv")

##########################################
library(ggplot2)
library(gridExtra )
library(glmmTMB )
library(car )
library(effects)
library(emmeans)
library(DescTools)
library(readr)

setwd("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/Projects/WP4/Deeplabcut/TinyDancer_PM")
granger_mus_F_range <- read_delim("granger_mus_F_range_upsample_110324.csv", 
                                  delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                  locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                  trim_ws = TRUE)
granger_mus_F_range$Lags <- as.numeric(granger_mus_F_range$Lags)
granger_mus_F_range$direction <- "music2move"
granger_mus_F_range_minuslags <- read_delim("granger_mus_F_minusrange_upsample_110324.csv", 
                                            delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                            locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                            trim_ws = TRUE)

granger_mus_F_range_minuslags$Lags <- as.numeric(granger_mus_F_range_minuslags$Lags)
granger_mus_F_range_minuslags$direction <- "move2music"



ds_final <- rbind(granger_mus_F_range_minuslags,granger_mus_F_range)
ds_final$PM <- as.factor(ds_final$PM)
# ds_final$X.w <- Winsorize(ds_final$X,na.rm=T)
ds_final$age <- factor(ds_final$age, levels=c('3m', '6m', '12m'))
ds_final$Lags <- as.numeric(ds_final$Lags)

##### ttest
# Assuming your data frame is named ds_final

# Define the age groups, PMs, and conditions
age_groups <- c("3m", "6m", "12m")
PMs <- 1:10
Lags <- unique(ds_final$Lags)
Categories<- c("Lower body","Hybrid","Upper body")

#### LMM
library(lme4)


# Assuming ds_final has a column named PM representing PM numbers
ds_final$Category <- rep(NA, nrow(ds_final))

# Identify Upper body, Lower body, and Hybrid PMs
upper_body_pms <- c(1, 2, 3, 5, 6)
lower_body_pms <- c(4, 7, 9, 10)
hybrid_pms <- c( 8)

# Assign categories based on PM numbers
ds_final$Category[ds_final$PM %in% upper_body_pms] <- "Upper body"
ds_final$Category[ds_final$PM %in% lower_body_pms] <- "Lower body"
ds_final$Category[ds_final$PM %in% hybrid_pms] <- "Hybrid"

# Convert Category to a factor
ds_final$Category <- factor(ds_final$Category, levels = c("Upper body", "Lower body", "Hybrid"))

ds_final$Lags.f <- factor(ds_final$Lags)

ds_final_check <- ds_final[!ds_final$Child==20,]
ds_final_check <- ds_final_check[!ds_final_check$Child==103,]
ds_final_check <- ds_final_check[!ds_final_check$Child==211,]
ds_final_check <- ds_final_check[!ds_final_check$Child==214,]

m1 <- lmer(X ~  age  * PM * direction + +(1|PM)+(1|Child),
           data=ds_final_check[ds_final_check$cond=="Music",])
Anova(m1,3)
summary(m1)
emmeans(m1, pairwise ~ direction, adjust="sidak")


plotdf <- aggregate(X ~ cond * age * Lags.f   * direction,  data=ds_final, FUN="mean", na.rm=T)
plotdf.sd <- aggregate(X ~ cond * age * Lags.f  * direction,  data=ds_final, FUN="sd", na.rm=T)
colnames(plotdf.sd) <- list("cond", "age", "Lags.f","direction","SD")
plotdf.f <- merge(plotdf,plotdf.sd, by=c("age", "Lags.f","direction", "cond"),all.X=T)
plotdf.f$Lags <- as.numeric(plotdf.f$Lags.f)
ggplot(data = plotdf.f[plotdf.f$Lags<20 &plotdf.f$cond=="Music",] , aes(x = Lags, y = X, group=interaction(direction,age), col=direction)) +
  # geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
  # geom_line(aes(y = sig_line),
  #           color = "black", size = 1) +
  labs(title = "", x = "Lags in frames (sr: 25 Hz)", y = "Granger Causality F-value")+
  theme_bw()+
  theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
  facet_grid(~age)


library(MKinfer)

library(dplyr)
ds_final.superPM <- aggregate(X ~ cond * age * Lags.f * Child * direction,  data=ds_final, FUN="mean", na.rm=T)
ds_final.music2move<-ds_final.superPM[ds_final.superPM$direction=="music2move",]
ds_final$cond <- factor(ds_final$cond)

ds_final_check <- ds_final.music2move[!ds_final.music2move$Child==20,]
ds_final_check <- ds_final_check[!ds_final_check$Child==103,]
ds_final_check <- ds_final_check[!ds_final_check$Child==211,]
ds_final.music2move <- ds_final_check[!ds_final_check$Child==214,]


##### Plots per PM

result <- ds_final.music2move %>%
  group_by(Lags.f, age) %>%
  summarize(t_test_p_value = boot.t.test(X ~ cond)$boot.p.value)

# Display the result
print(result)

result <- result %>%
  mutate(sig_line = ifelse(t_test_p_value < 0.005, 0.5, NA))
plotdf <- aggregate(X ~ cond * age * Lags.f  ,  data=ds_final, FUN="mean", na.rm=T)

plotdf.sd <- aggregate(X ~ cond * age * Lags.f  ,  data=ds_final, FUN="sd", na.rm=T)
colnames(plotdf.sd) <- list("cond", "age", "Lags.f","SD")

plotdf.f <- merge(plotdf,plotdf.sd, by=c("age", "Lags.f", "cond"),all.X=T)
plotdf.f <- merge(plotdf.f,result, by=c("age", "Lags.f"))

plotdf.f$Lags <- as.numeric(plotdf.f$Lags.f)
ggplot(data = plotdf.f[plotdf.f$Lags<15 ,], aes(x = Lags, y = X, group=interaction(cond,age), col=cond)) +
  # geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
  geom_line(aes(y = sig_line),
            color = "black", size = 1) +
  labs(title = "", x = "Lags in frames (sr: 25 Hz)", y = "Granger Causality F-value")+
  theme_bw()+
  theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
  facet_grid(1~age)

##### Plots for super PM

result <- ds_final.music2move %>%
  group_by(Lags.f, age) %>%
  summarize(t_test_t_value = boot.t.test(X ~ cond)$statistic,
            t_test_df_value = boot.t.test(X ~ cond)$parameter,
            t_test_p_value = boot.t.test(X ~ cond)$boot.p.value)

result <- result %>%
  mutate(sig_line = ifelse(t_test_p_value < 0.005, 0.5, NA))
plotdf <- aggregate(X ~ cond * age * Lags.f  ,  data=ds_final, FUN="mean", na.rm=T)

plotdf.sd <- aggregate(X ~ cond * age * Lags.f  ,  data=ds_final, FUN="sd", na.rm=T)
colnames(plotdf.sd) <- list("cond", "age", "Lags.f","SD")

plotdf.f <- merge(plotdf,plotdf.sd, by=c("age", "Lags.f", "cond"),all.X=T)
plotdf.f <- merge(plotdf.f,result, by=c("age", "Lags.f"))



plotdf.f$Lags <- as.numeric(plotdf.f$Lags.f)
ggplot(data = plotdf.f[plotdf.f$Lags<15 ,], aes(x = Lags, y = X, group=interaction(cond,age), col=cond)) +
  # geom_point()+
  geom_line() +
  geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
  # geom_line(aes(y = sig_line),
  #           color = "black", size = 1) +
  labs(title = "", x = "Lags in frames (sr: 25 Hz)", y = "Granger Causality F-value")+
  theme_bw()+
  theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
  facet_grid(~age)
##################################

ds_final.music2move <- ds_final[ds_final$direction=="music2move" ,]
ds_final.move2music <- ds_final[ds_final$direction=="move2music" ,]

ds.directioncontrast <- merge(ds_final.music2move,ds_final.move2music, by=c("Child","PM","Lags","cond"))
ds.directioncontrast <- ds.directioncontrast %>% select(Child, PM,Lags,age.x,cond, Category.x,X.x, X.y)
ds.directioncontrast$X <- ds.directioncontrast$X.x - ds.directioncontrast$X.y
ds.directioncontrast$Category <- ds.directioncontrast$Category.x
result <- ds.directioncontrast %>%
  group_by(Lags, Category,age.x) %>%
  summarize(t_test_p_value = t.test(X ~ cond)$p.value, alternative="less")
result <- result %>%
  mutate(sig_line = ifelse(t_test_p_value < 0.05, -1, NA))
# Display the result
print(result)


plotdf.directioncontrast <- aggregate(X ~ cond *  Lags * Category * age.x,  data=ds.directioncontrast, FUN="mean", na.rm=T)
plotdf.directioncontrast <- merge(plotdf.directioncontrast,result, by=c("age.x", "Category", "Lags"))


ggplot(data = plotdf.directioncontrast, aes(x = Lags, y = X, group=interaction(cond), col=interaction(cond))) +
  # geom_point()+
  geom_line() +
  # geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
  # geom_vline(xintercept=0,linetype="dotted")+
  # geom_vline(xintercept=1,linetype="dotted")+
  # geom_vline(xintercept=2,linetype="dotted")+
  # geom_vline(xintercept="0.5",linetype="dotted")+
  # geom_vline(xintercept="1.5",linetype="dotted")+
  geom_line(aes(y = sig_line),
            color = "red", size = 1) +
  labs(title = "", x = "Lags in seconds", y = "Granger Causality F-value")+
  theme_minimal()+
  theme(legend.position = "top")+
  facet_grid(age.x~Category)
#######################################################
granger_mus_F_range <- read_delim("granger_mus_F_range_upsample_110324.csv", 
                                  delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                  locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                  trim_ws = TRUE)
granger_mus_F_range$Lags <- as.numeric(granger_mus_F_range$Lags)
granger_mus_F_range$direction <- "music2move"
granger_mus_F_range_minuslags <- read_delim("granger_mus_F_minusrange_upsample_110324.csv", 
                                            delim = ";", escape_double = FALSE, col_types = cols(...1 = col_skip()), 
                                            locale = locale(decimal_mark = ",", grouping_mark = "."), 
                                            trim_ws = TRUE)

granger_mus_F_range_minuslags$Lags <- as.numeric(granger_mus_F_range_minuslags$Lags)
granger_mus_F_range_minuslags$direction <- "move2music"


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


ds_final <- rbind(granger_mus_F_range_minuslags,granger_mus_F_range,granger_freq_F_range_minuslags,granger_freq_F_range)
# ds_final <- granger_mus_F_range
ds_final <- rbind(granger_freq_F_range_minuslags,granger_freq_F_range)
# ds_final <- rbind(granger_mus_F_range_minuslags,granger_mus_F_range)


ds_final$PM <- as.factor(ds_final$PM)
# ds_final$X.w <- Winsorize(ds_final$X,na.rm=T)
ds_final$age <- factor(ds_final$age, levels=c('3m', '6m', '12m'))
ds_final$Lags <- as.numeric(ds_final$Lags)
ds_final$cond <- as.factor(ds_final$cond)

##### ttest
# Assuming your data frame is named ds_final

# Define the age groups, PMs, and conditions
age_groups <- c("3m", "6m", "12m")
PMs <- 1:10
Lags <- unique(ds_final$Lags)
Categories<- c("Lower body","Hybrid","Upper body")

#### LMM
library(lme4)
library(car)
library(emmeans)

# Assuming ds_final has a column named PM representing PM numbers
ds_final$Category <- rep(NA, nrow(ds_final))

# Identify Upper body, Lower body, and Hybrid PMs
upper_body_pms <- c(1, 2, 3, 5, 6)
lower_body_pms <- c(4, 7, 9, 10)
hybrid_pms <- c(8)

# Assign categories based on PM numbers
ds_final$Category[ds_final$PM %in% upper_body_pms] <- "Upper body"
ds_final$Category[ds_final$PM %in% lower_body_pms] <- "Lower body"
ds_final$Category[ds_final$PM %in% hybrid_pms] <- "Hybrid"

# Convert Category to a factor
ds_final$Category <- factor(ds_final$Category, levels = c("Upper body", "Lower body", "Hybrid"))

ds_final$Lags.f <- factor(ds_final$Lags)
xdata <- ds_final[ds_final$Lags >= 0.20 & ds_final$Lags <= 0.20 ,]
xdata <- aggregate(X ~ age * Child * PM * Category *direction  * cond,  data=xdata, FUN="mean", na.rm=T)
# 
# m0 <- lmer(X ~  direction + (1|Child),
#            data=ds_final)
# Anova(m0,2)
# summary(m0)
# emmeans(m0, pairwise ~ direction, adjust="sidak")
# 
# data <- data.frame(
#   direction = c("Music to Movement", "Movement to Music"),
#   mean = c(1.07, 1.04 ),
#   se = c(0.00470, 0.00643)
# )
# data <- data.frame(
#   direction = c("Music to Movement", "Movement to Music"),
#   mean = c(1.06, 1.05 ),
#   se = c(0.00546, 0.00546)
# )
# ggplot(data, aes(x = direction, y = mean, col = direction)) +
#   geom_point(stat = "identity", position = "dodge", width = 0.7) +
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.7), width = 0.25) +
#   labs(title = "Information flow direction",
#        x = "Direction",
#        y = "Mean and SE") +
#   theme_minimal()+
#   theme(legend.position = "none")+
#   ylim(1,1.1)
# 
# 
# library(datawizard)
# xdata$X.w <- winsorize(xdata$X, method="zscore", threshold=3)
# 
# xdata.f<-xdata[xdata$cond=="HighPitch" | xdata$cond=="LowPitch" ,]
# xdata.m<-xdata[xdata$cond=="Music" | xdata$cond=="Control" ,]
# 
# xdata.m$X.w2 <- winsorize(xdata.m$X, method="zscore", threshold=3)
# xdata.f$X.w2 <- winsorize(xdata.f$X, method="zscore", threshold=3)
# 
# hist(xdata.m$X)
# hist(xdata.m$X.w)
# hist(xdata.m$X.w2)
# 
# m1 <- lmer(X.w2 ~  age  * PM * cond + (1+cond|Child),
#            data=xdata.m[xdata.m$direction=="music2move",])
# Anova(m1,3)
# summary(m1)
# emmeans(m1, pairwise ~ cond|PM|age, adjust="fdr")
# emmeans(m1, pairwise ~ cond|age, adjust="fdr")
# 
# 
# hist(xdata.f$X)
# range(xdata.f$X)
# mean(xdata.f$X)+3*sd(xdata.f$X)
# 
# m2 <- lmer(X.w2 ~  age  * PM * cond + (1+cond|Child),
#            data=xdata.f[xdata.f$direction=="music2move",])
# Anova(m2,2)
# summary(m2)
# emmeans(m2, pairwise ~ cond, adjust="fdr")
# emmeans(m2, pairwise ~ cond|PM|age, adjust="fdr")
# emmeans(m2, pairwise ~ cond|age, adjust="fdr")


# 
# 
# # plot
# desired_order <- c("3m", "6m", "12m")
# xdata.m$age <- ordered(xdata.m$age, levels = desired_order)
# xdata.f$age <- ordered(xdata.f$age, levels = desired_order)
# 
# xdata.m_summary <- ddply(xdata.m,.(PM,cond,age, direction),summarize,avg_value=mean(X,na.rm=TRUE),se_value=sd(X,na.rm=TRUE)/length(X))
# xdata.f_summary <- ddply(xdata.f,.(PM,cond,age, direction),summarize,avg_value=mean(X,na.rm=TRUE),se_value=sd(X,na.rm=TRUE)/length(X))
# 
# xdata.m_summary_superPM <- ddply(xdata.m,.(cond,age,direction),summarize,avg_value=mean(X,na.rm=TRUE),se_value=sd(X,na.rm=TRUE)/length(X))
# xdata.f_summary_superPM <- ddply(xdata.f,.(cond,age,direction),summarize,avg_value=mean(X,na.rm=TRUE),se_value=sd(X,na.rm=TRUE)/length(X))
# 
# 
# ggplot(xdata.m_summary_superPM[xdata.m_summary_superPM$direction=="music2move",], aes(x = age, y = avg_value, fill = cond, group = cond)) +
#   geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
#   geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
#                 position="dodge",width=0.5, linewidth=0.5)+
#   # geom_boxplot() +
#   # geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.7), width = 0.25) +
#   labs(title = "GC Music vs Control",
#        x = "Conditions",
#        y = "Granger F-Values") +
#   theme_minimal()+
#   theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
#   ylim(0,1.5)
# 
# ggplot(xdata.m_summary[xdata.m_summary$direction=="music2move",], aes(x = age, y = avg_value, fill = cond, group = cond)) +
#   geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
#   geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
#                 position="dodge",width=0.5, linewidth=0.5)+
#   # geom_boxplot() +
#   # geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.7), width = 0.25) +
#   labs(title = "GC Music vs Control",
#        x = "Conditions",
#        y = "Granger F-Values") +
#   theme_minimal()+
#   theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
#   facet_grid(PM~"")
# 
# 
# ggplot(xdata.f_summary_superPM[xdata.m_summary_superPM$direction=="music2move",], aes(x = age, y = avg_value, fill = cond, group = cond)) +
#   geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
#   geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
#                 position="dodge",width=0.5, linewidth=0.5)+
#   # geom_boxplot() +
#   # geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.7), width = 0.25) +
#   labs(title = "GC HP vs LP",
#        x = "Conditions",
#        y = "Granger F-Values") +
#   theme_minimal()+
#   theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
#   ylim(0,1.5)
# 
# ggplot(xdata.f_summary[xdata.m_summary$direction=="music2move",], aes(x = age, y = avg_value, fill = cond, group = cond)) +
#   geom_col(position="dodge", width=0.6)+ # , position=position_jitter(5)
#   geom_errorbar(aes(ymax=avg_value+se_value,ymin=avg_value-se_value),
#                 position="dodge",width=0.5, linewidth=0.5)+
#   # geom_boxplot() +
#   # geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.7), width = 0.25) +
#   labs(title = "GC HP vs LP",
#        x = "Conditions",
#        y = "Granger F-Values") +
#   theme_minimal()+
#   theme(legend.position = "right",panel.grid = element_blank(),strip.text=element_blank())+
#   facet_grid(PM~"")
# 
# 
# 
# library(dplyr)
# result <- ds_final %>%
#   group_by(Lags.f, Category, age,direction) %>%
#   summarize(t_test_p_value = t.test(X ~ cond, alternative="two.sided")$p.value)
# 
# # Display the result
# print(result)
# 
# ##### Plots
# 
# 
# result <- result %>%
#   mutate(sig_line = ifelse(t_test_p_value < 0.005, 0, NA))
# 
# 
# # plotdf.f.wc <- ds_final[!ds_final$cond=="Control",]
# 
# plotdf <- aggregate(X ~ age * Lags.f * Category * direction * cond,  data=ds_final, FUN="mean", na.rm=T)
# 
# plotdf.sd <- aggregate(X ~ age * Lags.f * Category * direction * cond,  data=ds_final, FUN="SD", na.rm=T)
# colnames(plotdf.sd) <- list("age", "Lags.f","Category","direction","cond","SD")
# 
# plotdf.f <- merge(plotdf,plotdf.sd, by=c("age", "Category", "Lags.f","direction","cond"),all.X=T)
# plotdf.f <- merge(plotdf.f,result, by=c("age", "Category", "Lags.f","direction"))
# 
# plotdf.f$Lags <- as.numeric(plotdf.f$Lags)
# ggplot(data = plotdf.f[plotdf.f$direction=="music2move" & plotdf.f$Lags < 45,], aes(x = Lags, y = X, group=interaction(cond,Category), col=interaction(cond,Category))) +
#   # geom_point()+
#   geom_line() +
#   geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
#   # geom_vline(xintercept=4,linetype="dotted")+
#   # geom_vline(xintercept=5,linetype="dotted")+
#   # geom_vline(xintercept=12,linetype="dotted")+
#   # geom_vline(xintercept=19,linetype="dotted")+
#   # geom_vline(xintercept="3.52",linetype="dotted")+
#   # geom_line(aes(y = sig_line),
#   # color = "black", size = 1) +
#   labs(title = "", x = "Lags in seconds", y = "Granger Causality F-value")+
#   theme_minimal()+
#   theme(legend.position = "top")+
#   facet_grid(~age)
# 
# ggplot(data = plotdf.f[plotdf.f$direction=="music2move",], aes(x = Lags.f, y = X, group=interaction(cond), col=interaction(cond))) +
#   # geom_point()+
#   geom_line() +
#   geom_ribbon(aes(ymin = X - (SD/sqrt(24)), ymax = X + (SD/sqrt(24))), alpha = 0.2) +  
#   geom_vline(xintercept=0,linetype="dotted")+
#   geom_vline(xintercept="0.44",linetype="dotted")+
#   geom_vline(xintercept="0.88",linetype="dotted")+
#   geom_vline(xintercept="1.76",linetype="dotted")+
#   geom_vline(xintercept="3.52",linetype="dotted")+
#   # geom_line(aes(y = sig_line),
#             # color = "red", size = 1) +
#   labs(title = "", x = "Lags in seconds", y = "Granger Causality F-value")+
#   theme_minimal()+
#   theme(legend.position = "top")+
#   facet_grid(age~PM)


### combine with EEG and test with P50 or P200

xdata$Child <- factor(xdata$Child,levels = c("1","2","3","4","5","6","7",
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




ERP_all_auc <- read_excel("~/OneDrive - Fondazione Istituto Italiano Tecnologia/IIT_Postdoc/WP4/MATLAB/MUSICOM_R/ERP_all_auc.xls")
library(dplyr)

# Rename condition variable
ERP_all_auc<- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "baseline", "Music", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "control", "Control", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "lowbass", "LowPitch", condition))
ERP_all_auc <- ERP_all_auc %>%
  mutate(condition = ifelse(condition == "highvoice", "HighPitch", condition))

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

xdata$cond <- as.factor(xdata$cond)
xdata$age <- as.factor(xdata$age)

em.data <- merge(xdata,df, by=c("Child","cond", "age"), all.y=T)
em.data1d <- em.data[em.data$direction=="music2move",]

em.data.mus <- em.data1d[em.data1d$cond=="Music" | em.data1d$cond=="Control",]
em.data.freq <- em.data1d[em.data1d$cond=="HighPitch"|em.data1d$cond=="LowPitch",]

em.data.mus$amplitude.w.z <- scale(Winsorize(em.data.mus$amplitude, na.rm=T))


m1 <- lmer(X ~ amplitude  * cond * age + (1+amplitude|Child), #+ I(amplitude.w^2)
           data=em.data.freq, na.action=na.omit)
summary(m1)
Anova(m1,2)

plot(effect('amplitude:cond',m1))
plot(effect('amplitude:cond:age',m1))

emtrends(m1, pairwise~cond|age,var="amplitude" )
emtrends(m1, pairwise~cond,var="amplitude" )
emtrends(m1, pairwise~age,var="amplitude" )

## load assr for assr granger analyses

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
power_all$power.w.z <-DescTools::Winsorize(power_all$power.z)


power_mus <- power_all
power_mus$age <- factor(power_mus$age,levels=c("3","6","12","ad"), labels=c("3m","6m","12m","ad"))
power_mus$condition <- factor(power_mus$condition,levels=c("music","control", "highvoice", "lowbass"), labels=c("Music","Control", "HighPitch", "LowPitch"))

power_mus$cond <- as.factor(power_mus$condition)
power_mus$age <- as.factor(power_mus$age)
power_mus$Child <- as.factor(power_mus$ID)

power_mus.data <- merge(xdata,power_mus, by=c("Child","cond", "age"), all=T)

power_mus.data1d <- power_mus.data[power_mus.data$direction=="music2move",]

am.data.mus <- power_mus.data1d[power_mus.data1d$cond=="Music" | power_mus.data1d$cond=="Control",]
am.data.freq <- power_mus.data1d[power_mus.data1d$cond=="HighPitch"|power_mus.data1d$cond=="LowPitch",]



m1 <- lmer(X ~  cond * age*power.w.z + (1|Child), #+ I(amplitude.w^2)
           data=am.data.freq, na.action=na.omit)
summary(m1)
Anova(m1,3)

plot(effect('cond:power.w.z',m1))
plot(effect('amplitude:cond:age',m1))

emtrends(m1, pairwise~cond|age,var="amplitude" )
emtrends(m1, pairwise~cond,var="amplitude" )
emtrends(m1, pairwise~age,var="amplitude" )