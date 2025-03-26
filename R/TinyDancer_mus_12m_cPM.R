library(readr)
library(circular)
library(plotrix)
library(ggplot2)
library(gsignal)
library(dplyr)

# Replace "your_folder_path" with the actual path to your folder
age="12m"
folder_path <- paste0("//clnsd009/Users/fbigand/OneDrive - Fondazione Istituto Italiano Tecnologia/WORK/COLLABS/Trinh/Tiny-dancers/tinyDancers_commonAges/csv_for_trinh/age",age)

# Create an empty data frame to store results
baby_results_df <- data.frame()
baby_results <- list()

# Set the range of baby numbers
baby_numbers <- sprintf("%02d", 1:24)  # Formatting to ensure two-digit representation (e.g., 01, 02, ..., 24)
PM_numbers <- sprintf("%02d", 1:10)  # Formatting to ensure two-digit representation (e.g., 01, 02, ..., 24)

# Set condition
condition="WT_mus"

# Loop through each baby
for (baby_number in baby_numbers) {
  # Baseline Condition
  # List all CSV files in the folder that contain the word "baseline"
  pattern<-paste0(condition,".*\\.csv$")
  file_names <- list.files(paste0(folder_path,"/babe",baby_number,"/QoM-WT_pm",PM_numbers), pattern = pattern,recursive=T, full.names = TRUE)
  
  # Read, transpose, and stack files dynamically
  list_of_dfs <- lapply(file_names, function(file) {
    df <- read.csv(file)
    # colnames(df)<-list("period", "PM_music", "PM_control")
    # transposed_df <- t(df)
    # transposed_df <- transposed_df[-1, ]  # Exclude the first row
  })
  freq_data <- do.call(rbind, list_of_dfs)
  colnames(freq_data) <- list("period", "control", "music")
  freq_data$baby <- baby_number
  freq_data$PM <- rep(1:10, each = 81)
  # Convert the list to a data frame
  baby_results_df <- bind_rows(baby_results_df,freq_data)
  
}
write.csv(baby_results_df, "mus_12m_baby.csv", row.names = FALSE)


#######################################################
library(reshape2)
library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(emmeans)
library(permutes)

baby_results_df <- read_csv("mus_12m_baby.csv")
# Specify id.vars: the variables to keep but not split apart on
freq_Music_12m_df <- melt(baby_results_df, id.vars=c("baby","period","PM"))

freq_Music_12m_df$PM <- factor(freq_Music_12m_df$PM)
freq_Music_12m_df$condition <- factor(freq_Music_12m_df$variable)
#freq_Music_12m_df$period <- factor(freq_Music_12m_df$period)



period_bins <- read_csv("period_bins.csv")
colnames(period_bins) <- list("period", "beatinfo")
freq_Music_12m_df<-merge(freq_Music_12m_df,period_bins, by="period")

freq_Music_12m_df$period <- as.numeric(freq_Music_12m_df$period)
freq_Music_12m_df$beatinfo <- as.numeric(freq_Music_12m_df$beatinfo)

freq_Music_12m_df <- freq_Music_12m_df %>%
  filter(beatinfo <= 10.00)
period_bins <- period_bins %>%
  filter(beatinfo <= 10.00)

# extract each PM
freq_Music_12m_PM1_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 1, ]
freq_Music_12m_PM2_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 2, ]
freq_Music_12m_PM3_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 3, ]
freq_Music_12m_PM4_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 4, ]
freq_Music_12m_PM5_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 5, ]
freq_Music_12m_PM6_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 6, ]
freq_Music_12m_PM7_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 7, ]
freq_Music_12m_PM8_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 8, ]
freq_Music_12m_PM9_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 9, ]
freq_Music_12m_PM10_df <- freq_Music_12m_df[freq_Music_12m_df$PM == 10, ]


# only analyse PMs till 10x beat

m1 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM1_df,
                       series.var=~beatinfo,
                       type='anova')
m1
summary(m1)
indices <- m1$p.cluster_mass < 0.005
m1$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]



m2 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM2_df,
                       series.var=~beatinfo,
                       type='anova')
m2
summary(m2)
indices <- m2$p.cluster_mass < 0.005
m2$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]


m3 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM3_df,
                       series.var=~beatinfo,
                       type='anova')
m3
summary(m3)
indices <- m3$p.cluster_mass < 0.005
m3$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]


m4 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM4_df,
                       series.var=~beatinfo,
                       type='anova')
m4
summary(m4)
indices <- m4$p.cluster_mass < 0.005
m4$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

m5 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM5_df,
                       series.var=~beatinfo,
                       type='anova')
m5
summary(m5)
indices <- m5$p.cluster_mass < 0.005
m5$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

m6 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM6_df,
                       series.var=~beatinfo,
                       type='anova')
m6
summary(m6)
indices <- m6$p.cluster_mass < 0.005
m6$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

m7 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM7_df,
                       series.var=~beatinfo,
                       type='anova')
m7
summary(m7)
indices <- m7$p.cluster_mass < 0.005
m7$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

m8 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM8_df,
                       series.var=~beatinfo,
                       type='anova')
m8
summary(m8)
indices <- m8$p.cluster_mass < 0.005
m8$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]


m9 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                       data=freq_Music_12m_PM9_df,
                       series.var=~beatinfo,
                       type='anova')
m9
summary(m9)
indices <- m9$p.cluster_mass < 0.005
m9$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

m10 <- clusterperm.lmer(value ~ variable + (variable|baby), 
                        data=freq_Music_12m_PM10_df,
                        series.var=~beatinfo,
                        type='anova')
m10
summary(m10)
indices <- m10$p.cluster_mass < 0.005
m10$beatinfo[indices]
every_second_value <- indices[seq(2, length(indices), by = 2)]
period_bins$beatinfo[every_second_value]

library(ggplot2)
library(wesanderson)
library(gridExtra)

# Assuming perm_results contains the relevant information
PM1_plot <- ggplot() +
  geom_line(data = freq_Music_12m_PM1_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM1", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')


PM2_plot <-ggplot() +
  geom_rect(aes(xmin = 0.5000000, xmax = 0.5693788, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_rect(aes(xmin = 0.6771398, xmax = 0.9575792, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM2_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM2", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM3_plot <-ggplot() +
  geom_rect(aes(xmin = 5.186722, xmax = 7.336757, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM3_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM3", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM4_plot <-ggplot() +
  geom_rect(aes(xmin = 6.724950, xmax = 9.107468, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM4_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM4", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM5_plot <-ggplot() +
  geom_rect(aes(xmin = 1.915342, xmax = 8.354219, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM5_df, aes(x = beatinfo, y = value, color=variable), method="gam") +
  labs(title = "Periodicity in PM5", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM6_plot <-ggplot() +
  geom_rect(aes(xmin = 3.668379, xmax = 5.186722, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM6_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM6", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM7_plot <-ggplot() +
  geom_rect(aes(xmin = 0.5000000, xmax = 0.5693788, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM7_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM7", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM8_plot <-ggplot() +
  geom_rect(aes(xmin = 2.953337, xmax = 3.829950, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM8_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM8", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')


PM9_plot <-ggplot() +
  geom_rect(aes(xmin = 0.5000000, xmax = 0.7383888, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM9_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM9", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

PM10_plot <-ggplot() +
  geom_rect(aes(xmin = 0.5693788, xmax = 1.5422579, ymin = -Inf, ymax = Inf),
            fill = "#AABBA3", alpha = 0.3, position = "identity") +
  geom_line(data = freq_Music_12m_PM10_df, aes(x = beatinfo, y = value, color=variable)) +
  labs(title = "Periodicity in PM10", x = "Period (beat)", y = "Power (a.u.)")+
  theme_minimal()+
  scale_color_manual(values=wes_palette(name="Moonrise3"))+
  geom_vline(xintercept = c(1, 2, 4, 8), color = "black", linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 2, 4, 8), limits = c(0.5, 10),trans = 'log2')

# Arrange the plots in a 5x2 grid
grid.arrange(grobs = list(PM1_plot,PM2_plot,PM3_plot,PM4_plot,PM5_plot,
                          PM6_plot,PM7_plot,PM8_plot,PM9_plot,PM10_plot), ncol = 2)

