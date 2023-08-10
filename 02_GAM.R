# Header ####
rm(list = ls())
library(ggplot2)
library(GGally)
library(mgcv)
library(gratia)
library(dplyr)
library(xtable)
library(tibble)
setwd('./')

# 0) Load data ####
# - - - - - - - -

input_folder = './'
graphs_path = './graphs/'

worst_concurvity <- function(m){
  concurvity_matrix <- concurvity(m, full=F)$worst
  # Remove any interaction terms
  concurvity_matrix <- concurvity_matrix[,!grepl('^ti[(]', colnames(concurvity_matrix))]
  concurvity_matrix <- concurvity_matrix[!grepl('^ti[(]', rownames(concurvity_matrix)),]
  # Remove concurvities of parameters with themselves
  concurvity_matrix <- concurvity_matrix[which(concurvity_matrix < 1)]
  return(max(concurvity_matrix))
}

concurvity_list <- function(model_list){
  conc_list <- c()
  for (m in model_list){
    conc_list <- append(conc_list, worst_concurvity(m))
  }
  return(conc_list)
}

r2_list <- function(model_list){
  r2_list <- c()
  for (m in model_list){
    r2_list <- append(r2_list, summary(m)$r.sq)
  }
  return(r2_list)
}

# 1) Load & prepare data ####
# - - - - - - - - - - - - - -

df <- read.csv(paste0(input_folder,'data_full.csv'))
df$timestamp = as.POSIXct(strptime(df$timestamp, '%Y-%m-%d %H:%M:%S'))
names(df)

temp <- df
# Select daytime from above-canopy PAR
temp = temp[which(temp$PAR_above_canopy > 50),]
# Exclude CO flux, leaf temperature and transpiration outliers
temp = temp[which(temp$TL < 45),]
temp = temp[which(temp$Tr >= 0),]
temp = temp[which(temp$co.flux > -10),]

temp$treatment = as.factor(temp$treatment)
temp$season = as.factor(temp$season)

# Remove incomplete rows as the GAM can't be estimated otherwise
temp = temp[complete.cases(temp[,c('co.flux', 'Tr', 'PAR', 'TL', 'SWC')]),]

temp$time <- as.numeric(strftime(temp$timestamp, '%H')) + as.numeric(strftime(temp$timestamp, '%M'))/60
temp$doy <- as.numeric(strftime(temp$timestamp, '%j'))

# 2) Data exploration ####
# - - - - - - - - - - - -

plt <- ggpairs(temp %>% select(c(co.flux, Tr, PAR, TL, VPD, SWC)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) )) +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt
ggsave(paste0(graphs_path, 'pairs.jpg'), width=18, height=18, units = "cm", dpi = 600)
ggsave(paste0(graphs_path, 'pairs.pdf'), width=18, height=18, units = "cm", dpi = 600)


plt <- ggpairs(temp[which(temp$treatment == 'Irrigated'),] %>% select(c(co.flux, Tr, PAR, TL, VPD, SWC)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) )) +
  theme_bw() +
  ggtitle('Irrigated') + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt
ggsave(paste0(graphs_path, 'pairs_irr.jpg'), width=18, height=18, units = "cm", dpi = 600)

plt <- ggpairs(temp[which(temp$treatment == 'Droughted'),] %>% select(c(co.flux, Tr, PAR, TL, VPD, SWC)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) )) +
  theme_bw() +
  ggtitle('Droughted') + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt
ggsave(paste0(graphs_path, 'pairs_ctr.jpg'), width=18, height=18, units = "cm", dpi = 600)

# 3) GAMs ####
# - - - - - - 
# For tutorials on GAMs, check:
# https://noamross.github.io/gams-in-r-course/chapter2
# https://towardsdatascience.com/producing-insights-with-generalized-additive-models-gams-cf2b68b1b847

# Leaf T and PAR
m0 <- gam(co.flux ~ s(PAR, by=treatment) + s(TL, by=treatment) + treatment, data=temp, method='REML')
m1 <- gam(co.flux ~ s(PAR) + s(TL), data=temp, method='REML')
#m1 <- gam(co.flux ~ s(PAR, k=3) + s(TL), data=temp, method='REML', family=scat(link="identity")) # Test T-distribution

# Adding transpiration to m0 & m1
m2 <- gam(co.flux ~ s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m3 <- gam(co.flux ~ s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Adding SWC to m2 is significant in gam.check(), so shouldn't be done
# Adding SWC to m3: same problem

# Add VPD to m2: VPD is significant in gam.check()
# Add VPD to m3: VPD is significant in gam.check()

# Add interactions:
# Adding to m2 & m3: PAR with TL
m4 <- gam(co.flux ~ ti(PAR, TL, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m5 <- gam(co.flux ~ ti(PAR, TL) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Adding to m2 & m3: PAR with Tr
m6 <- gam(co.flux ~ ti(PAR, Tr, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m7 <- gam(co.flux ~ ti(PAR, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Adding to m2 & m3: TL with Tr
m8 <- gam(co.flux ~ ti(TL, Tr, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m9 <- gam(co.flux ~ ti(TL, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Adding to m2 & m3: TL with Tr, removing PAR
m10 <- gam(co.flux ~ ti(TL, Tr, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m11 <- gam(co.flux ~ ti(TL, Tr) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Adding to m8 & m9: SWC is significant in gam.check()
# Any interaction with SWC is significant in gam.check()

# Adding to m8 & m9: TL with VPD
m12 <- gam(co.flux ~ ti(TL, VPD, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m13 <- gam(co.flux ~ ti(TL, VPD) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)
# Adding to m8 & m9: Tr with VPD
m14 <- gam(co.flux ~ ti(Tr, VPD, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m15 <- gam(co.flux ~ ti(Tr, VPD, k=3) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Add 3-way interactions:
# m2 & m3 with interaction of PAR, TL and Tr
m16 <- gam(co.flux ~ ti(PAR, TL, Tr, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment), data=temp, method='REML', select=T)
m17 <- gam(co.flux ~ ti(PAR, TL, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# m2 & m3 with interaction of PAR, TL and Tr, removing PAR (only increases AIC slightly)
m18 <- gam(co.flux ~ ti(PAR, TL, Tr, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m19 <- gam(co.flux ~ ti(PAR, TL, Tr) + s(TL) + s(Tr), data=temp, method='REML', select=T)

# Other interactions
# m8 & m9 with interaction of VPD, TL and Tr
m20 <- gam(co.flux ~ ti(VPD, TL, Tr, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + s(PAR, by=treatment) + treatment, data=temp, method='REML', select=T)
m21 <- gam(co.flux ~ ti(VPD, TL, Tr) + s(TL) + s(Tr) + s(PAR), data=temp, method='REML', select=T)

# m8 & m9, with VPD instead of TL
m22 <- gam(co.flux ~ ti(VPD, Tr, by=treatment) + s(PAR, by=treatment) + s(VPD, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m23 <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(VPD) + s(Tr), data=temp, method='REML', select=T)

# m8 & m9, with VPD instead of TL, except in the s()
m24 <- gam(co.flux ~ ti(VPD, Tr, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment) + treatment, data=temp, method='REML', select=T)
m25 <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp, method='REML', select=T)

# Final model: m27
m26 <- gam(co.flux ~ ti(VPD, Tr, k=4, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment), data=temp, method='REML', select=T)
m27 <- gam(co.flux ~ ti(VPD, Tr, k=4) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T)
summary(m27)
gam.check(m27, rep=500)
plot(m27, scheme = 1, pages=1, too.far = 0.05)
vis.gam(m27, view = c('VPD','Tr'), plot.type = "contour", too.far = 0.05)

# All
comparison_df <- AIC(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,
                     m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$interact <- 0
comparison_df[which(comparison_df$model %in%
                      c('m4','m5','m6','m7','m8','m9','m10','m11','m12','m13','m14','m15','m22','m23','m24','m25','m26','m27')),]$interact <- 2
comparison_df[which(comparison_df$model %in%
                      c('m16','m17','m18','m19','m20','m21')),]$interact <- 3
comparison_df$r2 <- r2_list(list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,
                                 m17, m18, m19, m20, m21, m22, m23, m24, m25, m26,m27))
comparison_df$concurvity <- concurvity_list(list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16,
                                 m17, m18, m19, m20, m21, m22, m23, m24, m25, m26,m27))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# by plot
comparison_df <- AIC(m0, m2, m4, m6, m8, m10, m12, m14, m16,
                     m18, m20, m22, m24, m26)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$interact <- 0
comparison_df[which(comparison_df$model %in%
                      c('m4','m5','m6','m7','m8','m9','m10','m11','m12','m13','m14','m15','m22','m23','m24','m25','m26','m27')),]$interact <- 2
comparison_df[which(comparison_df$model %in%
                      c('m16','m17','m18','m19','m20','m21')),]$interact <- 3
comparison_df$r2 <- r2_list(list(m0, m2, m4, m6, m8, m10, m12, m14, m16,
                                 m18, m20, m22, m24, m26))
comparison_df$concurvity <- concurvity_list(list(m0, m2, m4, m6, m8, m10, m12, m14, m16,
                                                 m18, m20, m22, m24, m26))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# Complete, not by plot
comparison_df <- AIC(m1, m3, m5, m7, m9, m11, m13, m15, 
                     m17, m19, m21, m23, m25, m27)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$interact <- 0
comparison_df[which(comparison_df$model %in%
                      c('m4','m5','m6','m7','m8','m9','m10','m11','m12','m13','m14','m15','m22','m23','m24','m25','m26','m27')),]$interact <- 2
comparison_df[which(comparison_df$model %in%
                      c('m16','m17','m18','m19','m20','m21')),]$interact <- 3
comparison_df$r2 <- r2_list(list(m1, m3, m5, m7, m9, m11, m13, m15, m17, m19, m21, m23, m25, m27))
comparison_df$concurvity <- concurvity_list(list(m1, m3, m5, m7, m9, m11, m13, m15, m17, m19, m21, m23, m25, m27))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# 4) Rename models & save  ####
# - - - - - - - - - - - - - - -

# 4a) Latex table ####
# - - - - - - - - - - -
m1 <- gam(co.flux ~ s(PAR) + s(TL), data=temp, method='REML') # m1
m2 <- gam(co.flux ~ s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m3
m3 <- gam(co.flux ~ ti(PAR, TL) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m5
m4 <- gam(co.flux ~ ti(PAR, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m7
m5 <- gam(co.flux ~ ti(TL, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m9
m6 <- gam(co.flux ~ ti(TL, Tr) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m11
m7 <- gam(co.flux ~ ti(TL, VPD) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m13
m8 <- gam(co.flux ~ ti(Tr, VPD, k=3) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m15 
m9 <- gam(co.flux ~ ti(PAR, TL, Tr) + s(PAR) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m17
m10 <- gam(co.flux ~ ti(PAR, TL, Tr) + s(TL) + s(Tr), data=temp, method='REML', select=T) # m19
m11 <- gam(co.flux ~ ti(VPD, TL, Tr) + s(TL) + s(Tr) + s(PAR), data=temp, method='REML', select=T) # m21
m12 <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp, method='REML', select=T) # m25
m13 <- gam(co.flux ~ ti(VPD, Tr, k=4, by=treatment) + s(PAR, by=treatment) + s(TL, by=treatment) + s(Tr, by=treatment), data=temp, method='REML', select=T) # m26
m14 <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL), data=temp, method='REML', select=T)
# Create model comparison after renaming
comparison_df <- AIC(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$interact <- 0
comparison_df[which(comparison_df$model %in%
                      c('m3','m4','m5','m6','m7','m8','m12','m13','m14')),]$interact <- 2
comparison_df[which(comparison_df$model %in%
                      c('m9','m10','m11')),]$interact <- 3
comparison_df$r2 <- r2_list(list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14))
comparison_df$concurvity <- concurvity_list(list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# Save as latex
print.xtable(xtable(comparison_df, digits=c(0,0,2,0,0,2,2)), file = paste0(graphs_path, 'model_comparison.tex'),
             include.rownames=F)

# 5) Final model diagnostics plots ####
# - - - - - - - - - - - - - - - - - - -

# Run model
m_final <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL), data=temp,
               method='REML', select=T)
#m_final <- gam(co.flux ~ ti(TL, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp,
#               method='REML', select=T) # Alternative with additive Tr
summary(m_final)
gam.check(m_final)

# Q-Q plot
plt = ggplot(data.frame(Standardized_Residuals = residuals(m_final)),
                        aes(sample=Standardized_Residuals))
plt = plt + stat_qq(alpha = 0.25, size=0.25)
plt = plt + geom_abline(slope = 1, intercept = 0, color = "#808080", linetype = "dashed")
plt = plt + labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Standardized Residuals")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
qq_plot <- plt

# Histogram
plt = ggplot(data.frame(Standardized_Residuals = residuals(m_final)),
             aes(x=Standardized_Residuals))
plt = plt + geom_histogram(binwidth = 0.5, color = "black", fill = "#808080", size=0.25)
plt = plt + labs(title = "Histogram of Standardized Residuals",
       x = "Standardized Residuals",
       y = "Frequency")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
hist_plot <- plt

# Response
plot_df <- data.frame(Fitted_Values = fitted(m_final), Observed_Response = temp[which(!is.na(temp$VPD)),]$co.flux)
plt = ggplot(plot_df, aes(x = Fitted_Values, y = Observed_Response))
plt = plt + geom_point(color = "black", alpha = 0.25, size=0.25)
plt = plt + geom_abline(slope = 1, intercept = 0, color = "#808080", linetype = "dashed")
plt = plt + labs(title = "Response vs Fitted Values",
       x = "Fitted Values",
       y = "Observed Response")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
response_plot <- plt

# Residuals
plot_df <- data.frame(Linear_Predictor = fitted(m_final), Residuals = residuals(m_final))
plt = ggplot(plot_df, aes(x = Linear_Predictor, y = Residuals))
plt = plt + geom_point(color = "black", alpha = 0.25, size=0.25)
plt = plt + geom_hline(yintercept = 0, color = "#808080", linetype = "dashed")
plt = plt + labs(title = "Residuals vs Linear Predictor",
       x = "Linear Predictor",
       y = "Residuals")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
residuals_plot <- plt

# Combine them
library(gridExtra)
combined_plot <- grid.arrange(qq_plot, hist_plot, response_plot, residuals_plot,
                              ncol = 2, nrow = 2)
print(combined_plot)
ggsave(paste0(graphs_path, 'model_stats.jpg'), plot = combined_plot, width = 18, height = 12, units = "cm", dpi=600)
ggsave(paste0(graphs_path, 'model_stats.pdf'), plot = combined_plot, width = 18, height = 12, units = "cm", dpi=600)

# 5a) Output plots ####
# - - - - - - - - - - -

# Obtain all the subplots
predictor_names <- attr(terms(m_final), "term.labels")

plot_width      <- 8 # cm
plot_height     <- 8 # cm
plot_resolution <- 1200 # dpi
for(plot_nb in 2:length(predictor_names)){
  # Save as jpg
  jpeg(paste0(graphs_path, 'model_output_', plot_nb,'.jpg'), 
       width = plot_width, height = plot_height, units='cm', res = plot_resolution)
  par(mar = c(2.5,2.5,1,1))
  plot(m_final, select=plot_nb, scheme = 1, too.far = 0.05, family='serif', mgp=c(1.2,0.4,0))
  dev.off()
  # Save as PDF
  pdf(paste0(graphs_path, 'model_output_', plot_nb,'.pdf'), 
      width = plot_width/2.54, height = plot_height/2.54)
  par(mar = c(1,1,1,1))
  plot(m_final, select=plot_nb, scheme = 1, too.far = 0.05, family='serif')
  dev.off()
}

dev.off()
# Save as jpg
jpeg(paste0(graphs_path, 'model_output_1.jpg'), 
     width = plot_width, height = plot_height, units='cm', res = plot_resolution)
par(mar = c(1,1,1,1))
plot(m_final, select=1, scheme = 1, too.far = 0.05, family='serif')
text(x = 0.1, y = 0.9, '(a)', cex = 1.5, col = "black", font = 2)
dev.off()
# Save as PDF
pdf(paste0(graphs_path, 'model_output_1.pdf'), 
     width = plot_width/2.54, height = plot_height/2.54)
par(mar = c(1,1,1,1))
plot(m_final, select=1, scheme = 1, too.far = 0.05, family='serif')
dev.off()

# COMBINED
# - - - - -
# PDF
pdf(paste0(graphs_path, 'model_output_all.pdf'), 
    width = plot_width/2.54, height = plot_height*3/2.54)
# Layout
layout(matrix(c(1, 2, 3), nrow = 3,
              ncol = 1, byrow = TRUE))
# interaction
par(mar = c(1,1,1,1))
plot(m_final, select=1, scheme = 1, too.far = 0.05, family='serif')
# 3 other parameters
for(plot_nb in 2:length(predictor_names)){
  par(mar = c(2.5,2.5,1,1))
  plot(m_final, select=plot_nb, scheme = 1, too.far = 0.05, family='serif', mgp=c(1.2,0.4,0))
}
dev.off()

# JPG
jpeg(paste0(graphs_path, 'model_output_all.jpg'), 
     width = plot_width, height = plot_height*3, units='cm', res = plot_resolution)
# Layout
layout(matrix(c(1, 2, 3), nrow = 3,
              ncol = 1, byrow = TRUE))
# interaction
par(mar = c(1,1,1,1))
plot(m_final, select=1, scheme = 1, too.far = 0.05, family='serif')
# 3 other parameters
for(plot_nb in 2:length(predictor_names)){
  par(mar = c(2.5,2.5,1,1))
  plot(m_final, select=plot_nb, scheme = 1, too.far = 0.05, family='serif', mgp=c(1.2,0.4,0))
}
dev.off()

# Interactions
# - - - - - - -
# JPG
jpeg(paste0(graphs_path, 'model_output_interaction.jpg'), 
     width = plot_width, height = plot_height, units='cm', res = plot_resolution)
par(mar = c(2.5,2.5,2,1))
vis.gam(m_final, view = c('VPD','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0))
dev.off()
# PDF
pdf(paste0(graphs_path, 'model_output_interaction.pdf'), 
     width = plot_width/2.54, height = plot_height/2.54)
par(mar = c(2.5,2.5,2,1))
vis.gam(m_final, view = c('VPD','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0))
dev.off()

# 5b) Plots of parameter interactions ####
# - - - - - - - - - - - - - - - - - - - - 
# m_final <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp,
#                method='REML', select=T)

plot_width      <- 8 # cm
plot_height     <- 8 # cm
plot_resolution <- 1200 # dpi

# JPG
jpeg(paste0(graphs_path, 'model_output_interactions.jpg'), 
     width = plot_width*2, height = plot_height*2, units='cm', res = plot_resolution)
# Layout
layout(matrix(c(1, 2, 3, 4), nrow = 2,
              ncol = 2, byrow = TRUE))
par(mar = c(2.5,2.5,2,1))

vis.gam(m_final, view = c('VPD','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(a) VPD vs. Tr')

vis.gam(m_final, view = c('TL','PAR'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(b) TL vs. PAR')

vis.gam(m_final, view = c('TL','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(c) TL vs. Tr')

vis.gam(m_final, view = c('Tr','PAR'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(d) Tr vs. PAR')

dev.off()

# JPG
pdf(paste0(graphs_path, 'model_output_interactions.pdf'), 
     width = plot_width*2/2.54, height = plot_height*2/2.54)
# Layout
layout(matrix(c(1, 2, 3, 4), nrow = 2,
              ncol = 2, byrow = TRUE))
par(mar = c(2.5,2.5,2,1))

vis.gam(m_final, view = c('VPD','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(a) VPD vs. Tr')

vis.gam(m_final, view = c('TL','PAR'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(b) TL vs. PAR')

vis.gam(m_final, view = c('TL','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(c) TL vs. Tr')

vis.gam(m_final, view = c('Tr','PAR'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0), main='(d) Tr vs. PAR')

dev.off()
