# Header ####
rm(list = ls())
library(ggplot2)
library(GGally) # Comparing variables and data exploration
library(mgcv) # Library to fit gams
library(gratia) # Modern library for visualizing and assessing gams
library(visreg)
library(dplyr)
library(tidyr)
library(xtable)
library(tibble)
setwd('./')

# 0) Load data ####
# - - - - - - - -

input_folder = './'
graphs_path = 'graphs/'

# Colours for each treatment, suitable for colourblindness
cbPalette <- c("#000000", "#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

# 1) Loading data ####
# - - - - - - - - - -
# https://towardsdatascience.com/producing-insights-with-generalized-additive-models-gams-cf2b68b1b847

temp <- read.csv(paste0(input_folder,'all_summarised_daily.csv'))

# Rename columns for easier use
names(temp)[names(temp) == 'flux.co.ch_oc.nmol_m2_s']        <- 'co.flux'
names(temp)[names(temp) == 'flux.h2o.ch_oc.mmol_m2_s']       <- 'Tr'
names(temp)[names(temp) == 'flux.co2.ch_oc.umol_m2_s1']      <- 'co2.flux'
names(temp)[names(temp) == 'conc.h2o.mmol_mol.oc']           <- 'H2O'
names(temp)[names(temp) == 'par.current.chamber.umol_m2_s1'] <- 'PAR'
names(temp)[names(temp) == 'temp.leaf.current.chamber.c.oc'] <- 'TL'
names(temp)[names(temp) == 'VPD.Pa.oc']                      <- 'VPD'
names(temp)[names(temp) == 'swc_mean']                       <- 'swc'
names(temp)[names(temp) == 'conc_ci.co.nmol_mol']            <- 'COi'
names(temp)[names(temp) == 'conc_ci.h2o.mmol_mol']           <- 'H2Oi'

temp = temp[complete.cases(temp[,c('co.flux', 'Tr', 'PAR', 'TL', 'VPD')]),]
temp$timestamp <- as.POSIXct(strptime(temp$timestamp, format='%Y-%m-%d'))

# 2) Data exploration ####
# - - - - - - - - - - - -

# Used in model development, shown in Supplementary Figure SM8.1

plt <- ggpairs(temp %>% select(c(co.flux, Tr, PAR, TL, VPD, swc, g_tCO, COi)),
               columnLabels = c('f[CO]', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC', 'g["t,CO"]', 'c["CO,i"]'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) ),
               ggplot2::aes(colour = temp$treatment),
               upper = list(continuous = wrap("cor", size = 2.5))) +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt = plt + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
plt
ggsave(paste0(graphs_path, 'pairs_summarised.jpg'), width=18, height=18, units = "cm", dpi = 600)
ggsave(paste0(graphs_path, 'pairs_summarised.pdf'), width=18, height=18, units = "cm", dpi = 600)

# 3) GAMs ####
# - - - - - - 
# Tutorial: https://noamross.github.io/gams-in-r-course/chapter2

# Prepare categorical columns
testing_dataset <- temp
testing_dataset$treatment <- as.factor(testing_dataset$treatment)
testing_dataset$season <- as.factor(testing_dataset$season)

# Remove outliers: First mark outliers
testing_dataset$prob <- ''
testing_dataset[which((testing_dataset$VPD > 2000) & (testing_dataset$co.flux < -0.5) & (testing_dataset$treatment == 'irr')),]$prob  <- 'b' #'bad irr'
testing_dataset[which((testing_dataset$co.flux < -0.7) & (testing_dataset$treatment == 'irr')),]$prob  <- 'b' #'bad irr'
# Now remove outliers
testing_dataset <- testing_dataset[which(testing_dataset$prob != 'b'),]
# Remove negative calculate leaf internal CO concentration
testing_dataset[which(testing_dataset$COi < 0),]$COi <- NA

# Show pairwise dataset
plt <- ggpairs(testing_dataset %>% select(c(co.flux, Tr, PAR, TL, VPD, swc, g_tCO, COi)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC', 'g["t,CO"]', 'CO~int.'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) ),
               ggplot2::aes(colour = testing_dataset$treatment),
               upper = list(continuous = wrap("cor", size = 2.5))) +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt = plt + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
plt

# Production & transport in 2 separate models
m_coi <- gam(COi ~ s(TL, k=4, by=treatment) + s(PAR, k=4, by=treatment) + treatment, # Production GOOD
          data=testing_dataset, method='REML', select=T)
summary(m_coi)
AIC(m_coi)
concurvity(m_coi, full=F)$worst
gam.check(m_coi)

# COi vs PAR
p = visreg(m_coi,"PAR", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi, xvar = 'PAR', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(PAR, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("PAR [",mu,mol~m^{-2}~s^{-1},"]")),
                 y = expression(paste(c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste('PAR contribution to ','c'['CO,i'])))
plt_coi_par <- plt
plt_coi_par
ggsave(paste0(graphs_path, 'model2a_output_PAR.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model2a_output_PAR.pdf'), width=8, height=8, units = "cm", dpi = 1200)

# COi vs TL
p = visreg(m_coi,"TL", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi, xvar = 'TL', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(TL, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste('T'[L]," [°C]")),
                 y = expression(paste(c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste('T'['L']," contribution to ",'c'['CO,i'])))
plt_coi_tl <- plt
plt_coi_tl
ggsave(paste0(graphs_path, 'model2b_output_TL.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model2b_output_TL.pdf'), width=8, height=8, units = "cm", dpi = 1200)

# Flux model
library(gratia)
library(metR)

# Overall model
m_cof <- gam(co.flux ~ te(Tr, g_tCO, by=treatment, k=4),
                 data=testing_dataset, method='REML', select=T)
summary(m_cof)
AIC(m_cof)

# Droughted treatment
m_cof_dro <- gam(co.flux ~ te(Tr, g_tCO, k=4),
                 data=testing_dataset[which((testing_dataset$treatment == 'dro')),], method='REML', select=T)
plt_cof_dro <- draw(m_cof_dro, contour = T, n = 50, select=1)
plt_cof_dro[[1]]$labels$caption <- NULL # Remove
plt_cof_dro[[1]]$labels$title <- NULL # Remove
plt_cof_dro <- plt_cof_dro + scale_fill_gradient2(low='#2166AC', mid='#eaeaea', high='#9c1a16', na.value="#ffffff00",
                                                  limits = c(-1.5, 4))
plt_cof_dro <- plt_cof_dro + metR::geom_text_contour(aes(z = est), stroke = 0.15)
plt_cof_dro <- plt_cof_dro + theme_bw()
plt_cof_dro <- plt_cof_dro + theme(legend.position = "right",
                                   text=element_text(family="serif"),
                                   plot.title = element_text(hjust = 0.5))
plt_cof_dro <- plt_cof_dro + ggtitle(expression(paste("Droughted: ",'g'['t,CO'],' & Tr contrib. to f'['CO'])))
plt_cof_dro <- plt_cof_dro + labs(y=expression(paste(g[paste("t,CO")]," [mol ",mol^{-2}," ",s^{-1},"]")),
                                  x=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))
plt_cof_dro
ggsave(paste0(graphs_path, 'model3b_output_interaction_H2Oi_Tr.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model3b_output_interaction_H2Oi_Tr.pdf'), width=8, height=8, units = "cm", dpi = 1200)


# Irrigated treatment
m_cof_irr <- gam(co.flux ~ te(Tr, g_tCO, k=4),
                 data=testing_dataset[which((testing_dataset$treatment == 'irr')),], method='REML', select=T)
plt_cof_irr <- draw(m_cof_irr, contour = T, n = 50, select=1)
plt_cof_irr[[1]]$labels$caption <- NULL # Remove
plt_cof_irr[[1]]$labels$title <- NULL # Remove
plt_cof_irr <- plt_cof_irr + scale_fill_gradient2(low='#2166AC', mid='#eaeaea', high='#9c1a16', na.value="#ffffff00",
                                                  limits = c(-1.5, 4))
plt_cof_irr <- plt_cof_irr + metR::geom_text_contour(aes(z = est), stroke = 0.15)
plt_cof_irr <- plt_cof_irr + theme_bw()
plt_cof_irr <- plt_cof_irr + theme(legend.position = "right",
                                   text=element_text(family="serif"),
                                   plot.title = element_text(hjust = 0.5))
plt_cof_irr <- plt_cof_irr + ggtitle(expression(paste("Irrigated: ",'g'['t,CO'],' & Tr contrib. to f'['CO'])))
plt_cof_irr <- plt_cof_irr + labs(y=expression(paste(g[paste("t,CO")]," [mol ",mol^{-2}," ",s^{-1},"]")),
                                  x=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))

plt_cof_irr
ggsave(paste0(graphs_path, 'model3c_output_interaction_H2Oi_Tr.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model3c_output_interaction_H2Oi_Tr.pdf'), width=8, height=8, units = "cm", dpi = 1200)

# Combine plots
library(ggpubr)
library(cowplot)

legend <- get_legend(plt_cof_irr)
plt_cof_irr2 <- plt_cof_irr + theme(legend.position = "none")
plt_cof_dro2 <- plt_cof_dro + theme(legend.position = "none")

plts <- ggarrange(plt_coi_par, plt_coi_tl,
          plt_cof_irr2, plt_cof_dro2,
          labels = c("(a)", "(b)",
                     "(c)", "(d)"),
          ncol = 2, nrow = 2)

final_plot <- plot_grid(plts, legend, ncol = 2, rel_widths = c(1, 0.1))
final_plot
ggsave(paste0(graphs_path, 'model_output_all.jpg'), width = 21, height = 20, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model_output_all.pdf'), width = 21, height = 20, units = "cm", dpi = 1200)

# Supplementary figures ####
# - - - - - - - - - - - - - -

#  Sup: COi vs VPD (Fig. SM8.3) ####
# - - - - - - - - - - - - - - - - - -
m_coi_vpd <- gam(COi ~ s(VPD, k=4, by=treatment) + s(PAR, k=4, by=treatment) + treatment, # Production GOOD
                 data=testing_dataset, method='REML', select=T)
summary(m_coi_vpd)
AIC(m_coi_vpd)


p = visreg(m_coi_vpd,"PAR", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi_vpd, xvar = 'PAR', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(PAR, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("PAR [",mu,mol~m^{-2}~s^{-1},"]")),
                 y = expression(paste('Part. contribution to  ',c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
#plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste('PAR contribution to ','c'['CO,i'])))
plt_coi2_par <- plt
plt_coi2_par
ggsave(paste0(graphs_path, 'supplement/04_model2a_output_PAR.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'supplement/04_model2a_output_PAR.pdf'), width=8, height=8, units = "cm", dpi = 1200)

p = visreg(m_coi_vpd,"VPD", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi_vpd, xvar = 'VPD', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(VPD, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("VPD [Pa]")),
                 y = expression(paste('Part. contribution to  ',c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste("VPD contribution to ",'c'['CO,i'])))
plt_coi2_vpd <- plt
plt_coi2_vpd
ggsave(paste0(graphs_path, 'supplement/04_model2b_output_VPD.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'supplement/04_model2b_output_VPD.pdf'), width=8, height=8, units = "cm", dpi = 1200)

ggarrange(plt_coi2_par, plt_coi2_vpd, 
          labels = c("(a)",  "(b)"),
          ncol = 2, nrow = 1)
ggsave(paste0(graphs_path, 'supplement/04_model_output_all.jpg'), width = 20, height = 8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'supplement/04_model_output_all.pdf'), width = 20, height = 8, units = "cm", dpi = 1200)


# Sup: COi model comparison, Supplementary Methods S8 ####
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

# Table SM8.1
comparison_df <- AIC(m_coi, m_coi_vpd)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$r2 <- r2_list(list(m_coi, m_coi_vpd))*100
comparison_df$concurvity <- concurvity_list(list(m_coi, m_coi_vpd))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# Save as latex
#print.xtable(xtable(comparison_df, digits=c(0,2,2,0,1,2)), file = paste0(graphs_path, 'model_comparison_coi.tex'),
#             include.rownames=F)

# Sup: fCO model comparison ####
# - - - - - - - - - - - - - - -

# Table SM8.2

m3 <- gam(co.flux ~ te(Tr, g_tCO, by=treatment, k=4),
             data=testing_dataset, method='REML', select=T)
m4 <- gam(co.flux ~ te(g_tCO, Tr, by=treatment, k=4) + s(COi, by=treatment, k=4) + treatment,
             data=testing_dataset, method='REML', select=T)
m5 <- gam(co.flux ~ s(COi, by=treatment, k=4) + s(g_tCO, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m6 <- gam(co.flux ~ te(COi, g_tCO, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m7 <- gam(co.flux ~ te(COi, Tr, by=treatment, k=4) + s(g_tCO, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m8 <- gam(co.flux ~ te(COi, Tr, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m9 <- gam(co.flux ~ s(COi, by=treatment, k=4) + s(Tr, by=treatment, k=4) + s(COi, by=treatment, k=4) + treatment,
             data=testing_dataset, method='REML', select=T)
m10 <- gam(co.flux ~ s(g_tCO, by=treatment, k=4) + s(COi, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m11 <- gam(co.flux ~ s(g_tCO, by=treatment, k=4) + treatment, # Chosen model
          data=testing_dataset, method='REML', select=T)
m12 <- gam(co.flux ~ s(g_tCO, by=treatment, k=4) + s(COi, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
m13 <- gam(co.flux ~ s(g_tCO, by=treatment, k=4) + s(g_tCO, by=treatment, k=4) + treatment,
           data=testing_dataset, method='REML', select=T)
AIC(m13)
summary(m13)
concurvity(m3, full=F)$worst

comparison_df <- AIC(m3,m4, m5, m6,m7,m8,m9,m10,m11,m12,m13)
comparison_df <- rownames_to_column(comparison_df)
names(comparison_df)[names(comparison_df) == 'rowname'] <- 'model'
comparison_df$r2 <- r2_list(list(m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13))*100
comparison_df$concurvity <- concurvity_list(list(m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13))
comparison_df <- comparison_df[order(comparison_df$AIC, decreasing = F), ]
comparison_df

# Save as latex
print.xtable(xtable(comparison_df, digits=c(0,2,2,0,1,2)), file = paste0(graphs_path, 'supplement/model_comparison_cof.tex'),
             include.rownames=F)
