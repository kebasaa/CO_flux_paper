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

# Colours for each treatment
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

visreg2 <- function(mod, xvar, by, data){
  plt_plot <- visreg(mod, xvar = xvar, by=by, data = data, gg=T, method = "REML", overlay=T)
  plt_line <- visreg(mod, xvar = xvar, data = data, gg=T, method = "REML", overlay=T)
  plt_plot_data = ggplot_build(plt_plot)$data[[2]] # points and colours
  plt_line_data = ggplot_build(plt_line)$data[[3]] # line
  plt_std_data = ggplot_build(plt_line)$data[[1]] %>% # stddev around the line
    group_by(x) %>%
    mutate(y_type = if_else(y == min(y), "ymin", "ymax")) %>%
    spread(key = y_type, value = y)
  yminval = min(plt_plot_data$y)- (max(plt_plot_data$y) - min(plt_plot_data$y))/20
  plt_plot_data$grp <- NA
  for (grp_num in 1:length(levels(plt_plot$plot_env$bb))){
    current_label = levels(plt_plot$plot_env$bb)[grp_num]
    current_colour = plt_plot$plot_env$col[grp_num]
    plt_plot_data[which(plt_plot_data$colour == current_colour),]$grp <- current_label
  }
  plt = ggplot()
  plt = plt + geom_point(data=plt_plot_data, aes(x=x, y=y, color = grp), size=0.75)
  plt = plt + geom_ribbon(data=plt_std_data, aes(x=x, ymin=ymin, ymax=ymax), fill = 'gray85', alpha=0.5)
  plt = plt + geom_line(data=plt_line_data, aes(x=x, y=y), color = 'black')
  plt = plt + theme_bw()
  plt = plt + labs(x=plt_plot$labels$x, y=plt_plot$labels$y)
  #plt = plt + labs(x='Tr', y='CO flux', colour='Treatment')
  return(plt)
}

visreg3 <- function(mod, xvar, by, data, minrange=0.25, maxrange=0.75){
  plt_plot <- visreg(mod, xvar = xvar, by=by, data = data, gg=T, method = "REML", overlay=T)
  plt_line <- visreg(mod, xvar = xvar, data = data, gg=T, method = "REML", overlay=T)
  plt_plot_data = ggplot_build(plt_plot)$data[[2]] # points and colours
  plt_line_data = ggplot_build(plt_line)$data[[3]] # line
  plt_std_data = ggplot_build(plt_line)$data[[1]] %>% # stddev around the line
    group_by(x) %>%
    mutate(y_type = if_else(y == min(y), "ymin", "ymax")) %>%
    spread(key = y_type, value = y)
  yminval = min(ggplot_build(plt_line)$data[[1]]$y)
  plt_plot_data$grp <- NA
  plt_line_data$grp <- NA
  for (grp_num in 1:length(levels(plt_plot$plot_env$bb))){
    current_label = levels(plt_plot$plot_env$bb)[grp_num]
    current_colour = plt_plot$plot_env$col[grp_num]
    plt_plot_data[which(plt_plot_data$colour == current_colour),]$grp <- current_label
    # Create quantile range df
    xmin = quantile(plt_plot_data[which(plt_plot_data$group == grp_num),]$x, minrange)
    xmax = quantile(plt_plot_data[which(plt_plot_data$group == grp_num),]$x, maxrange)
    # Find median and corresponding y value for median
    medianxval = median(plt_plot_data[which(plt_plot_data$group == grp_num),]$x)
    next_higher <- min(plt_line_data[which(plt_line_data$x > medianxval),]$y, na.rm = T)
    next_lower  <- max(plt_line_data[which(plt_line_data$x < medianxval),]$y, na.rm = T)
    medianyval = (next_higher + next_lower)/2
    if(!exists('quantile_df')){
      quantile_df <- plt_line_data[which((plt_line_data$x >= xmin) & (plt_line_data$x < xmax)),]
      quantile_df$grp <- current_label
      quantile_df = quantile_df[which(!is.na(quantile_df$grp)),]
    }else{
      temp <- plt_line_data[which((plt_line_data$x >= xmin) & (plt_line_data$x < xmax)),]
      temp$grp <- current_label
      temp = temp[which(!is.na(temp$grp)),]
      quantile_df = rbind(quantile_df, temp)
    }
    if(!exists('median_df')){
      median_df = data.frame(grp = current_label, x = medianxval, y = medianyval)
    }else{
      median_df = rbind(median_df, data.frame(grp = current_label, x = medianxval, y = medianyval))
    }
  }
  print(median_df)
  plt = ggplot()
  plt = plt + geom_ribbon(data=quantile_df, aes(x=x, ymin=yminval, ymax=y, fill=grp), alpha=0.5)
  plt = plt + geom_ribbon(data=plt_std_data, aes(x=x, ymin=ymin, ymax=ymax), fill = 'gray85', alpha=0.5)
  plt = plt + geom_line(data=plt_line_data, aes(x=x, y=y), color = 'black')
  plt = plt + geom_point(data=median_df, aes(x=x, y=y, colour=grp), size=2)
  plt = plt + theme_bw()
  plt = plt + labs(x=plt_plot$labels$x, y=plt_plot$labels$y)
  return(plt)
}

# 1) Preparing data ####
# - - - - - - - - - - - -
# https://towardsdatascience.com/producing-insights-with-generalized-additive-models-gams-cf2b68b1b847

df <- read.csv(paste0(input_folder,'all_summarised_daily.csv'))

temp <- df
#temp = temp[which(temp$flux.h2o.ch_oc.mmol_m2_s >= 0),] # NEWLY REMOVED 20 MAY 2024

# Direction of flux
# temp$dir <- 'emission'
# temp[which(temp$flux.co.ch_oc.nmol_m2_s < 0.0),]$dir <- 'uptake'
# temp$dir <- factor(temp$dir)

# Rename columns
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

plt <- ggpairs(temp %>% select(c(co.flux, Tr, PAR, TL, VPD, swc, g_tCO, COi)),
               columnLabels = c('f[CO]', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC', 'g["t,CO"]', 'c["CO,i"]'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) ),
               ggplot2::aes(colour = temp$treatment),
               upper = list(continuous = wrap("cor", size = 2.5))) +
  #labs(subtitle = "Numeric variable exploration") +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt = plt + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
plt
ggsave(paste0(graphs_path, 'pairs_summarised.jpg'), width=18, height=18, units = "cm", dpi = 600)
ggsave(paste0(graphs_path, 'pairs_summarised.pdf'), width=18, height=18, units = "cm", dpi = 600)

# 3) GAMs ####
# - - - - - - 
# https://noamross.github.io/gams-in-r-course/chapter2

# Prep
# temp2 <- temp
# temp2$treatment <- as.factor(temp2$treatment)
# temp2$season <- as.factor(temp2$season)

testing_dataset <- temp
testing_dataset$treatment <- as.factor(testing_dataset$treatment)
testing_dataset$season <- as.factor(testing_dataset$season)

#testing_dataset <- temp2
testing_dataset$prob <- ''
testing_dataset[which((testing_dataset$VPD > 2000) & (testing_dataset$co.flux < -0.5) & (testing_dataset$treatment == 'irr')),]$prob  <- 'b' #'bad irr'
testing_dataset[which((testing_dataset$co.flux < -0.7) & (testing_dataset$treatment == 'irr')),]$prob  <- 'b' #'bad irr'

testing_dataset <- testing_dataset[which(testing_dataset$prob != 'b'),]

testing_dataset[which(testing_dataset$COi < 0),]$COi <- NA

mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + te(TL, Tr, by=treatment) + treatment,
          data=testing_dataset, method='REML', select=T)
AIC(mX)
summary(mX)
concurvity(mX, full=F)$worst
#plot(mX)

plt <- ggpairs(testing_dataset %>% select(c(co.flux, Tr, PAR, TL, VPD, swc, g_tCO, COi)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'SWC', 'g["t,CO"]', 'CO~int.'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) ),
               ggplot2::aes(colour = testing_dataset$treatment),
               upper = list(continuous = wrap("cor", size = 2.5))) +
  #labs(subtitle = "Numeric variable exploration") +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt = plt + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
plt

# Production & transport in 2 models!
m_coi <- gam(COi ~ s(TL, k=4, by=treatment) + s(PAR, k=4, by=treatment) + treatment, # Production GOOD
          data=testing_dataset, method='REML', select=T)
summary(m_coi)
AIC(m_coi)
concurvity(m_coi, full=F)$worst
gam.check(m_coi)
p = visreg(m_coi,"PAR", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi, xvar = 'PAR', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(PAR, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
#plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
#                     size = 5, hjust = 0, vjust = 1, family='serif')
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("PAR [",mu,mol~m^{-2}~s^{-1},"]")),
                 y = expression(paste(c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
#plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste('PAR contribution to ','c'['CO,i'])))
plt_coi_par <- plt
plt_coi_par
ggsave(paste0(graphs_path, 'model2a_output_PAR.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model2a_output_PAR.pdf'), width=8, height=8, units = "cm", dpi = 1200)

p = visreg(m_coi,"TL", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(m_coi, xvar = 'TL', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(TL, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
#plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
#                     size = 5, hjust = 0, vjust = 1, family='serif')
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste('T'[L]," [°C]")),
                 y = expression(paste(c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
#plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
                  plot.title = element_text(hjust = 0.5))
plt = plt + ggtitle( expression(paste('T'['L']," contribution to ",'c'['CO,i'])))
plt_coi_tl <- plt
plt_coi_tl
ggsave(paste0(graphs_path, 'model2b_output_TL.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model2b_output_TL.pdf'), width=8, height=8, units = "cm", dpi = 1200)

# FLUX
#m_cof <- gam(co.flux ~ te(H2Oi, Tr, by=treatment, k=4) + treatment,
#            data=testing_dataset, method='REML', select=T)
m_cof <- gam(co.flux ~ te(Tr, g_tCO, by=treatment, k=4) + treatment,
             data=testing_dataset, method='REML', select=T) # Dan wants me to drop H2Oi
summary(m_cof)
AIC(m_cof)
# p = visreg(m_cof,"COi", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
#            partial = F, rug = F)
# p2 = visreg(m_cof, xvar = 'COi', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
# plt = ggplot(p$fit, aes(COi, visregFit, linetype=treatment, fill=treatment, colour=treatment))
# plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
# plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
# plt = plt + geom_line()
# #plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
# #                     size = 5, hjust = 0, vjust = 1, family='serif')
# plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
# plt = plt + labs(x = expression(paste('c'['CO,i']," [",nmol~mol^{-1},"]")),
#                  y = expression(paste('Part. contribution to  ','f'['CO']," [",nmol~m^{-2}~s^{-1},"]")),
#                  colour='Treatment', fill='Treatment', linetype='Treatment')
# plt = plt + theme_bw()
# #plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
# plt = plt + theme(legend.position = c(0.18, 0.80), text=element_text(family="serif"),
#                   plot.title = element_text(hjust = 0.5))
# plt = plt + ggtitle( expression(paste('c'['CO,i'],' contribution to f'['CO'])))
# plt = plt + coord_cartesian(ylim= c(-0.5,3.5))
# plt_cof_coi <- plt
# plt_cof_coi
# ggsave(paste0(graphs_path, 'model3a_output_COi.jpg'), width=8, height=8, units = "cm", dpi = 1200)
# ggsave(paste0(graphs_path, 'model3a_output_COi.pdf'), width=8, height=8, units = "cm", dpi = 1200)


library(gratia)
library(metR)

m_cof <- gam(co.flux ~ te(Tr, g_tCO, by=treatment, k=4),
                 data=testing_dataset, method='REML', select=T)
summary(m_cof)
AIC(m_cof)

# m_cof_dro <- gam(co.flux ~ te(H2Oi, Tr, k=4),
#                  data=testing_dataset[which((testing_dataset$treatment == 'dro')),], method='REML', select=T)
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
# plt_cof_dro <- plt_cof_dro + ggtitle(expression(paste("Dro: ",'c'['H2O,i'],' & Tr contrib. to f'['CO'])))
# plt_cof_dro <- plt_cof_dro + labs(x=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
#                                   y=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))
plt_cof_dro <- plt_cof_dro + ggtitle(expression(paste("Droughted: ",'g'['t,CO'],' & Tr contrib. to f'['CO'])))
plt_cof_dro <- plt_cof_dro + labs(y=expression(paste(g[paste("t,CO")]," [mol ",mol^{-2}," ",s^{-1},"]")),
                                  x=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))
plt_cof_dro
ggsave(paste0(graphs_path, 'model3b_output_interaction_H2Oi_Tr.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model3b_output_interaction_H2Oi_Tr.pdf'), width=8, height=8, units = "cm", dpi = 1200)


# m_cof_irr <- gam(co.flux ~ te(H2Oi, Tr, k=4),
#                  data=testing_dataset[which((testing_dataset$treatment == 'irr')),], method='REML', select=T)
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
# plt_cof_irr <- plt_cof_irr + ggtitle(expression(paste("Irr: ",'c'['H2O,i'],' & Tr contrib. to f'['CO'])))
# plt_cof_dro <- plt_cof_dro + labs(x=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
#                                   y=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))
plt_cof_irr <- plt_cof_irr + ggtitle(expression(paste("Irrigated: ",'g'['t,CO'],' & Tr contrib. to f'['CO'])))
plt_cof_irr <- plt_cof_irr + labs(y=expression(paste(g[paste("t,CO")]," [mol ",mol^{-2}," ",s^{-1},"]")),
                                  x=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")))

plt_cof_irr
ggsave(paste0(graphs_path, 'model3c_output_interaction_H2Oi_Tr.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model3c_output_interaction_H2Oi_Tr.pdf'), width=8, height=8, units = "cm", dpi = 1200)

# library(grImport2)
# library(rsvg)
# rsvg::rsvg_svg("../graphs/model1_outline.svg", "../graphs/model1_outline2.svg")
# raw <- grImport2::readPicture("C:/Users/Jonathan/Documents/_research_statistics/2.6 - CO branch chambers/graphs/model1_outline2.svg")
# grob <- grImport2::pictureGrob(raw, just = c("left",'bottom'), x=unit(0, "cm"), y = unit(0, "cm"))

library(ggpubr)
library(cowplot)
# ggarrange(grob, plt_cof_coi,
#           plt_coi_par, plt_cof_irr,
#           plt_coi_tl, plt_cof_dro, 
#           labels = c("(1)",  "(3a)",
#                      "(2a)", "(3b)",
#                      "(2b)", "(3c)"),
#           ncol = 2, nrow = 3)
# ggsave(paste0(graphs_path, 'model_output_all.jpg'), width = 20, height = 27, units = "cm", dpi = 1200)
# ggsave(paste0(graphs_path, 'model_output_all.pdf'), width = 20, height = 27, units = "cm", dpi = 1200)

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


#  Sup: COi vs VPD ####
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
#plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
#                     size = 5, hjust = 0, vjust = 1, family='serif')
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
#plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
#                     size = 5, hjust = 0, vjust = 1, family='serif')
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("VPD [Pa]")),
                 y = expression(paste('Part. contribution to  ',c['CO,i']," [",nmol~mol^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
#plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
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

# Sup: COi vs PAR & TL####

plot_width      <- 8 # cm
plot_height     <- 8 # cm
plot_resolution <- 1200 # dpi

# Droughted, PAR and TL vs COi
m_coi_dro <- gam(COi ~ s(TL, k=4) + s(PAR, k=4), # Production GOOD
             data=testing_dataset[which(testing_dataset$treatment == 'dro'),], method='REML', select=T)
jpeg(paste0(graphs_path, 'supplement/model_coi_dro_PAR_TL.jpg'), 
     width = plot_width, height = plot_height, units='cm', res = plot_resolution)
par(mar = c(2.5,2.7,2,1),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(m_coi_dro, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(a) Droughted: PAR & ",'T'['L'],' effect on c'['CO,i'])),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()
pdf(paste0(graphs_path, 'supplement/model_coi_dro_PAR_TL.pdf'), 
    width = plot_width/2.54, height = plot_height/2.54)
par(mar = c(2.5,2.7,2,1),
    mgp = c(3, 1, 0),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(m_coi_dro, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(a) Droughted: PAR & ",'T'['L'],' effect on c'['CO,i'])),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()

# Irrigated, PAR and TL vs COi
m_coi_irr <- gam(COi ~ s(TL, k=4) + s(PAR, k=4), # Production GOOD
                 data=testing_dataset[which(testing_dataset$treatment == 'irr'),], method='REML', select=T)
jpeg(paste0(graphs_path, 'supplement/model_coi_irr_PAR_TL.jpg'), 
     width = plot_width, height = plot_height, units='cm', res = plot_resolution)
par(mar = c(2.5,2.7,2,1),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(m_coi_irr, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(b) Irrigated: PAR & ",'T'['L'],' effect on c'['CO,i'])),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()
pdf(paste0(graphs_path, 'supplement/model_coi_irr_PAR_TL.pdf'), 
    width = plot_width/2.54, height = plot_height/2.54)
par(mar = c(2.5,2.7,2,1),
    mgp = c(3, 1, 0),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(m_coi_irr, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(b) Irrigated: PAR & ",'T'['L'],' effect on c'['CO,i'])),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()

# Combined
pdf(paste0(graphs_path, 'supplement/model_coi_all.pdf'), 
    width = plot_width*2/2.54, height = plot_height/2.54)
# Layout
layout(matrix(c(1, 2), nrow = 1,
              ncol = 2, byrow = TRUE))
par(mar = c(2.5,2.5,1,1))
vis.gam(m_coi_dro, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(a) Droughted")),
        family='serif', mgp=c(1.5,0.4,0))
vis.gam(m_coi_irr, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(b) Irrigated")),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()
# JPG
jpeg(paste0(graphs_path, 'supplement/model_coi_all.jpg'), 
     width = plot_width*2, height = plot_height, units='cm', res = plot_resolution)
layout(matrix(c(1, 2), nrow = 1,
              ncol = 2, byrow = TRUE))
par(mar = c(2.5,2.5,1,1))
vis.gam(m_coi_dro, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(a) Droughted")),
        family='serif', mgp=c(1.5,0.4,0))
vis.gam(m_coi_irr, view = c('PAR','TL'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste("PAR [",mu,"mol ",m^{-2},s^{-1},"]")),
        ylab=expression(paste(T[L]," [°C]")),
        main = expression(paste("(b) Irrigated")),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()

# Sup: COi model comparison ####
# - - - - - - - - - - - - - -
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
# - - - - - - - - - - - - - -

#m3 <- gam(co.flux ~ te(H2Oi, Tr, by=treatment, k=4) + treatment, # Chosen model
#           data=testing_dataset, method='REML', select=T)
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

summary(m8)
concurvity(m4, full=F)$worst

conc_table <- concurvity(m4, full=F)$worst
conc_table[conc_table < 0.8] <- NA
conc_table[conc_table == 1] <- NA
conc_table

conc_table <- concurvity(m7, full=F)$worst
conc_table[conc_table < 0.8] <- NA
conc_table[conc_table == 1] <- NA
conc_table

# THE END ####







testing_dataset$COi_cs <- (testing_dataset$COi - mean(testing_dataset$COi, na.rm=T)) / sd(testing_dataset$COi, na.rm=T)
testing_dataset$gt_cs <- (testing_dataset$g_tCO - mean(testing_dataset$g_tCO, na.rm=T)) / sd(testing_dataset$g_tCO, na.rm=T)
m3 <- gam(co.flux ~ te(Tr, g_tCO, by=treatment, k=4),
          data=testing_dataset, method='REML', select=T)
AIC(m3)
summary(m3)
concurvity(m3, full=F)$worst




pdf(paste0(graphs_path, 'model_output_all.pdf'), 
    width = plot_width*2/2.54, height = plot_height*3/2.54)
# Layout
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3,
              ncol = 2, byrow = TRUE))
# interaction
plt_outline
print(plt_coi_par)
print(plt_coi_tl)
print(plt_cof_coi)
par(mar = c(1,1,1,1))

dev.off()
par(mfrow=c(1,1)) # undo layout


install.packages("gridExtra")
install.packages("cowplot")
install.packages("svglite")
install.packages("rsvg")
library(gridExtra)
library(cowplot)
library(svglite)
library(rsvg)
plot_cof_h2oi_tr_dro <- function() {
  par(mar = c(2.5,2.7,2,1),
      mgp = c(3, 1, 0),
      cex.axis = 0.8,   # Font size, axis labels
      cex.lab = 1.0,    # Font size, axis titles
      cex.main = 1.0)   # Font size, main title
  vis.gam(m_cof_dro, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
          xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
          ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
          main = expression(paste("(3b) Dro: ",'c'['H2O,i'],' & Tr contrib. to f'['CO'])),
          family='serif', mgp=c(1.5,0.4,0))
}
plot_cof_h2oi_tr_irr <- function() {
  library(ggplotify)
  m_cof_irr <- gam(co.flux ~ te(H2Oi, Tr, k=4) + s(COi, k=4), # Transport
                   data=testing_dataset[which((testing_dataset$treatment == 'irr')),], method='REML', select=T)
  par(mar = c(2.5,2.7,2,1),
      mgp = c(3, 1, 0),
      cex.axis = 0.8,   # Font size, axis labels
      cex.lab = 1.0,    # Font size, axis titles
      cex.main = 1.0)   # Font size, main title
  p <- vis.gam(m_cof_irr, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
          xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
          ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
          main = expression(paste("(3c) Irr: cH2O,i & Tr contrib. to CO flux")),
          family='serif', mgp=c(1.5,0.4,0))
  vis_plot_gg <- as.ggplot(~p)
  return(vis_plot_gg)
}
plot_cof_h2oi_tr_irr()

install.packages("magick")
install.packages("pdftools")
library(cowplot)
library(magick)
vis_image <- image_read_pdf(paste0(graphs_path, 'model3c_output_interaction_H2Oi_Tr.pdf'))
vis_plot_gg <- ggdraw() + draw_image(vis_image)
vis_plot_gg

library(patchwork)

# Read SVG file and convert to grob
svg_plot <- function(file_path) {
  rsvg::rsvg(file_path) %>%
    grid::rasterGrob(interpolate = TRUE)
}

svg_file <- "../graphs/model1_outline.svg"
plt_outline <- svg_plot(svg_file)
plt_cof_h2oi_tr_dro <- grid::grid.grabExpr(plot_cof_h2oi_tr_dro())
plt_cof_h2oi_tr_irr <- grid::grid.grabExpr(plot_cof_h2oi_tr_dro()) #plot_cof_h2oi_tr_irr()
print(plt_coi_par)
print(plt_coi_tl)
print(plt_cof_coi)


combined_plot <- gridExtra::grid.arrange(
  plt_outline, plt_cof_coi,
  plt_coi_par, plt_cof_h2oi_tr_dro,
  plt_coi_tl, plt_cof_h2oi_tr_irr, 
  ncol = 2, nrow = 3
)
ggsave(paste0(graphs_path, 'model_output_all.pdf'), combined_plot,
       width = 20, height = 27, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model_output_all.jpg'), combined_plot,
       width = 20, height = 27, units = "cm", dpi = 1200)

# PDF
pdf(paste0(graphs_path, 'model_output_all.pdf'), 
    width = plot_width*2/2.54, height = plot_height*3/2.54)
# Layout
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3,
              ncol = 2, byrow = F))
# interaction
print(plt_coi_par)
print(plt_coi_tl)
print(plt_cof_coi)

par(mar = c(2.5,2.7,2,1),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
m_cof_irr <- gam(co.flux ~ te(H2Oi, Tr, k=4) + s(COi, k=4), # Transport
                 data=testing_dataset[which((testing_dataset$treatment == 'irr')),], method='REML', select=T)
vis.gam(m_cof_irr, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(3c) Irr: cH2O,i & Tr contrib. to CO flux")),
        family='serif', mgp=c(1.5,0.4,0))

par(mar = c(2.5,2.7,2,1),
    mgp = c(3, 1, 0),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(m_cof_irr, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(3c) Irr: cH2O,i & Tr contrib. to CO flux")),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()












# TEST start
testing_dataset2 <- testing_dataset[which(testing_dataset$g_tCO > 0.01),]
#testing_dataset2 <- testing_dataset
#testing_dataset2 <- testing_dataset[which((testing_dataset$TL > 25) & (testing_dataset$PAR > 500)),]

mX <- gam(co.flux ~ ti(H2Oi, g_tCO, k=4,by=treatment) + s(H2Oi, k=4,by=treatment) + s(g_tCO, k=4,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ te(H2Oi, g_tCO, k=4,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ te(H2Oi, Tr, k=4,by=treatment) + treatment,
          data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ s(H2Oi, k=4,by=treatment) + s(Tr, k=4,by=treatment) + treatment,
          data=testing_dataset, method='REML', select=T)
AIC(mX)
summary(mX)
vis.gam(mX, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0),
        cond=list(treatment='dro'))
vis.gam(mX, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0),
        cond=list(treatment='irr'))
plot(mX)

# Separate treatments: IRR
mX <- gam(co.flux ~ s(H2Oi, k=4) + s(Tr, k=4) + s(COi, k=4), # Transport
          data=testing_dataset[which(testing_dataset$treatment == 'irr'),], method='REML', select=T)
mX <- gam(co.flux ~ te(H2Oi, Tr, k=4) + s(COi, k=4), # Transport
          data=testing_dataset[which(testing_dataset$treatment == 'irr'),], method='REML', select=T)
cur = testing_dataset[which((testing_dataset$treatment == 'irr') & (testing_dataset$TL > 10)),]
cur = testing_dataset[which((testing_dataset$treatment == 'irr')),]
mX <- gam(COi ~ s(TL, k=4) + s(PAR, k=4), # Production
          data=testing_dataset[which((testing_dataset$treatment == 'irr') & (testing_dataset$TL > 10)),], method='REML', select=T)
AIC(mX)
summary(mX)
plot(mX)
gam.check(mX)
vis.gam(mX, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0))
# Separate treatments: DRO
mX <- gam(co.flux ~ s(H2Oi, k=4) + s(Tr, k=4) + s(COi, k=4), # Transport
          data=testing_dataset[which(testing_dataset$treatment == 'dro'),], method='REML', select=T)
mX <- gam(co.flux ~ te(H2Oi, Tr, k=4) + s(COi, k=4), # Transport
          data=testing_dataset[which(testing_dataset$treatment == 'dro'),], method='REML', select=T)
mX <- gam(COi ~ s(TL, k=4) + s(PAR, k=4), # Production
          data=testing_dataset[which(testing_dataset$treatment == 'dro'),], method='REML', select=T)

AIC(mX)
summary(mX)
plot(mX)
vis.gam(mX, view = c('H2Oi','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(c[paste(H[2],"O,i")]," [mmol ",mol^{-1},"]")),
        ylab=expression(paste(Tr," [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0))


mX <- gam(H2Oi ~ s(PAR, k=4,by=treatment) + s(TL, k=4, by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
plot(mX)

testing_dataset2$doy <- as.numeric(strftime(testing_dataset2$timestamp, '%j'))
mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + s(TL, by=treatment) + te(g_tCO, COi, by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + te(g_tCO, TL,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + te(VPD, TL,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR, k=4,by=treatment) + s(TL, k=4,by=treatment) + s(VPD, k=4,by=treatment) + s(COi, k=4,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
mX <- gam(co.flux ~ te(PAR, TL, k=4,by=treatment) + s(H2Oi, k=4,by=treatment) + treatment,
          data=testing_dataset2, method='REML', select=T)
AIC(mX)
summary(mX)
concurvity(mX, full=F)$worst
plot(mX)
vis.gam(mX, view = c('g_tCO','TL'), plot.type = 'contour', too.far = 0.10,
        #xlab=expression(paste(T[L]," [",degree,C,"]")),
        #ylab=expression(paste("Tr [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0),
        cond=list(treatment='dro'))



mX <- gam(co.flux ~ s(PAR, k=4) + s(TL, k=4) + s(H2Oi, k=4) + s(g_tCO, k=4),
          data=testing_dataset[(testing_dataset$treatment == 'irr'),], method='REML', select=T)
AIC(mX)
summary(mX)
plot(mX)
mX <- gam(co.flux ~ s(PAR, k=4, by=treatment) + s(TL, k=4, by=treatment) + s(H2Oi, k=4, by=treatment) + s(COi, k=4, by=treatment) + treatment,
          data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR, k=4, by=treatment) + s(TL, k=4, by=treatment) + te(COi, g_tCO, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ s(VPD, by=treatment, k=4) + s(H2Oi, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
AIC(mX)
summary(mX)
plot(mX)
conc <- concurvity(mX, full=F)$worst
conc[conc < 0.8] <- NA
conc

plt = ggplot(testing_dataset)
plt = plt + geom_point(aes(x=PAR, y=co.flux, colour=treatment), size=1)
plt = plt + theme_bw()
plt
fit <- lm(formula = co.flux ~ H2Oi, data = testing_dataset[which((testing_dataset$treatment == 'dro')),])
summary(fit)

mX <- gam(H2Oi ~ s(VPD, by=treatment, k=4) + s(PAR, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T)
AIC(mX)
summary(mX)
plot(mX)

mX <- gam(co.flux ~ s(VPD, by=treatment, k=4) + s(PAR, by=treatment, k=4) + s(H2Oi, by=treatment, k=4) + s(g_tCO, by=treatment, k=4) +
            ti(g_tCO, H2Oi, by=treatment, k=4) + ti(VPD, H2Oi, by=treatment, k=4) + ti(PAR, H2Oi, by=treatment, k=4) +
            treatment,
          data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ ti(g_tCO, H2Oi, by=treatment, k=4) + s(H2Oi, by=treatment, k=4) + s(g_tCO, by=treatment, k=4) + treatment,
          data=testing_dataset, method='REML', select=T) # Transport mechanism
mX <- gam(COi ~ s(PAR, k=4) + s(TL, k=4) + s(swc, k=4),
          data=testing_dataset[which(testing_dataset$treatment == 'dro'),], method='REML', select=T) # Production of CO
AIC(mX)
summary(mX)
plot(mX)
# TEST end

# mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + s(TL, by=treatment, k=4) + treatment,
#           data=testing_dataset, method='REML', select=T)
# AIC(mX)
# summary(mX)

p = visreg(mX,"PAR", by="treatment", data = testing_dataset,  method = "REML", overlay=T, plot = F,
           partial = F, rug = F)
p2 = visreg(mX, xvar = 'PAR', by='treatment', data = testing_dataset, gg=T, method = "REML", overlay=T)
plt = ggplot(p$fit, aes(PAR, visregFit, linetype=treatment, fill=treatment, colour=treatment))
plt = plt + geom_point(data=p2$data, aes(x=x, y=y, colour=treatment), size=0.75)
plt = plt + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.5, colour=NA)
plt = plt + geom_line()
#plt = plt + annotate("text", x = min(p2$data$x), y = max(p2$data$y) + 1, label = "(a)",
#                     size = 5, hjust = 0, vjust = 1, family='serif')
plt = plt + scale_colour_manual(values=cbPalette)  + scale_fill_manual(values=cbPalette)
plt = plt + labs(x = expression(paste("PAR [",mu,mol~m^{-2}~s^{-1},"]")),
                 y = expression(paste("CO flux [",nmol~m^{-2}~s^{-1},"]")),
                 colour='Treatment', fill='Treatment', linetype='Treatment')
plt = plt + theme_bw()
#plt = plt + theme(legend.position="bottom", text=element_text(family="serif"))
plt = plt + theme(legend.position = c(0.18, 0.82), text=element_text(family="serif"))
plt
ggsave(paste0(graphs_path, 'model_output_PAR.jpg'), width=8, height=8, units = "cm", dpi = 1200)
ggsave(paste0(graphs_path, 'model_output_PAR.pdf'), width=8, height=8, units = "cm", dpi = 1200)

plot_width      <- 8 # cm
plot_height     <- 8 # cm
plot_resolution <- 1200 # dpi
jpeg(paste0(graphs_path, 'model_output_interaction_TL_Tr.jpg'), 
     width = plot_width, height = plot_height, units='cm', res = plot_resolution)
par(mar = c(2.5,2.7,2,1),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(mX, view = c('TL','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(T[L]," [",degree,C,"]")),
        ylab=expression(paste("Tr [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()
# PDF
pdf(paste0(graphs_path, 'model_output_interaction_TL_Tr.pdf'), 
    width = plot_width/2.54, height = plot_height/2.54)
par(mar = c(2.5,2.7,2,1),
    mgp = c(3, 1, 0),
    cex.axis = 0.8,   # Font size, axis labels
    cex.lab = 1.0,    # Font size, axis titles
    cex.main = 1.0)   # Font size, main title
vis.gam(mX, view = c('TL','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(T[L]," [",degree,C,"]")),
        ylab=expression(paste("Tr [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0))
dev.off()

# 4) Test multiple options ####
# - - - - - - - - - - - - -

create_combination_formulas <- function(y_var, variable_list, by=NULL, k=NULL, remove_repeated=TRUE){
  # Generate all possible combinations of variables, with max. 2 interacting variables
  combinations <- lapply(1:2, function(n) combn(variable_list, n, simplify = FALSE))
  
  # Initialize an empty list to store formulas
  formula_list <- list()
  
  if(is.null(by)){
    by_string = ''
  }else{
    by_string = paste0(', by=',by)
  }
  if(is.null(k)){
    k_string = ''
  }else{
    k_string = paste0(', k=',k)
  }
  
  # Loop through combinations of variables to create existing combinations of single variables
  for (i in seq_along(combinations)) {
    n_combos = 0
    for (j in seq_along(combinations[[i]])) {
      combo <- paste(combinations[[i]][[j]], collapse = ",")
      if(grepl(",", combo)){
        combo = paste0('te(', combo, by_string, k_string, ')')
        n_combos = n_combos + 1
      }else{
        combo = paste0('s(', combo, by_string, k_string, ')')
      }
      #if(n_combos > 1){
      #  next
      #}
      formula_list[[length(formula_list) + 1]] <- combo
    }
  }
  
  # Generate all possible combinations of variables, with max. 3 sets of variables/interactions
  combinations_final <- unique(do.call(c, lapply(1:length(variable_list), function(n) combn(unlist(formula_list), n, simplify = FALSE))))
  #combinations_final <- do.call(c, lapply(1:length(variable_list), function(n) combn(unlist(formula_list), n, simplify = FALSE)))
  
  # Initialize an empty list to store formulas
  formula_list_final <- list()
  
  # Create the formulas as text
  for (i in seq_along(combinations_final)) {
    #print(combinations_final[[i]])
    combo <- paste(combinations_final[[i]], collapse = " + ")
    if(!is.null(by)){
      combo <- paste0(combo, ' + ', by)
    }
    formula_list_final[[length(formula_list_final) + 1]] <- paste(y_var, '~', combo)
  }
  
  # Now check all the formulas:
  # - If there are repetitions, either remove them or change "te" to "ti":
  check_repetitions <- function(input_string) {
    # Clean formulas: remove "s(" and "te(" from the string
    cleaned_string <- gsub("\\b(s|te)\\(|\\)", "", input_string)
    cleaned_string <- gsub(by_string, "", cleaned_string)
    cleaned_string <- gsub(k_string, "", cleaned_string)
    # Extract variables separated by " + ", "," or " ~ "
    substrings <- unlist(strsplit(cleaned_string, "\\s*\\+\\s*|,| ~ ", perl=TRUE))
    # Check if any substring is repeated
    any(duplicated(substrings))
  }
  # Remove formulas with repetitions
  remove_formulas_with_repetitions <- function(strings) {
    # Apply the check_repeated_substring function to each string in the vector
    # Keep only the strings that do not have repeated substrings
    filtered_strings <- strings[!sapply(strings, check_repetitions)]
    
    # Return the filtered vector of strings
    return(filtered_strings)
  }
  # Adjust formulas with repetitions from "te" to "ti"
  adjust_formulas_with_repetitions <- function(input_string) {
    # Check for repeated substrings
    if (check_repetitions(input_string)) {
      # Replace "te(" with "ti("
      replaced_string <- gsub("te\\(", "ti(", input_string)
      return(replaced_string)
    } else {
      return(input_string)
    }
  }
  
  if(remove_repeated){
    formula_list_out <- remove_formulas_with_repetitions(unlist(formula_list_final))
  }else{
    formula_list_out <- adjust_formulas_with_repetitions(unlist(formula_list_final))
  }
  
  return(formula_list_out)
}
variables <- c('PAR', 'TL', 'Tr', 'g_tCO')
formulas1 <- create_combination_formulas('co.flux', variables, remove_repeated = T, by='treatment')

variables <- c('PAR', 'VPD', 'g_tCO')
formulas2 <- create_combination_formulas('co.flux', variables, remove_repeated = T, by='treatment')
all_formulas = c(formulas1, formulas2)
all_formulas

# Run all the model formulas
models <- list()
result_df <- data.frame()
i = 1
for(formula in all_formulas){
  models[[i]] <- gam(as.formula(formula), data=testing_dataset, method='REML', select=T)
  cat(paste(i, '\t', AIC(models[[i]]), '\t', formula, '\n'))
  current_result <- c(i, AIC(models[[i]]), summary(models[[i]])$dev.expl, deparse(formula))
  result_df <- rbind(result_df, current_result)
  i=i+1
}
colnames(result_df) <- c('i','AIC','R2','formula')
result_df$AIC <- as.numeric(result_df$AIC)
result_df$R2 <- as.numeric(result_df$R2)
result_df <- result_df[order(result_df$AIC), ]

# Remove unlogical interaction terms. PAR shouldn't interact with anything directly
results <- result_df[which(!grepl("te\\(PAR,g_tCO", result_df$formula)),]
#results <- result_df[which(!grepl("te\\(PAR", result_df$formula)),]
#results <- result_df[which( (!grepl("te\\(PAR", result_df$formula) | (grepl("te\\(PAR,g_tCO", result_df$formula)) )),]
View(results)


# 4b) Show test model output for some best options ####
# - - - - - - - - - - - - - - - - - - - - - - - -

mX <- gam(co.flux ~ s(Tr, by=treatment) + te(TL,PAR, by=treatment) + treatment, data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ te(VPD,PAR, by=treatment) + treatment, data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR, by=treatment) + te(TL,Tr, by=treatment) + treatment, data=testing_dataset, method='REML', select=T)
mX <- gam(co.flux ~ s(PAR) + te(TL,Tr), data=testing_dataset, method='REML', select=T)
AIC(mX)
summary(mX)
plot(mX)
#concurvity(mX, full=F)$worst

vis.gam(mX, view = c('TL','Tr'), plot.type = 'contour', too.far = 0.1,
        xlab=expression(paste(T['L']," [°C]")),
        ylab=expression(paste("PAR [",mu,mol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0))

plot(mX)

# 5) Test lag ####
# - - - - - - - - 
library(dplyr)
library(tidyr)

complete_daily_timestamps <- function(temp_df){
  # Create a dataframe with all possible combinations of timestamps and treatments
  template_df <- expand.grid(
    date = seq.Date(as.Date(min(temp_df$timestamp)), as.Date(max(temp_df$timestamp)), by = "day"),
    treatment = unique(temp_df$treatment)
  )
  
  temp_df$date <- as.Date(temp_df$timestamp)
  # Left join the original dataframe with the template dataframe
  complete_df <- left_join(template_df, temp_df, by = c("date", "treatment"))
  
  # Sort the complete dataframe by timestamp and treatment
  complete_df <- complete_df %>%
    arrange(date, treatment)
  
  complete_df$doy <- as.numeric(strftime(complete_df$date, format = "%j"))
  
  return(complete_df)
}

# Create complete timestamps for irr & dro
test1 <- complete_daily_timestamps(testing_dataset[which((testing_dataset$prob != 'b') & (testing_dataset$treatment == 'irr')),])
test2 <- complete_daily_timestamps(testing_dataset[which((testing_dataset$prob != 'b') & (testing_dataset$treatment == 'dro')),])

# Create lag
test1$lagged_Tr <- lag(test1$Tr, 40)
test2$lagged_Tr <- lag(test2$Tr, 40)

# Combine
testing = rbind(test1, test2)

# Show resulting df
plt <- ggplot(testing)
plt <- plt + geom_point(aes(x=lagged_Tr, y=TL, colour=treatment))
plt

# Test model
testing$treatment <- as.factor(testing$treatment)
mX <- gam(co.flux ~ s(PAR, by=treatment, k=4) + te(TL, lagged_Tr, by=treatment, k=4) + treatment, data=testing, method='REML', select=T)
AIC(mX)
summary(mX)
vis.gam(mX, view = c('TL','lagged_Tr'), plot.type = 'contour', too.far = 0.1)
plot(mX)

plt <- ggpairs(testing %>% select(c(co.flux, lagged_Tr, PAR, TL, VPD, g_tCO, doy)),
               columnLabels = c('CO~flux', 'Tr', 'PAR', 'T[L]', 'VPD', 'g_tCO', 'doy'), labeller='label_parsed',
               lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                            combo = wrap("dot", alpha = 0.4,            size=0.2) ),
               ggplot2::aes(colour = testing$treatment),
               upper = list(continuous = wrap("cor", size = 2.5))) +
  #labs(subtitle = "Numeric variable exploration") +
  theme_bw() + 
  theme(text=element_text(family="serif"), axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))
plt = plt + scale_colour_manual(values=cbPalette) + scale_fill_manual(values=cbPalette)
plt

# 6) Final model diagnostics plots ####
# - - - - - - - - - - - - - - - - - - -

# Run model
m_final <- gam(co.flux ~ s(PAR, by=treatment, k=4) + te(TL, Tr, by=treatment) + treatment,
               data=testing_dataset, method='REML', select=T)
# m_final <- gam(co.flux ~ s(PAR, by=treatment) + s(g_tCO, by=treatment) + te(TL,Tr, by=treatment) + treatment,
#                data=testing_dataset, method='REML', select=T)
summary(m_final)
gam.check(m_final)
AIC(m_final)

res <- resid(m_final)
qqnorm(res)
qqline(res) 

# Q-Q plot
plt = ggplot(data.frame(Standardized_Residuals = residuals(m_final)),
                        aes(sample=Standardized_Residuals))
plt = plt + stat_qq(alpha = 0.25, size=0.25)
plt = plt + geom_qq_line(color = "#808080", linetype = "dashed")
plt = plt + labs(title = "Q-Q Plot", x = "Theoretical Quantiles", y = "Standardized Residuals")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
qq_plot <- plt
qq_plot

# Histogram
plt = ggplot(data.frame(Standardized_Residuals = residuals(m_final)),
             aes(x=Standardized_Residuals))
plt = plt + geom_histogram(binwidth = 0.5, color = "black", fill = "#808080", linewidth=0.25)
plt = plt + labs(title = "Histogram of Standardized Residuals",
       x = "Standardized Residuals",
       y = "Frequency")
plt = plt + theme_bw()
plt = plt + theme(plot.title = element_text(hjust = 0.5),
                  text=element_text(size=10, family='serif'))
hist_plot <- plt

# Response
plot_df <- data.frame(Fitted_Values = fitted(m_final), Observed_Response = testing_dataset$co.flux)
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


# END ----------------------------------

# Combine both treatments
library(patchwork)
vis.gam(mX, view = c('TL','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(T[L]," [",degree,C,"]")),
        ylab=expression(paste("Tr [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0),
        cond=list(treatment='irr'))

p2 <- wrap_elements(panel = ~vis.gam(mX, view = c('TL','Tr'), plot.type = 'contour', too.far = 0.10,
        xlab=expression(paste(T[L]," [",degree,C,"]")),
        ylab=expression(paste("Tr [",mmol~m^{-2}~s^{-1},"]")),
        main = expression(paste("(b) Part. effect on CO flux [",nmol~m^{-2}~s^{-1},"]")),
        family='serif', mgp=c(1.5,0.4,0),
        cond=list(treatment='dro')))
ggsave("fig.png", p1 + p2)


# 5a) Output plots ####
# - - - - - - - - - - -

# Obtain all the subplots
predictor_names <- attr(terms(m_final), "term.labels")

plot_width      <- 8 # cm
plot_height     <- 8 # cm
plot_resolution <- 1200 # dpi
for(plot_nb in 1:length(predictor_names)){
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

# Inserted new here ####

m_final <- gam(co.flux ~ ti(VPD, Tr) + s(PAR, by=treatment, k=5) + s(TL), data=temp2,
               method='REML', select=T)

visreg3(m_final, xvar = "TL", by='treatment', data = temp2)
visreg3(m_final, xvar = "PAR", by='treatment', data = temp2, minrange=0.05, maxrange=0.95)
#visreg2d(mX, "PAR", "TL", by='treatment')
visreg(m_final, xvar = "PAR", by='treatment', overlay = T, data=temp2, scale="response", points=list(cex=0.5, pch=1), xlab="PAR", ylab='CO flux')
visreg(m_final, xvar = "TL", overlay = T, data=temp2, scale="response", points=list(cex=0.5, pch=1), xlab="PAR", ylab='CO flux')


# End new here ####

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
# m_final <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp2,
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

# TESTS ####
# - - - - - -

# Need to bin co flux data, then check again
temp2 <- temp2 %>%
  mutate(co_bin = cut_width(co.flux, width = 0.1, center = 0, labels=F))
temp3 <- temp2 %>%
  group_by(co_bin) %>%
  summarise(
    co.flux = median(co.flux),
    TL = median(TL),
    Tr = median(Tr),
    PAR = median(PAR)
  )
hist(temp3$co.flux)

mX <- gam(co.flux ~ ti(TL, Tr) + s(sqrt(PAR)) + s(sqrt(TL)) + s(sqrt(Tr)), data=temp3,
          method='REML', select=T)#, family=scat(link="identity"))
summary(mX)
gam.check(mX)
plot(mX, pages = 1, scheme = 1)
vis.gam(mX, view = c('TL','Tr'), plot.type = "contour", too.far = 0.05,
        family='serif', mgp=c(1.2,0.4,0))
AIC(mX)

mX <- gam(co.flux ~ s(par) + s(t.leaf) + s(h2o.flux) + ti(par, co2.flux) + s(co2.flux), data=temp, method='REML', select=T, family=scat(link="identity"))
summary(mX)
gam.check(mX)
plot(mX, pages = 1, scheme = 1)

library(MASS)
boxcox(lm(temp2$co.flux ~ 1))

# # Load required libraries
library(e1071)
# library(ggplot2)

# Function to calculate skewness and kurtosis for each variable
calculate_skew_kurt <- function(x) {
  skewness_value <- skewness(x)
  kurtosis_value <- kurtosis(x)
  return(data.frame(Skewness = skewness_value, Kurtosis = kurtosis_value))
}

# Calculate skewness and kurtosis for each variable
skew_kurt_results <- lapply(temp3[,c('co.flux','Tr','PAR','TL')], calculate_skew_kurt)

# Print the results
print(skew_kurt_results)
# 
# # Function to create histograms for each variable
# create_histogram <- function(x) {
#   ggplot(data.frame(x = x), aes(x)) +
#     geom_histogram(binwidth = 5, fill = "lightblue", color = "black") +
#     labs(title = paste("Histogram of", names(x)),
#          x = "Value",
#          y = "Frequency") +
#     theme_minimal()
# }
# 
# # Create histograms for each variable
# histograms <- lapply(df, create_histogram)
# 
# # Print the histograms
# print(histograms)



# Testing T distribution for heavy-tailed data:
# m_final <- gam(co.flux ~ ti(VPD, Tr) + s(PAR) + s(TL) + s(Tr, k=3), data=temp2,
#                method='REML', select=T, family=scat(link="identity"))
# summary(m_final)
# gam.check(m_final)