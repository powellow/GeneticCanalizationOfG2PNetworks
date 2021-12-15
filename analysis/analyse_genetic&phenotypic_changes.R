## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------

### Libraries ----
#library(data.table) #For quickly collapsing a list of data.frames
# library(ggplot2) #For making plots
# library(grid)
# library(ggthemes)
# require(gridExtra)
# library(ggridges)
# library(tidyr)
#library(tidyverse)
# library(patchwork)
# library(plyr) #For summarizing multiple reps

#remove.packages('plyr')


### Custom ggplot theme ----
#Paper Scaling
theme_Publication <- function(base_size = 20,base_family = "serif") {
  (theme_foundation(base_size = base_size, base_family = base_family)
   + theme(plot.title = element_text(face = "bold",size = rel(1.25),margin = margin(b = 20,t = 20,r = 20, l = 20),hjust=0.5),
           text=element_text(family="sans"),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold", size = rel(1.4)),
           axis.title.y = element_text(angle = 90),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(size = rel(1.4)),
           axis.line = element_line(colour = NA, size = rel(0.75)),
           axis.line.x = element_line(colour = 'black'),
           axis.line.y = element_line(colour = 'black'),
           axis.ticks = element_line(),
           axis.title.y.right = element_text(angle = 90,vjust = 2),
           panel.grid.major.y = element_line(colour = "#f0f0f0"),
           panel.grid.major.x = element_blank(),
           panel.grid.minor = element_blank(),
           panel.spacing = unit(3, "lines"),
           legend.key = element_rect(colour = NA),
           legend.text = element_text(size = rel(1)),
           legend.position = "none",
           legend.direction = "vertical",
           legend.key.size = unit(1, "cm"),
           #legend.margin = unit(margin(c(0, 0, 0, 0)), "mm"),
           legend.title = element_text(face = "bold", size = rel(1)),
           plot.margin = unit(c(1, 5, 5, 1), "mm"),
           strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
           strip.text = element_text(face = "bold", size = rel(1.4)),
           legend.background=element_blank(),
           #legend.spacing.y = unit(0.5, "mm")
   )
  )
}

theme_no_ticks <- function(base_size = 20,base_family = "serif") {
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(plot.title = element_text(face = "bold",size = rel(1),margin = margin(b = 20,t = 20,r = 20, l = 20),hjust=0.5),
            text = element_text(family="sans"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(0.75)),
            axis.title.y = element_text(angle = 90),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(0.75)),
            axis.text.y = element_blank(),
            axis.line = element_line(colour = NA, size = rel(0.75)),
            axis.line.x = element_line(colour = 'black'),
            axis.line.y = element_line(colour = 'black'),
            axis.ticks = element_blank(),
            panel.grid.major.y = element_line(colour = "#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size = rel(0.75)),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(1, "cm"),
            #legend.margin = unit(margin(c(0, 0, 0, 0)), "mm"),
            legend.title = element_text(face = "bold", size = rel(0.75)),
            plot.margin = unit(c(1, 2, 3, 1), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold", size = rel(0.75)),
            legend.background=element_blank(),
            #legend.spacing.y = unit(0.5, "mm")
    )
  )
}



## ----Specify Simulation Scenarios-----------------------------------------------------------------------------------------------------

### Color Palletes & Vectors ----
colorVec = c("black","blue","green","")
#fb_colour_vec=c("Auxin"="#CC6666","Strigolactones"="#E69F00","Sucrose"="#009E73","Cytokinins"="#0072B2","Integrator Signal"="#F0E442","Time to Bud Outgrowth"="#000000")
fb_colour_vec=c("Auxin"="#CC6666","Strigolactones"="#E69F00","Sucrose"="#009E73","Cytokinins"="#0072B2","Integrator Signal"="#000000")
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#colorVec = c("black",colorRampPalette(c("lightblue", "darkblue"))(2),colorRampPalette(c("lightgreen", "darkgreen"))(2))

### Select Scenarios ----

CENTER = TRUE #Places genetic variance, genic variance, and mean to zero in year 0
Scenarios=c(target.trait)


## ----Read Data, echo=FALSE------------------------------------------------------------------------------------------------------------

###  Genetic Summary Data ----

GenRawData = vector("list",length(nreps)*length(Scenarios))
i = 0L
for(SCENARIO in Scenarios){
  for(REP in 1:nreps){
    i = i+1L
    FILE = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("GeneticData_Rep",REP,".rds"))
    tmp = data.frame(readRDS(FILE),
                     stringsAsFactors=FALSE)
    tmp$Rep = tmp$rep
    tmp$`Physiological Trait` = tmp$trait
    tmp$Scenario = SCENARIO
    tmp$cycle = as.numeric(tmp$cycle)
    tmp$Mean = as.numeric(tmp$Mean)
    tmp$Variance = as.numeric(tmp$Variance)
    GenRawData[[i]] = tmp
  }
}
GenRawData = rbindlist(GenRawData, fill = TRUE)
GenRawData <- GenRawData %>% select(-c("rep","trait"))

## Phenotypic Summary Data ----

PhenRawData = vector("list",length(nreps)*length(Scenarios))
i = 0L
for(SCENARIO in Scenarios){
  for(REP in 1:nreps){
    i = i+1L
    FILE = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("PhenotypicData_Rep",REP,".rds"))
    tmp = data.frame(readRDS(FILE),
                     stringsAsFactors=FALSE)
    tmp$Rep = tmp$rep
    tmp$`Physiological Trait` = tmp$trait
    tmp$Scenario = SCENARIO
    tmp$cycle = as.numeric(tmp$cycle)
    tmp$Mean = as.numeric(tmp$Mean)
    tmp$Variance = as.numeric(tmp$Variance)
    PhenRawData[[i]] = tmp
  }
}
PhenRawData = rbindlist(PhenRawData, fill = TRUE)
PhenRawData <- PhenRawData %>% select(-c("rep","trait"))
PhenRawData$log_Mean <- log(PhenRawData$Mean);PhenRawData$log_Mean[PhenRawData$log_Mean=="NaN"] <- -11.51293
#PhenRawData$log_Mean <- PhenRawData$log_Mean - min(PhenRawData$log_Mean);
summary(PhenRawData$log_Mean)

#GenRawData = GenRawData %>% filter(.,`Physiological Trait` != target.trait)
# rawData <- PhenRawData %>% filter(`Physiological Trait` == target.trait) %>%
#   rbind(.,GenRawData)

### Setwd

dir.create(paste(cwd,"output",scenario,"Plots",sep="/"))
setwd(paste(cwd,"output",scenario,"Plots",sep="/"))

## ----Estimate Summary Statistics, echo=FALSE------------------------------------------------------------------------------------------

# Value used to transform the data
coeff <- 10

# tmp = ddply(
#   rawData,
#   c("cycle","`Physiological Trait`","Metric"),
#   summarize,
#   mean = mean(Mean),
#   pheno_mean = mean/coeff,
#   SE = sd(Mean) / sqrt(nreps),
#   SD = sd(Mean)
# )
#
# tmp <- tmp %>% group_by(`Physiological Trait`) %>%
#   mutate(MinMean = min(mean), MaxMean = max(mean), Norm_Mean = ((mean - MinMean)/(MaxMean - MinMean)))
# tmp[tmp$`Physiological Trait` == target.trait,"Norm_Mean"] <- tmp[tmp$`Physiological Trait` == target.trait,"mean"]
# tmp[tmp$Metric == "Genetic Value","pheno_mean"] <- NA
# tmp[tmp$Metric == "Phenotypic Value","mean"] <- NA

## ----Plot - Mean Across Reps, echo=FALSE----------------------------------------------------------------------------------------------

# StandAddGenRawData = GenRawData %>%
#   filter(.,grepl('Normalised', `Physiological Trait`)) %>%
#   group_by(cycle,`Physiological Trait`) %>%
#   mutate(.,overall_mean = mean(Mean))
#
# ### Replace trait names with prettier names --
# StandAddGenRawData[StandAddGenRawData$`Physiological Trait` == "Normalised.Auxin.AdditiveGeneticValue","`Physiological Trait`"] <- "Auxin"
# StandAddGenRawData[StandAddGenRawData$`Physiological Trait` == "Normalised.Cytokinin.AdditiveGeneticValue","`Physiological Trait`"] <- "Cytokinin"
# StandAddGenRawData[StandAddGenRawData$`Physiological Trait` == "Normalised.Sucrose.AdditiveGeneticValue","`Physiological Trait`"] <- "Sucrose"
# StandAddGenRawData[StandAddGenRawData$`Physiological Trait` == "Normalised.Strigolactone.AdditiveGeneticValue","`Physiological Trait`"] <- "Strigolactone"
# #
# # tmp = ddply(
# #   AddGenRawData,
# #   c("cycle","`Physiological Trait`","Metric"),
# #   summarize,
# #   mean = mean(Mean),
# #   SE = sd(Mean) / sqrt(nreps),
# #   SD = sd(Mean)
# # )
#
# StandAddPopMean_Line <- ggplot(StandAddGenRawData,aes(group=interaction(Rep,`Physiological Trait`))) +
#   geom_line(aes(x=cycle,y=Mean,colour=`Physiological Trait`),
#             size=0.1,alpha=0.1) +
#   geom_line(aes(x=cycle,y=overall_mean,colour=`Physiological Trait`),size=1.5,alpha=0.9) +
#   guides(alpha = FALSE) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.title = element_text(size = rel(1.5)),
#     axis.text = element_text(size = rel(1.25)),
#   ) +
#   scale_x_continuous("Cycle",limits = c(0, n.cycles)) +
#   scale_y_continuous("(Normalised) Additive Genetic Value",limits = c(-0.2,1.2),breaks = seq(0,1,0.2)) +
#   scale_color_manual(values = fb_colour_vec) +
#   ggtitle("Selection Trajectories of Network") +
#   #scale_fill_manual(values = fb_colour_vec) +
#   theme_Publication()

# AddGenRawData = GenRawData %>%
#   filter(.,!grepl('Normalised', `Physiological Trait`)) %>%
#   filter(.,grepl('Additive', `Physiological Trait`)) %>%
#   group_by(cycle,`Physiological Trait`) %>%
#   mutate(.,overall_mean = mean(Mean))

# ### Replace trait names with prettier names --
# AddGenRawData[AddGenRawData$`Physiological Trait` == "Auxin.AdditiveGeneticValue","`Physiological Trait`"] <- "Auxin"
# AddGenRawData[AddGenRawData$`Physiological Trait` == "Cytokinin.AdditiveGeneticValue","`Physiological Trait`"] <- "Cytokinin"
# AddGenRawData[AddGenRawData$`Physiological Trait` == "Sucrose.AdditiveGeneticValue","`Physiological Trait`"] <- "Sucrose"
# AddGenRawData[AddGenRawData$`Physiological Trait` == "Strigolactone.AdditiveGeneticValue","`Physiological Trait`"] <- "Strigolactone"
#
# AddPopMean_Line <- ggplot(AddGenRawData,aes(group=interaction(Rep,`Physiological Trait`))) +
#   geom_line(aes(x=cycle,y=Mean,colour=`Physiological Trait`),
#             size=0.1,alpha=0.1) +
#   geom_line(aes(x=cycle,y=overall_mean,colour=`Physiological Trait`),size=1.5,alpha=0.9) +
#   guides(alpha = FALSE) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.title = element_text(size = rel(1.5)),
#     axis.text = element_text(size = rel(1.25)),
#   ) +
#   scale_x_continuous("Cycle",limits = c(0, n.cycles)) +
#   scale_y_continuous("Additive Genetic Value",limits = c(-1.2,5.2),breaks = seq(-1,5,1)) +
#   scale_color_manual(values = fb_colour_vec) +
#   #ggtitle("Selection Trajectories of Network") +
#   #scale_fill_manual(values = fb_colour_vec) +
#   theme_Publication()

#GenRawData$Mean[GenRawData$Mean<=0] <- 0.0001
#GenRawData$Mean[GenRawData$`Physiological Trait`=="TimeToBudOutgrowth.TotalGeneticValue"] <- log(GenRawData$Mean[GenRawData$`Physiological Trait`=="TimeToBudOutgrowth.TotalGeneticValue"])

TotalGenRawData = GenRawData %>%
  filter(.,grepl('Normalised', `Physiological Trait`)) %>%
  filter(.,grepl('Total', `Physiological Trait`)) %>%
  #filter(.,!grepl('TimeToBudOutgrowth', `Physiological Trait`)) %>%
  # filter(.,!grepl(0, cycle)) %>%
  group_by(cycle,`Physiological Trait`) %>%
  mutate(.,overall_mean = mean(Mean),SE = sd(Mean)/sqrt(nreps),SD = sd(Mean))

TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.IntegratorSignal.TotalGeneticValue","Physiological Trait"] <- "Integrator Signal"

# TotalGenRawData = GenRawData %>% filter(.,!grepl('Total', `Physiological Trait`))  %>%
#   filter(.,!grepl('TimeToBudOutgrowth', `Physiological Trait`)) %>%
#   ddply(.,
#   c("cycle","`Physiological Trait`","Metric"),
#   summarize,
#   mean = mean(Mean),
#   SE = sd(Mean) / sqrt(nreps),
#   SD = sd(Mean)
# )


### Replace trait names with prettier names --
TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.Auxin.TotalGeneticValue","Physiological Trait"] <- "Auxin"
TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.Cytokinin.TotalGeneticValue","Physiological Trait"] <- "Cytokinins"
TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.Sucrose.TotalGeneticValue","Physiological Trait"] <- "Sucrose"
TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.Strigolactone.TotalGeneticValue","Physiological Trait"] <- "Strigolactones"
#TotalGenRawData[TotalGenRawData$`Physiological Trait` == "Normalised.TimeToBudOutgrowth.TotalGeneticValue","Physiological Trait"] <- "Time to Bud Outgrowth"

TotalGenRawData = TotalGenRawData %>% filter(.,!grepl('Total', `Physiological Trait`))

TotalPopMean_Line <- ggplot(data = TotalGenRawData) +
  geom_line(aes(x=cycle,y=Mean,group=interaction(Rep,`Physiological Trait`),colour=`Physiological Trait`),
            size=0.5,alpha=0.1) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  geom_line(aes(x=cycle,y=overall_mean,colour=`Physiological Trait`),
            size=2) +
  # geom_ribbon(aes(x = cycle,
  #                 ymin = overall_mean - SD,
  #                 ymax = overall_mean + SD,
  #                 fill = `Physiological Trait`
  # ),alpha=0.3,show.legend=FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.25)),
  ) +
  scale_x_continuous("Selection Cycle",limits = c(1, n.cycles),breaks=c(1,seq(10,n.cycles,10))) +
  scale_y_continuous("Population Genetic Mean", limits = c(-0.05,1.05),breaks=seq(0,1,0.1)) +
  scale_color_manual(values = fb_colour_vec) +
  #ggtitle("Selection Trajectories") +
  scale_fill_manual(values = fb_colour_vec) +
  theme_Publication()

#AddPopMean_Line + TotalPopMean_Line

# PopMean_LineBar <- ggplot(tmp,aes(x=cycle)) +
#   geom_bar(aes(y=pheno_mean),stat="identity",
#             size=.1,colour="black",alpha=0.4) +
#   geom_line(aes(y=mean,colour=`Physiological Trait`),
#             size=2) +
#    geom_ribbon(aes(
#     x = cycle,
#     ymin = mean - SD,
#     ymax = mean + SD,
#     fill = `Physiological Trait`
#   ),
#   alpha = 0.2,show.legend=FALSE) +
#   guides(alpha = FALSE) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.title = element_text(size = rel(1.5)),
#     axis.text = element_text(size = rel(1.25)),
#   ) +
#   scale_x_continuous("Cycle",limits = c(0, n.cycles)) +
#   scale_y_continuous("Endogenous Signal Synthesis",limits = c(-0.2,1.2),breaks = seq(0,1,0.2),
#                      sec.axis = sec_axis(~.*coeff, name=paste(target.trait,"(Days)"),breaks = seq(0,10,2))) +
#   scale_color_manual(values = fb_colour_vec) +
#   ggtitle("Trajectories of Network") +
#   #scale_fill_manual(values = fb_colour_vec) +
#   #theme_Publication()

#PopMean_Line
#PopMean_LineBar

# tmp = ddply(
#   PhenRawData,
#   c("cycle","`Physiological Trait`","Metric"),
#   summarize,
#   mean = mean(Mean),
#   pheno_mean = mean(log_Mean)/coeff,
#   SE = sd(Mean) / sqrt(nreps),
#   SD = sd(Mean)
# )

# Pheno_PopMean_Line <- ggplot(tmp) +
#   geom_ribbon(aes(
#     x = cycle,
#     ymin = mean - SD,
#     ymax = mean + SD,
#     fill = `Physiological Trait`
#   ),
#   alpha = 0.2,show.legend=FALSE) +
#   geom_line(aes(x=cycle,y=mean,colour=`Physiological Trait`),
#             size=2) +
#   guides(alpha = FALSE) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.title = element_text(size = rel(1.5)),
#     axis.text = element_text(size = rel(1.5)),
#   ) +
#   scale_x_continuous("Generations",limits = c(0, n.cycles)) +
#   scale_y_continuous("Phenotypic Mean",limits = c(-0.2,1.2),breaks = seq(0,1,0.2)) +
#   #scale_color_manual(values = fb_colour_vec) +
#   ggtitle("Selection Trajectories of Shoot Branching Network") +
#   #scale_fill_manual(values = fb_colour_vec) +
#   theme_Publication()


## ----Plot Final Figure, echo=FALSE----------------------------------------------------------------------------------------------------

ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AGVs.pdf",sep="_"),plot = TotalPopMean_Line, width = 12, height = 7, units = "in")
ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"ES_vs_AF.pdf",sep="_"),plot = AF_Seagull, width = 12, height = 7, units = "in")



#merged <- (TotalPopMean_Line/AF_Seagull) + plot_layout(heights = unit(c(10,10, 5), c('in','in','null')),guides = 'collect')

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("ConstantH2=",H2),"Figure2.pdf",sep="_"),plot = merged, width = 25, height = 25, units = "in")

#merged2 <- (TotalPopMean_Line/ES_Sel) + plot_layout(heights = unit(c(10,10, 5), c('in','in','null')),guides = 'collect')

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("ConstantH2=",H2),"Figure2_SC.pdf",sep="_"),plot = merged2, width = 25, height = 25, units = "in")

#merged3 <- (TotalPopMean_Line/AF_Hist)  + plot_layout(heights = unit(c(10,10, 5), c('in','in','null')),guides = 'collect')

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("ConstantH2=",H2),"Figure2_AFHist.pdf",sep="_"),plot = merged3, width = 25, height = 25, units = "in")

