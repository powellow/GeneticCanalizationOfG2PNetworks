## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------

### Libraries ----
library(data.table) #For quickly collapsing a list of data.frames
library(ggplot2) #For making plots
library(grid)
library(ggthemes)
require(gridExtra)
library(ggridges)
library(tidyr)
library(tidyverse)
library(patchwork)
library(plyr) #For summarizing multiple reps

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
            legend.position = "none",
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
fb_colour_vec=c("Auxin"="#CC6666","Strigolactone"="#E69F00","Sucrose"="#009E73","Cytokinin"="#0072B2","Integrator Signal"="#F0E442",target.trait="#000000")
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#colorVec = c("black",colorRampPalette(c("lightblue", "darkblue"))(2),colorRampPalette(c("lightgreen", "darkgreen"))(2))

### Select Scenarios ----

CENTER = TRUE #Places genetic variance, genic variance, and mean to zero in year 0
Scenarios=c(target.trait)


## ----Read Data, echo=FALSE------------------------------------------------------------------------------------------------------------

rawData = vector("list",length(nreps)*length(Scenarios))
i = 0L
for(SCENARIO in Scenarios){
  for(REP in 1:nreps){
    i = i+1L
    FILE = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("QtlData_Rep",REP,".rds"))
    tmp = data.frame(readRDS(FILE),
                     stringsAsFactors=FALSE)
    #tmp$Rep = REP
    tmp$Scenario = SCENARIO
    tmp$gen= as.numeric(tmp$gen)
    tmp$freq = as.numeric(tmp$freq)
    tmp$sel.coef = as.numeric(tmp$sel.coef)
    tmp$qtl = as.factor(tmp$qtl)
    rawData[[i]] = tmp
  }
}
rawData = rbindlist(rawData, fill = TRUE)

rawData <- rawData %>%
  transmute(Rep = rep,
   Cycle = as.numeric(gen),
   Frequency = as.numeric(freq),
   Selection.Coefficient = as.numeric(sel.coef),
   Trait = as.factor(gene.type),
   QTL = as.factor(qtl),
   Scenario = Scenario,
   Effect.Size = as.numeric(effect.size)) %>%
  #group_by(Trait) %>%
  transmute(Rep,Cycle,Frequency,Selection.Coefficient,Trait,QTL,Effect.Size,
   Min.Effect.Size = min(Effect.Size),
   Max.Effect.Size = max(Effect.Size),
   Norm.Effect.Size = ((Effect.Size - Min.Effect.Size)/(Max.Effect.Size - Min.Effect.Size)))

dir.create(paste(cwd,"output",scenario,"Plots",sep="/"))
setwd(paste(cwd,"output",scenario,"Plots",sep="/"))



## ----Estimate Summary Statistics, echo=FALSE------------------------------------------------------------------------------------------

tmp = ddply(
  rawData,
  c("Cycle","Trait"),
  summarize,
  Mean = mean(Frequency),
  SE = sd(Frequency) / sqrt(nreps),
  SD = sd(Frequency)
)



## ----Plot - Mean Across Reps, echo=FALSE----------------------------------------------------------------------------------------------

AF_Average <- ggplot(tmp) +
  geom_ribbon(aes(
    x = Cycle,
    ymin = Mean - SD,
    ymax = Mean + SD,
    fill = Trait
  ),
  alpha = 0.2,show.legend=FALSE) +
  geom_line(aes(x=Cycle,y=Mean,colour=Trait),
            size=2) +
  guides(alpha = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.25)),
  ) +
  scale_x_continuous("Cycle",limits = c(0,n.cycles)) +
  scale_y_continuous("Allele Frequency",limits = c(-0.2,1.2),breaks = seq(0,1,0.2)) +
  scale_color_manual(values = fb_colour_vec) +
  #ggtitle("Allele Frequency Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_Publication()

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.Average.pdf",sep="_"),plot = AF_Average, width = 20, height = 12, units = "in")


## ----Plot All Runs, echo=FALSE--------------------------------------------------------------------------------------------------------

AF_AR <- ggplot(rawData) +
  geom_point(aes(x=Cycle,y=Frequency,colour=Trait),
            size=2) +
  guides(alpha = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.25)),
  ) +
  scale_x_continuous("Cycle",limits = c(0, n.cycles)) +
  scale_y_continuous("Allele Frequency",limits=c(-0.2,1.2),breaks = seq(0,1,0.2)) +
  scale_color_manual(values = fb_colour_vec) +
  ggtitle("Allele Freq. Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_Publication()

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.AllReps.pdf",sep="_"),plot = AF_AR, width = 20, height = 12, units = "in")


## ----Ridged Density Plot, echo= F-----------------------------------------------------------------------------------------------------

cycles_of_interest <- c(1,5,10,30)
#cycles_of_interest <- c(1,25,50,75,100)

AF <- rawData %>% filter(Cycle %in% cycles_of_interest) %>% group_by(Trait) %>%
  mutate(QTL = fct_reorder(QTL, desc(Effect.Size)),Trait = factor(Trait,levels = c("Sucrose","Cytokinin","Strigolactone","Auxin")))

AF$QTL.Trait <- interaction(AF$QTL, AF$Trait)

AF_Distr <- ggplot(AF,aes(x=Frequency,y=QTL.Trait,fill=Trait)) +
  geom_density_ridges2(scale=0.8) +
  #geom_density_ridges2(stat = "binline", binwidth = 0.1,scale=0.95) +
  scale_x_continuous("Allele Frequency",limits = c(-0.2,1.2), breaks = seq(0,1,0.5)) +
  scale_y_discrete("Causal Genetic Loci") +
  scale_colour_manual(values = fb_colour_vec) +
  #ggtitle("Allele Freq. Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_no_ticks() +
  facet_grid(cols = vars(Cycle))

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.Distributions.pdf",sep="_"),plot = AF_Distr, width = 20, height = 12, units = "in")

AF_Hist <- ggplot(AF,aes(x=Frequency,y=QTL.Trait,fill=Trait)) +
  #geom_density_ridges2(scale=0.95) +
  geom_density_ridges2(stat = "binline", binwidth = 0.1,scale=0.95) +
  scale_x_continuous("Allele Frequency",limits = c(-0.1,1.1),breaks = seq(0,1,0.5)) +
  scale_y_discrete("Causal Genetic Loci",labels = "") +
  scale_colour_manual(values = fb_colour_vec) +
  #ggtitle("Allele Freq. Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_no_ticks() +
  facet_grid(cols = vars(Cycle))

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.Histograms.pdf",sep="_"),plot = AF_Hist, width = 20, height = 12, units = "in")

## ----AF vs Eff. Size------------------------------------------------------------------------------------------------------------------
#install.packages('ggbeeswarm')
library(ggbeeswarm)

AF_ES = ddply(
  rawData,
  c("Cycle","Trait","QTL"),
  summarize,
  Mean = mean(Frequency),
  SE = sd(Frequency) / sqrt(nreps),
  SD = sd(Frequency)
)

cycles_of_interest <- c(1,5,10,30)
#cycles_of_interest <- c(seq(1,10,1),15,20)

AF_ES <- rawData %>% filter(Cycle %in% cycles_of_interest)

AF_Seagull <- ggplot(AF_ES,aes(x=Frequency,y=Norm.Effect.Size,colour=Trait)) +
  geom_point(size=2) +
  #geom_quasirandom() +
  scale_x_continuous("Allele Frequency",limits = c(0,1), breaks = seq(0,1,0.5)) +
  scale_y_continuous("Genetic Effect Size",limits = c(-0.1,1.1), breaks = seq(0,1,0.2)) +
  scale_colour_manual(values = fb_colour_vec) +
  #ggtitle("Allele Freq. Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_Publication() +
  theme(legend.position = "none") +
  facet_wrap(. ~ Cycle,nrow = 1)

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.ES.Seagull.pdf",sep="_"),plot = AF_Seagull, width = 20, height = 12, units = "in")

## ----AF vs Selection Coeff------------------------------------------------------------------------------------------------------------------

library(ggbeeswarm)

# ES_SC = ddply(
#   rawData,
#   c("Cycle","Trait","QTL"),
#   summarize,
#   Mean = mean(Selection.Coefficent),
#   SE = sd(Selection.Coefficent) / sqrt(nreps),
#   SD = sd(Selection.Coefficent)
# )

cycles_of_interest <- c(1,5,10,30)
#cycles_of_interest <- c(seq(1,10,1),15,20)

ES_SC <- rawData %>% filter(Cycle %in% cycles_of_interest)

ES_Sel <- ggplot(ES_SC,aes(x=Selection.Coefficient,y=Norm.Effect.Size,colour=Trait)) +
  geom_point(size=2) +
  #geom_quasirandom() +
  scale_x_continuous("Selection Coefficient",limits = c(0,1), breaks = seq(0,1,0.2)) +
  scale_y_continuous("Genetic Effect Size",limits = c(-0.1,1.1), breaks = seq(0,1,0.2)) +
  scale_colour_manual(values = fb_colour_vec) +
  #ggtitle("Allele Freq. Changes of Network") +
  scale_fill_manual(values = fb_colour_vec) + theme_Publication() +
  theme(legend.position = "none") +
  facet_wrap(. ~ Cycle,nrow = 1)

#ggsave(file=paste(selection.direction,target.trait,Nutrients,paste0("H2=",H2),"AF.ES.Seagull.pdf",sep="_"),plot = AF_Seagull, width = 20, height = 12, units = "in")


## -------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)

#AF_Distr + AF_Average
#AF_Hist + AF_Average

#AF_Distr/AF_Average
#AF_Hist/AF_Average


