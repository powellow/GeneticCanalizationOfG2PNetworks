rm(list=ls())

conflictRules("MASS", exclude = "select")
#install.packages('hapsim',exclude="select")
library('hapsim')
library(plyr)
library('AlphaSimR')
library('tidyverse')
library('ggpubr')
library('ggridges')
library(purrr)        # Functional programming
library(dplyr)        # Data wrangling
library(tidyr)        # Tidy-ing data
library(broom)        # List columns within tibbles
library(corrplot)
library(gridExtra)
library(data.table) #For quickly collapsing a list of data.frames
library(ggplot2)
library(ggthemes)
library(ggforce)
require('Ghat')
require('rrBLUP')
library(latex2exp)
library(patchwork)
library(scales)
library(Surrogate)
library(beepr)

options(scipen = 1000000)

# suppressMessages(require('data.table'))
# suppressMessages(require('jsonlite'))
#suppressMessages(require('asreml'))
#asreml.options(fail="soft",ai.sing=TRUE)
# #asreml.license.activate()

source('./code/Functions.R')
source('./code/network_functions.R')
source('./analysis/specify_parameters.R')


for(H2 in heritabilities){
  set.seed(231)
  for(rep_count in 1:nreps){
    ### SetUp Folder Hierarchy
    cwd <- getwd()
    dir.create(paste(cwd,"output",scenario,sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,selection.or.drift,sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,sep="/"))
    dir.create(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"))

    source('./analysis/branching_network.R')
  }
}

beep()
