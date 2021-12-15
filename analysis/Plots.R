rm(list=ls())

### Libraries ----
library(plyr)
library(ggthemes)
library(ggforce)
library(patchwork)
library(tidyr)
library(beepr)
library(data.table)
library(dplyr)

source('./code/network_functions.R')
source('./analysis/specify_parameters.R')

for(H2 in heritabilities){
    cwd <- getwd()
    source('./analysis/analyse_all_freq_changes.R')
    source('./analysis/analyse_genetic&phenotypic_changes.R')
}
beep()
