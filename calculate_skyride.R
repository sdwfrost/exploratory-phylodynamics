library(ape)
library(INLA)
library(magrittr)
library(ggplot2)
source("skyride.R")

lasv.chronos <- read.tree("LASV_chronos.nwk")
lasv.sr <- calculate.heterochronous.skyride(lasv.chronos)
lasv.tipdates <- strsplit(lasv.chronos$tip.label,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double
lasv.sr$year <- max(lasv.tipdates)-lasv.sr$time
ggplot(data=lasv.sr,aes(x=year))+geom_line(aes(y=sr.median))+ylab(expression(N[e]))+xlab("Time since present")+scale_y_log10(limits=c(10,1000))+theme(strip.text.x=element_text(size = rel(0.6)))
