library(ape)
library(magrittr)
library(ggplot2)
library(scales)
source("phylo2ggphy.R")
source("plot_ggphy.R")

lasv.chronos <- read.tree("LASV_chronos.nwk")
lasv.chronos <- ladderize(lasv.chronos)
lasv.host <- rep("Human",length(lasv.chronos$tip.label))
lasv.host[grep("Josiah",lasv.chronos$tip.label)] <- "Lab"
lasv.host[grep("LM",lasv.chronos$tip.label)] <- "Mastomys"
lasv.host[grep("ZO",lasv.chronos$tip.label)] <- "Mastomys"
lasv.annotations <- data.frame(taxa=lasv.chronos$tip.label,Host=lasv.host)
lasv.tipdates <- strsplit(lasv.chronos$tip.label,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double
lasv.tipdates <- as.Date(paste(lasv.tipdates,"-01-01",sep=""))
lasv.ggphy <- phylo2ggphy(lasv.chronos,tip_dates=lasv.tipdates,branch_unit="year")
plot_ggphy(lasv.ggphy,tip_labels=F,tip_attribute=lasv.annotations,var_tip_labels="taxa",var_tip_colour="Host")
