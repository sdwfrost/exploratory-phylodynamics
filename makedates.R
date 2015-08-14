library(ape)
library(magrittr)

lasv <- read.dna("LASV.phy",format="sequential")
lasv.tipnames <- row.names(lasv)
lasv.tipdates <- lasv.tipnames %>% strsplit(.,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double
lasv.df <- data.frame(c(length(lasv.tipnames),lasv.tipnames),c("",lasv.tipdates))
write.table(lasv.df,"LASV.tipdates",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
