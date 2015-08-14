library(ape)
library(magrittr)
library(treeImbalance)

lasv.chronos <- read.tree("LASV_chronos.nwk")
lasv.tipdates <- strsplit(lasv.chronos$tip.label,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double

lasv.cherries.obs <- ct(lasv.chronos)
lasv.sackins.obs <- snt(lasv.chronos)

treelist <- list()
nsims <- 1000
for(i in 1:nsims){
  treelist[[i]] <- getSimTree(lasv.chronos)
}

lasv.cherries.sim <- lapply(treelist,ct)
lasv.sackins.sim <- lapply(treelist,snt)

par(mfrow=c(1,2),pty="s")
plot(max(lasv.tipdates)-lasv.cherries.obs[[1]],lasv.cherries.obs[[2]],col="red",xlab="Year",ylab="Cherries",las=1,ylim=c(0,30),type="n",main="Cherries")
for(i in 1:nsims){
  lines(max(lasv.tipdates)-lasv.cherries.sim[[i]][[1]],lasv.cherries.sim[[i]][[2]],type="s",col="gray")
}
lines(max(lasv.tipdates)-lasv.cherries.obs[[1]],lasv.cherries.obs[[2]],type="s",col="red",lwd=2)
plot(max(lasv.tipdates)-lasv.sackins.obs[[1]],lasv.sackins.obs[[2]],col="red",xlab="Year",ylab="Sackins",las=1,type="n",main="Sackin's index")
for(i in 1:nsims){
  lines(max(lasv.tipdates)-lasv.sackins.sim[[i]][[1]],lasv.sackins.sim[[i]][[2]],type="s",col="gray")
}
lines(max(lasv.tipdates)-lasv.sackins.obs[[1]],lasv.sackins.obs[[2]],type="s",col="red",lwd=2)
