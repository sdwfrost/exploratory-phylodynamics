#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: rtt.R <tree> <outname>')
}

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(adephylo))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(pander))

trName <- args[1]
outName <- args[2]

tr <- read.tree(trName)

dates <- strsplit(tr$tip.label,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double

tr.rtt <- rtt(tr,dates)
write.tree(tr.rtt,paste(outName,".nwk",sep=""))

# Now summarise results

rd <- distRoot(tr.rtt)
td <- dates[match(tr.rtt$tip.label,tr$tip.label)]
rtt.lm <- lm(rd~td)
root.time <- unname(-as.double(coef(rtt.lm)[1])/coef(rtt.lm)[2])
results <- data.frame(TMRCA=root.time,SubstRate=as.double(coef(rtt.lm)[2]))
write.table(results,paste(outName,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
pander(results)
