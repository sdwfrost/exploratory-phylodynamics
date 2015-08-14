#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: treble2.R <tree> <outname>')
}

#suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(adephylo))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(inline))
suppressPackageStartupMessages(require(pander))

trName <- args[1]
outName <- args[2]

tr <- read.tree(trName)

dates <- strsplit(tr$tip.label,"-",fixed=TRUE) %>% lapply(.,tail,1) %>% unlist %>% as.double

# Cophenetic distance
dst <- cophenetic(tr)
dvar <- dst

code <- "
	int i = 0, j = 0, k = 0;
int num_seq = INTEGER(dim)[0];

double r_ij = 0, r = 0, w = 0, w_ij = 0, tl = 0, te = 0;

SEXP output;
PROTECT(output = allocVector(REALSXP,1));

for (i = 0; i<num_seq; i++)
{
  for (j = 0; j<num_seq;j++)
  {
  if (i != j)
  {
  for(k = 0; k<num_seq;k++)
  {
  if ((k !=i)&&(k!=j)&&(REAL(variance)[i*num_seq+j]>0)&&(REAL(dist)[i*num_seq+k] > REAL(dist)[i*num_seq+j])&&(REAL(dist)[j*num_seq+k] > REAL(dist)[i*num_seq+j]))
  {
  r_ij = 0;
  w_ij = 0;
  tl = 0;	te = 0;
  if (REAL(dates)[i] != REAL(dates)[j])
  {
  r_ij = (REAL(dist)[i*num_seq+k] - REAL(dist)[j*num_seq+k])/(REAL(dates)[i]-REAL(dates)[j]);
  w_ij =  pow((REAL(dates)[i]-REAL(dates)[j]),2)/(REAL(variance)[i*num_seq+k]); 
  }
  else
  {
  r_ij = fabs(REAL(dist)[i*num_seq+k] - REAL(dist)[j*num_seq+k])*2;///CHECK THIS AGAIN
  w_ij = 0.25/(REAL(variance)[i*num_seq+k]);
  }
  if (r_ij >0)  //FURTHER CONDITIONS NEEDED
  {
  tl = 0.5*(REAL(dates)[i]+REAL(dates)[j] - REAL(dist)[i*num_seq+j]/r_ij);
  te = 0.25*(REAL(dates)[i] + REAL(dates)[k] - REAL(dist)[i*num_seq+k]/r_ij) + 0.25*(REAL(dates)[k] + REAL(dates)[j] - REAL(dist)[k*num_seq+j]/r_ij);						
  }
  if ((r_ij >0)&&(te < tl)&&(tl < REAL(dates)[i])&&(tl < REAL(dates)[j]))
  {
  r += r_ij*w_ij;
  w += w_ij;	
  }
  }
  }
  }
  }
}
r = r/w;
REAL(output)[0]=r;
UNPROTECT(1);
return(output);
"
calcRate <- cfunction(signature(dist="matrix", dates="array", dim="integer", variance="matrix"), code)

r <- calcRate(dist=dst,dates=dates,dim=length(dates),variance=dvar)
dmrca <- max(distRoot(tr))
root.time <- max(dates)-dmrca/r
results <- data.frame(TMRCA=root.time,SubstRate=r)
write.table(results,paste(outName,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
pander(results)
