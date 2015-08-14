library(ape)

lasv.tree <- read.tree("ExaML_result.LASV")
lasv.tree <- unroot(lasv.tree)
lasv.tree.binary <- multi2di(lasv.tree)
lasv.tree.binary$edge.length <- lasv.tree.binary$edge.length+1e-6
write.tree(lasv.tree.binary,"LASV.nwk")
