###################################
## Load drug targets from STITCH ##
###################################
source("code/vars.R")

saveGraphs = T

## Load the stitch networks
drugs = c("Fenofibrate", "T0901317")
names(drugs) = drugs

## Map human targets to mouse genes and convert to igraph
allsymbols = as.list(org.Mm.egSYMBOL2EG)
names(allsymbols) = toupper(names(allsymbols))
drugs.graphs = lapply(drugs, function(d) {
  d.table = read.delim(paste(dataPath, d, ".txt", sep=""), sep="\t", as.is=T)
  prots = unique(unlist(d.table[,c(1,2)]))
  prots = setdiff(toupper(prots), toupper(d))
  ids = allsymbols[prots]
  ids[toupper(d)] = d
  d.table[,1] = as.character(ids[toupper(d.table[,1])])
  d.table[,2] = as.character(ids[toupper(d.table[,2])])
  g = graph.data.frame(d.table, directed=F)
})

## Add attributes
library(org.Mm.eg.db)
drugs.graphs = lapply(drugs.graphs, function(g) {
  V(g)[!(name %in% drugs)]$name = paste("L:", V(g)[!(name %in% drugs)]$name, sep="")
  
  labels = as.character(as.list(org.Mm.egSYMBOL)[sub("L:", "", V(g)$name)])
  labels[labels == "NULL"] = sub("L:", "", V(g)$name[labels == "NULL"])
  V(g)$label = labels
  
  g
})
save(drugs, drugs.graphs, file = paste(outPath, "stitch.graphs.RData", sep=""))