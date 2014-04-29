library(igraph)

source("code/vars.R")
myOutPath = paste0(outPath, 'networks/');
netPath = paste0(dataPath, 'networks/')

cacheFile = paste(myOutPath, "graphs.RData", sep="")
if(file.exists(cacheFile)) {
  load(cacheFile)  
} else {
  ## Read network
  files = listFiles(netPath, "gml")
  graphs = lapply(files, read.graph, "gml")
  names(graphs) = files
  graph = mergeGraphs(graphs)
  
  ## Add gene symbols
  library(org.Mm.eg.db)
  library(GO.db)
  ids = sub("L:", "", V(graph)$name)
  V(graph)$entrez = ids
  V(graph)$label = as.character(as.list(org.Mm.egSYMBOL)[ids])
  
  ## Add GO annotation
  goids = lapply(as.list(org.Mm.egGO)[ids], names)
  id2term = Term(GOTERM)
  id2ont = Ontology(GOTERM)
  go.bp = lapply(goids, function(g) {
    g = setdiff(g, "GO:0008150")
    id2term[g[id2ont[g] == "BP"]]
  })
  V(graph)$goids = sapply(goids, paste, collapse=";")
  V(graph)$goBP = sapply(go.bp, paste, collapse=";")
  
  graph.full = graph
  graph = as.undirected(graph.full, "each")
  graph = simplify(graph)
  save(graph.full, graph, file=cacheFile)
}
