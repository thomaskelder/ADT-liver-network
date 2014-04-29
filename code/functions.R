require(devtools)

source_url("https://raw.github.com/thomaskelder/myRScripts/master/network-analysis/igraph.functions.R")
source_url("https://raw.github.com/thomaskelder/myRScripts/master/util/utils.R")
source_url("https://raw.github.com/thomaskelder/myRScripts/master/network-analysis/wgcna.functions.R")

performGoORA = function(set, reference = NULL, ontology = "BP", pvalueCutoff = 0.05, conditional = T) {
  if(is.null(reference)) {
    if(is.null(allMmEG_)) {
      allMmEG_ <<- unique(unlist(as.list(org.Mm.egGO2EG)))
      reference = allMmEG_
    } else {
      reference = allMmEG_
    }
  }
  
  params = new("GOHyperGParams", 
               geneIds = set, universeGeneIds = reference,
               ontology = ontology, pvalueCutoff = pvalueCutoff,
               conditional = conditional, annotation = "org.Mm.eg.db"
  )
  
  hasTerm = sapply(as.list(org.Mm.egGO)[set], function(id) {
    if(length(id)) sum(sapply(id, function(go) go["Ontology"] == ontology), na.rm = T) > 0
    else F
  })
  if(sum(hasTerm) > 0) {
    hyperGTest(params)
  } else {
    NULL
  }
}

getGeneInfo = function(entrezIds) {
  sym = as.character(as.list(org.Mm.egSYMBOL)[entrezIds])
  names(sym) = entrezIds
  nm = as.character(as.list(org.Mm.egGENENAME)[entrezIds])
  names(nm) = entrezIds
  list(symbols = sym, names = nm)
}

saturateValues = function(values, satMin = -10, satMax = 10) {
  values.sat = values
  t(apply(values, 1, function(x) {
    x[x < satMin] = satMin
    x[x > satMax] = satMax
    x
  }))
}
