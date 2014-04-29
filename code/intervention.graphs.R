####################################################################
## Intervention-specific networks for drugs treatments, including ##
## interactions from STITCH and WGCNA                             ##
####################################################################
source("code/vars.R")
source("code/loadNetworks.R")

## Remove unused attributes
vattr = list.vertex.attributes(graph)
for(d in c("09wkHighFat", "Atorvastatin", "Glibenclamide", 
					 "Metformine", "Pioglitazone", "Rosiglitazone",
					 "Salicylate", "Sitagliptin", "Vioxx")) {
	for(n in grep(d, vattr, value=T)) graph = remove.vertex.attribute(graph, n)
}

## Create a network for LifeStyle and enriched TFs
## Add LifeStyle as node to response graph and link to enriched TFs
tfs = read.delim(paste0(dataPath, "/IPA_LifeStyle_TFs.txt"), as.is=T, skip = 2)
rownames(tfs) = toupper(gsub(" \\(includes .+\\)$", "", tfs[,3]))
tfs = tfs[tfs[,"p.value.of.overlap"] < 0.001,]

allsymbols = as.list(org.Mm.egSYMBOL2EG)
names(allsymbols) = toupper(names(allsymbols))
allsymbols["CYP2C19"] = "13096"
allsymbols["CYP2F1"] = "13107"
allsymbols["CYP3A4"] = "13112"
allsymbols["CYP3A43"] = "56388"
allsymbols["AGTR1"] = "11607"
allsymbols["CYP2C18"] = "72082"

tf.targets = apply(tfs, 1, function(x) {
	tgts = strsplit(x["Target.molecules.in.dataset"], ",")[[1]]
	tgts = toupper(gsub(" \\(includes .+\\)$", "", tgts))
	ids = as.character(allsymbols[tgts])
	names(ids) = tgts
	ids
})
tf.ids = sapply(names(tf.targets), function(x) as.character(allsymbols[x]))
names(tf.targets) = tf.ids
tf.targets = tf.targets[names(tf.targets) != "NULL"]
tf.table = c("", "")
for(tf in names(tf.targets)) {
	for(tgt in tf.targets[[tf]]) {
		tf.table = rbind(tf.table, paste("L:", c(tf, tgt), sep=""))
	}
	tf.table = rbind(tf.table, c("LifeStyle", paste("L:", tf, sep="")))
}
tf.table = tf.table[2:nrow(tf.table),]
graph.lifestyle = graph.data.frame(as.data.frame(tf.table))
labels = as.character(as.list(org.Mm.egSYMBOL)[sub("L:", "", V(graph.lifestyle)$name)])
labels[labels == "NULL"] = sub("L:", "", V(graph.lifestyle)$name[labels == "NULL"])
V(graph.lifestyle)$label = labels

## Add the WGCNA results and STITCH interactions to the main graph
load(paste0(outPath, "wgcna/igraph.wgcna.RData"))
load(paste0(outPath, "stitch.graphs.RData"))
V(graph.wgcna)$name = paste("L:", V(graph.wgcna)$name, sep="")

## Only keep black, red, and yellow modules (based on previous selection)
graph.wgcna = induced.subgraph(graph.wgcna, V(graph.wgcna)[module %in% kwalksModules])

## Only keep edges within module
mods = V(graph.wgcna)$module
mods.edge = apply(get.edgelist(graph.wgcna, names=F), 1, function(e) {
  mods[e[1]] == mods[e[2]]
})
graph.wgcna = delete.edges(graph.wgcna, E(graph.wgcna)[!mods.edge])

tomerge = list(know = removeLonelyNodes(graph), coexp = graph.wgcna)
tomerge = append(tomerge, drugs.graphs)
tomerge = append(tomerge, list(LifeStyle = graph.lifestyle))

graph.resp = mergeGraphs(tomerge, T, T)
graph.resp = as.undirected(graph.resp)
E(graph.resp)$weight[is.na(E(graph.resp)$weight)] = 1
E(graph.resp)$wgcnaWeight = as.numeric(E(graph.resp)$weight)
graph.resp = remove.edge.attribute(graph.resp, "weight")
V(graph.resp)$module[is.na(V(graph.resp)$module)] = ""

graph.resp = dataToNodes(graph.resp, data.xref[,c(treatments.p, treatments.fdr, treatments.logfc)])

labels = as.character(as.list(org.Mm.egSYMBOL)[sub("L:", "", V(graph.resp)$name)])
labels[labels == "NULL"] = sub("L:", "", V(graph.resp)$name[labels == "NULL"])
write.table(cbind(id = V(graph.resp)$name, label = labels), paste(outPath, "gene.labels.txt", sep=""),
            row.names=F, sep="\t", quote=F)

save(graph.resp, file = paste0(outPath, "graph.combined.RData"))

responseGraph = function(attr, g, p.cutoff) {
  val = get.vertex.attribute(g, attr)
  val[is.na(val)] = 1
  
  ## Keep all nodes below p-value cutoff
  keep = val < p.cutoff
  
  ## Keep all drug interactions in case this is response network for drug
  d = sub("^.+_LDLR_", "", attr)
  d = sub("\\.LDLR_.+$", "", d)
  if(d %in% drugs) {
    keep.stitch = V(graph.resp)$name %in% V(drugs.graphs[[d]])$name
    keep = keep | keep.stitch
  }
  if(d == "LifeStyle") {
  	keep = keep | V(graph.resp)$name %in% c("LifeStyle", V(graph.lifestyle)[nei("LifeStyle")]$name)
  }
  subgraph(g, V(g)[keep])
}

## Response networks
filterAttr = treatments.p
rgraphs = lapply(filterAttr, responseGraph, graph.resp, p.max)
names(rgraphs) = filterAttr

## Network size
rgraphs.size = sapply(names(rgraphs), function(n) {
  g = rgraphs[[n]]
  gc = removeLonelyNodes(g)
  fc = get.vertex.attribute(gc, gsub("P.Value", "logFC", n))
  degcount = sum(data.xref[,n] <= p.max, na.rm = T)
  list(vcount = vcount(g), ecount = ecount(g),
       vcount.conn = vcount(gc), ecount.conn = ecount(gc),
       degcount = degcount, vup = sum(fc > 0, na.rm = T), vdown = sum(fc < 0, na.rm = T)
  )
})
write.table(rgraphs.size, file=paste(outPath, "response.graphs.sticth.wgcna.sizes.txt", sep=""), sep="\t", quote=F, col.names=NA)

rgraphs.treatment = rgraphs[2:length(rgraphs)]
rgraph.hf = rgraphs[[1]]

save(
  treatments, treatments.p, treatments.logfc, treatments.fdr, drugs, rgraphs,
  file = paste(outPath, "response.wgcna.stitch.RData", sep="")
  )

## Create google chart urls that show both HF and treatment FC in pie chart
## http://chart.apis.google.com/chart?chs=300x300&cht=p&chco=0000FF|FF0000&chd=t:50,50&chp=1.571
chartUrlBase = "http://chart.apis.google.com/chart?chs=100x100&cht=p&chd=t:50,50&chp=1.571&chf=bg,s,FFFFFF00"
chartColors = colorRamp(c("blue", "white", "red"))
lim = 2

charts = sapply(names(rgraphs.treatment), function(n) {
  print(sub("P\\.Value", "LogFC", n))
  values = cbind(
    get.vertex.attribute(graph.resp, logfcCols[[2]]),
    get.vertex.attribute(graph.resp, sub("P\\.Value", "logFC", n))
    )
  
  apply(values, 1, function(x) {
    xco = x
    xco[which(x > lim)] = lim
    xco[which(x < -lim)] = -lim
    if(is.na(sum(xco))) {
      ""
    } else {
      xco = rgb(chartColors((xco + lim) / (2*lim)), max=255)
      xco = gsub("#", "", xco)
      chco = paste(xco, collapse="|")
      paste(chartUrlBase, "&chco=", chco, sep="")
    }
  })
})
charts = cbind(V(graph.resp)$name, charts)
colnames(charts) = c("ID", sub("P\\.Value", "hf_treatment_chart", names(rgraphs.treatment)))
write.table(charts, paste(outPath, "/hf.vs.treatment.chart.txt", sep=""), sep="\t", quote=F, row.names=F)

## Find out if response networks are significantly enriched with disease associated genes in human
load(paste0(outPath, "associations.RData"))
all = rownames(data.xref) %in% paste("L:", names(associations.mouse),sep="")
names(all) = rownames(data.xref)
assoc = t(sapply(rgraphs, function(g) {
  set = all[intersect(V(g)$name, names(all))]
  p = performORA(set, all)
  c(all = 100 * sum(all) / length(all), resp = 100 * sum(set) / length(set), p = p)
}))
write.table(assoc, paste0(outPath, 'graph.disease.associations.enr.txt'), sep="\t", quote=F, col.names=NA)