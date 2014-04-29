####################################################################
## Perform kWalks algorithm on intervention-specific networks     ##
####################################################################
source("code/vars.R")
source("code/kwalks.functions.R")

saveGraphs = T
myOutPath = paste0(outPath, 'kwalks/')

load(paste0(outPath, "response.wgcna.stitch.RData"))
load(paste0(outPath, "wgcna/wgcna.modmembership.RData"))

useWeights = T

kwalks.treatments = c(LifeStyle = "LifeStyle", drugs)

iterations = 3
## Run kwalks for each dugs - module combination
kwalks.graphs = list()
for(d in kwalks.treatments) {
  message(d)
  g = rgraphs[[grep(d, names(rgraphs))]]
  V(g)$module[is.na(V(g)$module)] = ""
  
  ## Get largest connected component
  comp = decompose.graph(g, "weak", min.vertices = 4)
  gc = comp[order(sapply(comp, vcount), decreasing=T)][[1]]
  gc = simplify(gc, edge.attr.comb=list(wgcnaWeight="sum", "concat"))
  if(useWeights) E(gc)$weight = E(gc)$wgcnaWeight
  
  modules = unique(V(g)$module)
  modules = modules[modules %in% kwalksModules]
  
  smx.w = kwalkMatrix(gc, attr = "weight")

  from = V(gc)[d]
  
  kwalks.graphs[[d]] = list()
  for(m in modules) {
    message(m)
    to =  V(gc)[module == m]
    
    kw.attr = paste("kwalks", d, m, sep=".")
    g.kw = runKwalk(
      gc, smx.w, from, to, attr = kw.attr, iterations = iterations,
      tmpPath = paste("/tmp/", ifelse(useWeights, "w_", "nw_"), sep=""), scriptPath = 'code/')
    kwalks.graphs[[d]][[m]] = list(graph = g.kw)
  }
}
w = ifelse(useWeights, "weighted.", "")
save(kwalks.graphs, file = paste(myOutPath, "kwalks.", w, "graphs.RData", sep=""))

kwalks.graphs.vis = list()
for(d in names(kwalks.graphs)) {
	message(d)
	mgraphs = kwalks.graphs[[d]]
	kwalks.graphs.vis[[d]] = list()
	for(m in names(mgraphs)) {
		message("\t", m)
		g.kw = mgraphs[[m]]$graph
		kw.attr = paste("kwalks", d, m, sep=".")
		
		erel = get.edge.attribute(g.kw, kw.attr)
		incr = seq(max(erel, na.rm=T), 0.005, -0.005)
		incr = c(incr, seq(min(incr) - 0.0001, 0, -0.0001))
		# Maximum cutoff that keeps kwalks.treatments and at least x of module
		# nodes connected
		inmodule = V(g.kw)[module == m]$name
		#minInMod = c(length(inmodule), 5:1)
		minInMod = 5:1
		for(i in minInMod) {
			message("Trying to find cutoff for ", i, " module nodes")
			cutoff = findEdgeRelevanceCutoff(g.kw, V(g.kw)[d]$name, inmodule, i, attr = kw.attr, incr = incr)
			cutoff.value = cutoff$cutoff
			if(cutoff.value > 0) break
			minInMod = minInMod - 1
		}
		g.cut = cutoff$graph
		if(cutoff.value <= 0) {
			cutoff.value = quantile(erel, probs = 0.999)
			g.cut = doCutoff(g.kw, erel, cutoff.value)
			message("Could not reach module nodes with cutoff, using fixed cutoff")
		}
		
		## Create cutoff network including top module genes + edges
		rg = rgraphs[[grep(d, names(rgraphs))]]
		rg = induced.subgraph(rg, V(rg)[module == m])
    mm = modMembershipGenes[gsub("L:", "", V(rg)$name), "modMembershipR"]
    mm.keep = paste("L:", names(mm)[order(abs(mm), decreasing=T)[1:10]], sep="")
    
    keep.nodes = union(mm.keep, V(g.cut)$name)
    keep.nodes = intersect(keep.nodes, V(g.kw)$name)
    
		g.cut.mod = subgraph(g.kw, V(g.kw)[keep.nodes])
		erel.mod = get.edge.attribute(g.cut.mod, kw.attr)
		edges.rm = as.numeric(E(g.cut.mod)[erel.mod < cutoff.value])
		retain = apply(get.edgelist(g.cut.mod)[edges.rm,], 1, function(x) {
			sum(x %in% mm.keep) > 0
		})
		g.cut.mod = delete.edges(g.cut.mod, edges.rm[!retain])
		
		## Add attribute to distinguish added edges from kwalks edges
		g.cut.mod = set.edge.attribute(
			g.cut.mod, gsub("kwalks", "extra", kw.attr),
			value = get.edge.attribute(g.cut.mod, kw.attr) < cutoff.value
		)
		
		## Add attribute to mark kwalks.treatments and first neighbors
		V(g.cut.mod)$isDrug = "false"
		V(g.cut.mod)[d]$isDrug = "true"
		V(g.cut.mod)$isDrugNeighbor = "false"
		V(g.cut.mod)[nei(d)]$isDrugNeighbor = "true"
		
    V(g.cut.mod)[module == m]$moduleMembership = modMembershipGenes[gsub("L:", "", V(g.cut.mod)[module == m]$name), "modMembershipR"]
				
		message("cutoff: ", cutoff.value)
		message("nodes: ", vcount(g.cut))
		message("edges: ", ecount(g.cut))
		message("nodes + mod: ", vcount(g.cut.mod))
		message("edges + mod: ", ecount(g.cut.mod))
		
		g.cut = pruneBranches(g.cut, c(inmodule, V(g.kw)[d]$name))
		g.cut.mod = pruneBranches(g.cut.mod, c(inmodule, V(g.kw)[d]$name))
		message("nodes after pruning: ", vcount(g.cut))
		message("nodes after pruning + mod: ", vcount(g.cut.mod))
		
		kwalks.graphs.vis[[d]][[m]] = list(
			cutoff = cutoff.value, graph = g.cut,
			graph.mod = g.cut.mod
		)
    	saveGML(g.cut.mod, paste(myOutPath, "kwalks.", w, d, ".", m, ".cutoff.gml", sep=""), paste("kwalks", d, m, cutoff.value, sep="."))
	}
}

save.image(paste(myOutPath, "kwalks.total.RData", sep=""))

## Attribute file for included module nodes relevance/kwalks scores
for(d in names(kwalks.graphs.vis)) {
	r.d = kwalks.graphs.vis[[d]]
	for(m in names(r.d)) {
		r = r.d[[m]]
		kw.attr = paste("kwalks", d, m, sep=".")
		cutoff = r$cutoff
		g.cut = r$graph
		g.mod = r$graph.mod
		attr = paste("kwalksFromModule", d, m, sep=".")
		values = rep("module.none", vcount(g.mod))
		names(values) = V(g.mod)$label
		values[V(g.mod)$module == m] = "module.above.cutoff"
		
		values[!(V(g.mod)$name %in% V(g.cut)$name)] = "module.below.cutoff"
		
		write.table(c(attr, paste(names(values), "=", values)), paste(myOutPath, attr, ".txt", sep=""), col.names=F, row.names=F, quote=F)
	}
}

## Create table of node relevance scores
kwalks.tbls = list()

kwalksModules = c("yellow", "red", "black") #Different order for figure

nodes = unique(unlist(sapply(rgraphs, function(rg) V(rg)$name)))
node2mod = unlist(sapply(rgraphs, function(rg) V(rg)$module))
names(node2mod) = unlist(sapply(rgraphs, function(rg) V(rg)$name))
nodeInfo = getGeneInfo(sub("L:", "", nodes))

drugNbs = lapply(kwalks.treatments, function(d) {
	V(rgraphs[[grep(d, names(rgraphs))]])[nei(d)]$name
})

for(d in kwalks.treatments) {
	g = rgraphs[[grep(d, names(rgraphs))]]
	
	tbl = sapply(kwalksModules, function(m) {
		gk = kwalks.graphs[[d]][[m]]$graph
		kw.attr = paste("kwalks", d, m, sep=".")
		
		nrel = get.vertex.attribute(gk, kw.attr)
		r = rep(NA, length(nodes))
		names(r) = nodes
		r[V(gk)$name] = nrel
		r
	})
	rownames(tbl) = nodes
	
	kwalks.tbls[[d]] = tbl

	tbl = cbind(
		entrez = nodes,
		symbol = nodeInfo$symbols, name = nodeInfo$names,
		tbl)
	
	na.rows = rowSums(is.na(tbl)) == ncol(tbl) - 3 
	write.table(tbl[!na.rows,], 
		paste(myOutPath, "kwalks.table.", d, ".txt", sep=""), 
		sep="\t", row.names=F, quote=F)
}

## Create heatmap to summarize top nodes
tnr = 10
tn = unique(unlist(lapply(kwalks.tbls, extractTop, nr = tnr)))

kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.pdf")
kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.pdf", colorByDEG = T)

kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.group.pdf", colorByDEG = T, firstGroup = unique(unlist(extractTop(kwalks.tbls[["LifeStyle"]], tnr))))

kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.DLI.pdf", colorByDEG = T, firstGroup = unique(unlist(extractTop(kwalks.tbls[["LifeStyle"]], tnr))))
kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.FF.pdf", colorByDEG = T, firstGroup = unique(unlist(extractTop(kwalks.tbls[["Fenofibrate"]], tnr))))
kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.T09.pdf", colorByDEG = T, firstGroup = unique(unlist(extractTop(kwalks.tbls[["T0901317"]], tnr))))

kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.DLI.ranked.pdf", colorByDEG = T, rankSort = T, firstGroupCol = "LifeStyle", firstGroup = unique(unlist(extractTop(kwalks.tbls[["LifeStyle"]], tnr))))
kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.FF.ranked.pdf", colorByDEG = T, rankSort = T, firstGroupCol = "Fenofibrate", firstGroup = unique(unlist(extractTop(kwalks.tbls[["Fenofibrate"]], tnr))))
kwalksHeatmap.bygroup(tn, kwalks.tbls, "heatmap.kwalks.path.color.T09.ranked.pdf", colorByDEG = T, rankSort = T, firstGroupCol = "T0901317", firstGroup = unique(unlist(extractTop(kwalks.tbls[["T0901317"]], tnr))))

tn.target = unique(unlist(lapply(kwalks.tbls, function(tbl) {
	top = list()
	for(cn in colnames(tbl)) {
		mod = names(which(node2mod == cn, useNames = T))
		col = tbl[rownames(tbl) %in% mod, cn]
		top[[cn]] = names(col)[order(col, decreasing=T)[1:5]]
	}
	top
})))

#Sort by module
tn.target = tn.target[order(node2mod[tn.target])]
kwalksHeatmap(tn.target, kwalks.tbls, "heatmap.kwalks.target.pdf", cluster_rows = F, cluster_cols = F)

## Number of nodes not a treatment target or module gene
rby = c("black", "red", "yellow")
ntn = matrix(NA, nrow = length(kwalks.graphs.vis), ncol = length(rby))
colnames(ntn) = rby
rownames(ntn) = names(kwalks.graphs.vis)
ntn.r = ntn

for(m in rby) {
	for(tr in names(kwalks.graphs.vis)) {
		g = kwalks.graphs.vis[[tr]][[m]]$graph.mod
		tn = c(tr, V(g)[nei(tr)]$name)
		mn = V(g)[module == m]$name
		r = setdiff(V(g)$name, c(tn, mn))
		ntn[tr, m] = length(r)
		ntn.r[tr, m] = length(r) / length(tn)
	}
}