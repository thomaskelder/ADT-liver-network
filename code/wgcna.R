###########################################
## WGCNA on 10OAD liver transcriptiomics ##
###########################################
source("code/vars.R")
source("code/treatmentEffect.R")

myOutPath = paste0(outPath, "wgcna/", sep="")

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

## Calculate on data from all mice, except wk9
mice = names(mouse.trt[mouse.trt %in% c("chow", "HF - 16wk", "Lifestyle (chow)", "Fenofibrate", "T0901317")])
datExpr = expr.liver.data[, colnames(expr.liver.data) %in% mice]
datExpr = avereps(datExpr, ID=expr.liver.data.genes[,3])
datExpr = datExpr[!is.na(rownames(datExpr)),]
datExpr = t(datExpr)

## Pick a soft threshold
Rs = 0.9
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))
#sft = pickSoftThreshold(datExpr, verbose = 5, powerVector = powers, RsquaredCut = Rs)

## Evaluate soft threshold
#plotSoftThreshold(sft, Rs, file = paste(myOutPath, "softThreshold.pdf", sep=""))

softPower = 4
minModuleSize = 20

## Create the co-expression network
tomFileBase = paste0(myOutPath, "wgcna.tom")
net = blockwiseModules(datExpr, power = softPower,
   minModuleSize = minModuleSize, pamStage = FALSE,
   numericLabels = TRUE, pamRespectsDendro = FALSE,
   saveTOMs = TRUE, maxBlockSize = ncol(datExpr),
   saveTOMFileBase = tomFileBase,
   verbose = 5)
   tryCatch({
     plotNetAsDendro(net, file = paste0(myOutPath, "dendro.pdf"))
   },  error = function(e) warning(e))


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = net$MEs
geneTree = net$dendrograms[[1]]

moduleCounts = sort(table(moduleColors))
write.table(moduleCounts, paste0(myOutPath, "module.counts.txt"), quote=F, sep="\t")

## Check PCA scores of eigengens
MEdata = moduleEigengenes(datExpr, moduleColors, softPower = softPower)
varExpl = t(MEdata$varExplained)
rownames(varExpl) = colnames(MEdata$eigengenes)
write.table(varExpl, file = paste(myOutPath, "eigengene.PC1.scores.txt", sep=""), sep = "\t", quote = F, col.names=F)
minVarExpl = 0.5
validMods = rownames(varExpl)[varExpl >= minVarExpl]

## Correlate module eigengenes with disease endpoints
datTraits = t(ypar.data[inclPars, rownames(datExpr)])

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0[,validMods])
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

moduleTraitCor.sp = moduleTraitCor
moduleTraitPvalue.sp = moduleTraitPvalue
for(mod in colnames(MEs)) {
	for(trait in colnames(datTraits)) {
		r = cor.test(MEs[,mod], datTraits[,trait], method="spearman")
		moduleTraitCor.sp[mod, trait] = r$estimate
		moduleTraitPvalue.sp[mod, trait] = r$p.value
	}
}

tryCatch({
     plotTraitCorrelations(moduleTraitCor, moduleTraitPvalue, MEs, datTraits, 
       file = paste(myOutPath, "traitCor.pdf", sep=""))
   },  error = function(e) warning(e))

tryCatch({
     plotTraitCorrelations(moduleTraitCor.sp, moduleTraitPvalue.sp, MEs, datTraits, 
       file = paste(myOutPath, "traitCor.sp.pdf", sep=""))
   },  error = function(e) warning(e))

plotTraitMod = function(trait, mod, r = "", title = "", xlab = mod, ylab = trait) {
	library(ggplot2)
	x = datTraits[,trait]
	y = MEs[,mod]
	group = mouse.trt[rownames(datTraits)]
	dta = data.frame(
		mouse = rownames(datTraits), trait = x, eg = y, group = group
	)
	
	p = ggplot(dta, aes(eg, trait))
	p = p + geom_point(aes(colour = group), size = 4)
	#p = p + geom_text(aes(label = mouse), hjust=0, vjust=0, size=3)
	p = p + opts(title = title)
	p = p + xlab(xlab) + ylab(ylab)
	
	pdf(paste(myOutPath, "cor.", r, ".", trait, ".", mod, ".pdf", sep=""))
	print(p)
	dev.off()
}

topCors = c("trait", "module", "r", "p")
for(mod in colnames(MEs)) {
	for(trait in colnames(datTraits)) {
		r = moduleTraitCor[mod, trait]
    p = moduleTraitPvalue[mod, trait]
		if(abs(r) > 0.75) {
			plotTraitMod(trait, mod, r = round(r, 2), title = paste("r =", r))
			topCors = rbind(topCors, c(trait, sub("ME", "", mod), r, p))
		}
	}
}
colnames(topCors) = topCors[1,]
topCors = topCors[2:nrow(topCors),]
write.table(topCors, paste(myOutPath, "topTraitModuleCors.txt", sep=""), sep="\t", quote=F, row.names=F)
plotTraitMod("glucose", "MEdarkturquoise", r = "-0.54")
plotTraitMod("glucose", "MEdarkgreen", r = "-0.51")
plotTraitMod("insulin", "MEpurple", r = "-0.56")
plotTraitMod("QUICKI", "MEblue", r = "-0.44")
plotTraitMod("QUICKI", "MEgreenyellow", r = "0.35")
plotTraitMod("cholesterol", "MEpurple", r = "-0.56")

## Calculate module membership score for each gene (correlation to module eigengene)
modMembership = as.matrix(cor(datExpr, MEs, use = 'p'))
modMembershipP = as.matrix(corPvalueStudent(modMembership, nrow(datExpr)))
names(moduleColors) = colnames(datExpr)
modMembershipGenes = t(sapply(rownames(modMembership), function(g) {
  m = paste("ME", moduleColors[g], sep="")
  if(!(m %in% c("MENA", "MEgrey")) & m %in% colnames(modMembership)) {
    c(modMembership[g, m], modMembershipP[g, m])
  } else {
    c(NA,NA)
  }
}))
colnames(modMembershipGenes) = c("modMembershipR", "modMembershipP")
write.table(
  cbind(symbol = as.character(mget(rownames(modMembershipGenes), org.Mm.egSYMBOL, ifnotfound = NA)), modMembershipGenes, module = moduleColors),
  file = paste(myOutPath, "wgcna.modmembership.txt", sep=""), sep="\t", quote=F, col.names=NA)
save(modMembershipGenes, file = paste(myOutPath, "wgcna.modmembership.RData", sep=""))

## Perform GO enrichment analysis
go = GOenrichmentAnalysis(moduleColors, colnames(datExpr), organism="mouse")
go.tab = go$bestPTerms[[4]]$enrichment
write.table(go.tab, file = paste(myOutPath, "wgcna.go.txt", sep=""), sep = "\t", quote = F, row.names = F)

library(GO.db)
term.info = as.list(GOTERM)
goRef = unique(unlist(as.list(org.Mm.egGO2EG)))
goRef = intersect(goRef, colnames(datExpr))
go.res = lapply(unique(moduleColors), function(m) {
  mg = colnames(datExpr)[moduleColors == m]
  print(mg)
  bp = performGoORA(mg, reference = goRef, ontology = "BP", pvalueCutoff = 0.01, conditional = F)
  mf = performGoORA(mg, reference = goRef, ontology = "MF", pvalueCutoff = 0.01, conditional = F)
  list(bp = bp, mf = mf)
})
names(go.res) = unique(moduleColors)
processGORes = function(r, correct = T) {
  fdr = p.adjust(pvalues(r), method = "BH")
  cbind(
    GOID = names(pvalues(r)),
    term = sapply(term.info[names(pvalues(r))], Term),
    count = geneCounts(r),
    size = universeCounts(r),
    expected = expectedCounts(r),
    pvalue = pvalues(r),
    fdr.pvalue = fdr
  )
}

go.res.tab = NULL
showMin = 5
showMin.p = 1E-4
for(n in names(go.res)) {
  bp = processGORes(go.res[[n]]$bp)
  bp = bp[order(bp[,"fdr.pvalue"]),]
  bp.fdr = as.numeric(bp[,"fdr.pvalue"])
  if(sum(bp.fdr < showMin.p) > showMin) bp = bp[bp.fdr < showMin.p,]
  else bp = bp[1:5,]
  
  mf = processGORes(go.res[[n]]$mf)
  mf.fdr = as.numeric(mf[,"fdr.pvalue"])
  if(sum(mf.fdr < showMin.p) > showMin) mf = mf[mf.fdr < showMin.p,]
  else mf = mf[1:5,]
  
  tab = rbind(bp, mf)
  tab = cbind(module = rep(n, nrow(tab)), moduleSize = rep(sum(moduleColors == n), nrow(tab)), tab)
  if(is.null(go.res.tab)) go.res.tab = tab
  else go.res.tab = rbind(go.res.tab, tab)
}
write.table(go.res.tab, file = paste(myOutPath, "wgcna.go.gostats.txt", sep=""), sep = "\t", quote = F, row.names = F)

## Consolidate into single row per module
bestGO = t(sapply(names(go.res), function(mod) {
  print(mod)
  modGO = which(go.res.tab[,'module'] == mod)
  if(length(modGO) == 0) {
    c("", "", -1)
  } else {
    sel = go.res.tab[modGO, ]
    sig = rbind(sel[as.numeric(sel[,"fdr.pvalue"]) < 1E-2,])
    
    c(
      paste(sig[,"term"], collapse = "; "), 
      paste(sig[,"count"], collapse = "; "), 
      as.numeric(nrow(sig) > 0)
      )
  }
}))
colnames(bestGO) = c("GO Term", "Nr module genes in term", "Significant annotation")
write.table(bestGO[sub("ME", "", validMods),], paste(myOutPath, "wgcna.best.go.txt", sep=""), sep="\t", quote=F, col.names=NA)

save.image(paste(myOutPath, "wgcna.RData", sep=""))

## Print summary table containing module memberships and DEGs
wsummary = cbind(
  id = rownames(modMembershipGenes),
  symbol = as.character(mget(rownames(modMembershipGenes), org.Mm.egSYMBOL, ifnotfound = NA)), 
  module = moduleColors,
  modMembershipGenes
)

for(trt in treatments) {
  sig = data.xref[, grep(paste("P\\.Value_LDLR_", trt, sep=""), colnames(data.xref))] < p.max
  fc = data.xref[, grep(paste("logFC_LDLR_", trt, sep=""), colnames(data.xref))]
  sel = paste("L:", wsummary[,'id'], sep="")
  d = cbind(fc[sel], sig[sel])
  colnames(d) = paste(trt, c("log2FC", "significant"), sep="_")
  wsummary = cbind(wsummary, d)
}
write.table(wsummary, file=paste(myOutPath, "module.deg.summary.txt", sep=""), sep="\t", quote=F, row.names=F)

save.image(paste(myOutPath, "wgcna.RData", sep=""))

## Convert adjacency / TOM into igraph network
library(igraph)

inclModules = gsub("ME", "", validMods)
inclGenes = moduleColors %in% inclModules

nodeIDs = colnames(datExpr)[inclGenes]
id2label = expr.liver.data.genes[,2]
names(id2label) = expr.liver.data.genes[,3]
nodeLabels = id2label[nodeIDs]

adj = adjacency(datExpr[,inclGenes], power = softPower)
graph = graph.adjacency(adj, mode = "undirected", diag = F, weighted = T)
V(graph)$name = nodeIDs
V(graph)$label = nodeLabels
V(graph)$module = moduleColors[inclGenes]

## Save the graph for use in further network analysis
graph.wgcna = removeLonelyNodes(graph)
save(graph.wgcna, file = paste(myOutPath, "igraph.wgcna.RData", sep=""))

## Remove low weight edges to reduce network file size
graph.wgcna.filt = delete.edges(graph.wgcna, E(graph.wgcna)[weight < 0.25])

## Add data to graph
V(graph.wgcna.filt)$entrez = V(graph.wgcna.filt)$name
V(graph.wgcna.filt)$name = paste("L:", V(graph.wgcna.filt)$name, sep = "")
graph.wgcna.filt = dataToGraph(graph.wgcna.filt, data, treatments.p, "xref", edges = F, nodes = T)
graph.wgcna.filt = dataToGraph(graph.wgcna.filt, data, treatments.fdr, "xref", edges = F, nodes = T)
graph.wgcna.filt = dataToGraph(graph.wgcna.filt, data, treatments.logfc, "xref", edges = F, nodes = T)
saveGML(graph.wgcna.filt, paste(myOutPath, "wgcna.modules.gml", sep=""), title = "WGCNA modules")

## Find out which modules are overrepresented with DEGs
deg.sets = getAmplifiedReversed(p.max, "P.Value", p.max)

module.deg = NULL
for(m in unique(V(graph.wgcna)$module)) {
  modGenes = V(graph.wgcna)[module == m]$name
  
  ora = t(sapply(treatments, function(trt) {
    ref = data.xref[, grep(paste("P\\.Value_LDLR_", trt, sep=""), colnames(data.xref))] < p.max
    set = ref[paste("L:", modGenes, sep="")]
    rev = amp = add = ""
    if(trt != "16wkHighFat") {
    col = paste("LDLR_", trt, ".LDLR_16wkHighFat", sep="")
      amp = sum(deg.sets$amp[names(set), col])
      rev = sum(deg.sets$rev[names(set), col])
      add = sum(deg.sets$add[names(set), col])
    }
    c(
      treatment = trt,
      moduleSize = length(set),
      moduleDeg = sum(set),
      pvalue = performORA(set, ref),
      reversed = rev,
      amplified = amp,
      additional = add,
      `reversed %` = ifelse(rev == "", "", 100*rev / sum(set)),
      `amplified %` = ifelse(amp == "", "", 100*amp / sum(set)),
      `additional %` = ifelse(add == "", "", 100*add / sum(set))
      )
  }))
  
  mrow = c(m, rep("", 9))
  if(is.null(module.deg)) module.deg = rbind(module = mrow, ora)
  else module.deg = rbind(module.deg, module = mrow, ora)
}
write.table(module.deg, file=paste(myOutPath, "module.deg.txt", sep=""), sep="\t", quote=F, row.names=F)
save.image(paste(myOutPath, "wgcna.all.RData", sep=""))
