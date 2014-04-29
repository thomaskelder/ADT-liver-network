library(digest)
source("code/vars.R")

load(paste0(outPath, "kwalks/kwalks.weighted.graphs.RData"))
load(paste0(outPath, "response.wgcna.stitch.RData"))
load(paste0(outPath, "associations.RData"))
source("code/loadNetworks.R")

## All identifiers possible to rank (all nodes in the reference network)
allIds = gsub("L:", "", intersect(V(graph)$name, rownames(data.xref)))

areaUnderCurve = function(y,x) {
  y = y - y[1]
  idx = 2:length(x)
  as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2
}

calcRoc = function(allIds, rankedIds) {
  as.data.frame(t(sapply(1:length(rankedIds), function(i) {
    pos = rankedIds[1:i]
    
    truePosIds = intersect(allIds, referenceIds)
    trueNegIds = setdiff(allIds, truePosIds)
    
    neg = setdiff(allIds, pos)
    hit = pos %in% truePosIds
    miss = !(truePosIds %in% pos)
    tp = sum(hit)
    fn = sum(miss)
    fp = sum(!hit)
    tn = sum(neg %in% trueNegIds)
    c(TPR = tp / (tp + fn), FPR = fp / (fp + tn))
  })))
}

calcCoverage = function(allIds, rankedIds) {
  as.data.frame(t(sapply(1:length(rankedIds), function(i) {
    pos = rankedIds[1:i]
    c(
      PercentDiseaseGenes = 100 * sum(pos %in% intersect(allIds, referenceIds)) / length(referenceIds), 
      NumberGenes = i
    )
  })))
}

plotBenchmarks = function(intervention, scores, modname = "", disease = "", filePrefix = paste0(outPath, "/roc/")) {
  rankedIds.kwalks = names(scores)[order(scores, decreasing=T, na.last=T)]
  deg.pvals = data.xref[paste0("L:", allIds), paste0("adj.P.Val_LDLR_", intervention, ".LDLR_16wkHighFat")]
  names(deg.pvals) = sub("L:", "", names(deg.pvals))
  
  rankedIds.deg = allIds[order(deg.pvals)]
  rankedIds.deg = rankedIds.deg[1:length(rankedIds.kwalks)]
  
  allIds.deg = rankedIds.deg
  allIds.kwalks = rankedIds.kwalks
  
  message("max TP for DEG: ", length(intersect(allIds.deg, referenceIds)))
  message("max TP for kwalks: ", length(intersect(allIds.kwalks, referenceIds)))
  
  rnd.n = 10
  
  roc = list()
  roc$kwalks = calcRoc(allIds.kwalks, rankedIds.kwalks)
  roc$deg = calcRoc(allIds.deg, rankedIds.deg)
  roc.rnd = as.data.frame(lapply(1:rnd.n, function(x) {
    rnd = sample(allIds.kwalks, length(rankedIds.kwalks))
    calcRoc(allIds.kwalks, rnd)
  }))
  roc$random = cbind(
    TPR = rowMeans(roc.rnd[,grep("TPR", colnames(roc.rnd))]),
    FPR = rowMeans(roc.rnd[,grep("FPR", colnames(roc.rnd))])
  )
  
  cov = list()
  cov$kwalks = calcCoverage(allIds.kwalks, rankedIds.kwalks)
  cov$deg = calcCoverage(allIds.deg, rankedIds.deg)
  cov.rnd = as.data.frame(lapply(1:rnd.n, function(x) {
    rnd = sample(allIds, length(rankedIds.kwalks))
    calcCoverage(allIds.kwalks, rnd)
  }))
  cov$random = cbind(
    PercentDiseaseGenes = rowMeans(cov.rnd[,grep("PercentDiseaseGenes", colnames(cov.rnd))]),
    NumberGenes = rowMeans(cov.rnd[,grep("NumberGenes", colnames(cov.rnd))])
  )
  
  library(reshape)
  library(ggplot2)
  
  rocTbl = melt(as.data.frame(roc))
  rocTbl = cbind(rocTbl[grep("TPR", rocTbl$variable),], rocTbl[grep("FPR", rocTbl$variable), "value"])
  colnames(rocTbl) = c("Ranking", "TPR", "FPR")
  rocTbl$Ranking = gsub("\\..+PR$", "", rocTbl$Ranking)
  
  area = sapply(unique(rocTbl$Ranking), function(x) {
    areaUnderCurve(rocTbl[rocTbl$Ranking == x, 2], rocTbl[rocTbl$Ranking == x, 3])
  })
  
  svg(paste0(filePrefix, "roc.", disease, ".", intervention, ".svg"))
  p = ggplot(rocTbl, aes(FPR, TPR))
  p = p + geom_point(aes(colour = factor(Ranking)))
  p = p + ggtitle(paste0("ROC for identifying genes\nassociated with ", disease, "\nfor signature linking module '", modname, "' and ", intervention))
  p = p + annotate(
    "text", x = 0.1, y = 0.9, hjust = 0, size = 4,
    label = paste(c("AUC", paste(names(area), round(area, 2), sep=": ")), collapse="\n"))
  print(p)
  dev.off()
  
  covTbl = melt(as.data.frame(cov))
  covTbl = cbind(covTbl[grep("PercentDiseaseGenes", covTbl$variable),], covTbl[grep("NumberGenes", covTbl$variable), "value"])
  colnames(covTbl) = c("Ranking", "PercentDiseaseGenes", "NumberGenes")
  covTbl$Ranking = gsub("\\.PercentDiseaseGenes$", "", covTbl$Ranking)
  test <<- covTbl
  svg(paste0(filePrefix, "coverage.", disease, ".", intervention, ".svg"))
  p = ggplot(covTbl, aes(NumberGenes, PercentDiseaseGenes))
  p = p + geom_point(aes(colour = factor(Ranking)))
  p = p + ggtitle(paste0("Coverage of genes\nassociated with ", disease, "\n vs size of signature for subgraph linking module '", modname, "' and ", intervention))
  p = p + xlab("Number of genes") + ylab("% known disease genes")
  print(p)
  dev.off()
}

singleScores = function(intervention, module) {
  g = kwalks.graphs[[intervention]][[module]]$graph
  scores = get.vertex.attribute(g, paste0("kwalks.", intervention, ".", module))
  names(scores) = gsub("L:", "", V(g)$name)
  scores = scores[names(scores) != intervention]
  scores = scores[!is.na(scores)]
  #scores[is.na(scores)] = 0
  scores
}

combinedScores = function(intervention, modules) {
  ids = unique(unlist(lapply(kwalks.graphs[[intervention]], function(x) gsub("L:", "", V(x$graph)$name))))
  scoreTable = sapply(modules, function(m) {
    g = kwalks.graphs[[intervention]][[m]]$graph
    scores = get.vertex.attribute(g, paste0("kwalks.", intervention, ".", m))
    names(scores) = gsub("L:", "", V(g)$name)
    v = rep(NA, length(ids))
    names(v) = ids
    v[names(scores)] = scores
    v
  })
  maxScores = apply(scoreTable, 1, max)
  #maxScores = maxScores[names(maxScores) != intervention]
  maxScores = maxScores[!is.na(maxScores)]
  maxScores[is.na(maxScores)] = 0
  maxScores
}

diseases = list(
  t2dm = names(associations.mouse)[grep("t2dm or metsyn", associations.mouse)],
  t2dm.liver = names(associations.mouse)[grep("t2dm and liver", associations.mouse)],
  genetic = names(associations.mouse)[grep("genetic", associations.mouse)],
  all = names(associations.mouse)
)
for(dn in names(diseases)) {
  for(intervention in c("LifeStyle", "T0901317", "Fenofibrate")) {
    referenceIds = diseases[[dn]]
    message(intervention)
    plotBenchmarks(intervention, combinedScores(intervention, c("red", "black", "yellow")), "all", dn)
  }
}

## Find out if signatures are significantly enriched with disease associated genes
doOra = function(set, ref) {
  p = performORA(set, ref)
  r = c(ref = 100 * sum(ref) / length(ref), sig = 100 * sum(set) / length(set), p = p)
  r['fold'] = r['sig'] / r['ref']
  r
}

sig.ora = list()
for(dn in names(diseases)) {
  sig.ora[[dn]] = list()
  for(intervention in c("LifeStyle", "T0901317", "Fenofibrate")) {
    referenceIds = diseases[[dn]]
    scores = combinedScores(intervention, c("red", "black", "yellow"))
    
    ref = gsub("L:", "", rownames(data.xref)) %in% referenceIds
    names(ref) = gsub("L:", "", rownames(data.xref))
    r.kwalks = doOra(ref[intersect(names(ref), names(scores))], ref)

    deg.pvals = data.xref[paste0("L:", allIds), paste0("adj.P.Val_LDLR_", intervention, ".LDLR_16wkHighFat")]
    names(deg.pvals) = sub("L:", "", names(deg.pvals))
    r.deg = doOra(ref[allIds[order(deg.pvals)[1:length(scores)]]], ref)
    
    sig.ora[[dn]][[intervention]] = list(kwalks = r.kwalks, deg = r.deg)
  }
}
sig.ora = t(as.data.frame(sig.ora))
write.table(sig.ora, paste0(outPath, 'roc/sig.ora.txt'), quote=F, sep="\t")