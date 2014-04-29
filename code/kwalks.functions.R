pruneBranches = function(g, ignore = c()) {
	## Iteratively remove all nodes with only one connection
	## until there are none left
	checkDegree = function(g) {
		d = igraph::degree(g)
		dns = V(g)[d == 1]$name
		setdiff(dns, ignore)
	}
	rm = checkDegree(g)
	while(length(rm) > 0) {
		g = subgraph(g, V(g)[!(name %in% rm)])
		rm = checkDegree(g)
	}
	g
}

doCutoff = function(g, er, cutoff, removeNodes = T) {
	gc = delete.edges(g, E(g)[er < cutoff])
	if(removeNodes) removeLonelyNodes(gc)
	else gc
}

findEdgeRelevanceCutoff = function(g.kw, toInclude, optional, minOptional = 1, attr = "kwalks", incr = 0.001) {
  # Determine minimum edge relevance so that
  # both the drugs and module nodes are in same
  # component
  er = get.edge.attribute(g.kw, attr)
  
  g.kw.cut = NULL
  if(length(incr) == 1) incr = seq(max(er, na.rm=T), 0, -incr)
  for(cutoff in incr) {
    g.kw.cut = doCutoff(g.kw, er, cutoff, F)
    cl = clusters(g.kw.cut, mode = "weak")$membership
    names(cl) = V(g.kw.cut)$name
    cl.incl = cl[toInclude]
    cl.opt = cl[optional]
    if(length(unique(cl.incl)) == 1) {
    	if(sum(cl.opt == unique(cl.incl)) >= minOptional) break
    }
  }
  g.kw.cut = removeLonelyNodes(g.kw.cut)
  list(graph = g.kw.cut, cutoff = cutoff)
}

## Convert graph to a format readable for kwalks
kwalkMatrix = function(g, file, ...) {
  adj = get.adjacency(g, ...)
  smx = sapply(1:nrow(adj), function(n) {
    nbs = which(adj[n,] > 0)
    if(length(nbs > 0)) {
      e = paste(nbs, sprintf("%.32f", adj[n,nbs]), sep=":")
      paste(n, paste(e, collapse = " "))
    } else {
      n
    }
  })
  smx = c(vcount(g), smx)
}

runKwalk = function(g, smx, a, b, lmax = 50, tmpPath = "/tmp/", scriptPath = getwd(), attr = "kwalks", iterations = 1) {
  message("Running kwalks from ", a, " (", V(g)[a]$name, ") to ", b, " (", V(g)[b]$name, ")")
  
  # Write the graph file
  fin = paste(tmpPath, "kwalks.graph.txt", sep="")
  write.table(smx, fin, quote=F, row.names = F, col.names=F)
  
  fout = paste(tmpPath, "kwalks.out", sep="")

  # Run kwalks
  for(i in seq_len(iterations)) {
	  rel = paste(paste(a, collapse=":"), paste(b, collapse=":"), sep="#")
	  cmd = paste(scriptPath, "/lkwalk", sep="")
	  cmd = paste(cmd, "-g", fin, "-l", lmax, "-k", rel, "-o", fout)
	  message(cmd)
	  system(cmd)
	  
	  fin = paste(fout, ".dif", sep="")
  }
  message("Kwalks finished!")
  # Parse results
  g.w = g
  
  fout.dif = paste(fout, ".dif", sep="")
  dif = unlist(read.delim(fout.dif, sep="\t", stringsAsFactor = F, header=F))
  dif = dif[2:length(dif)]
  
  edges = c()
  values = c()
  
  weights = matrix(0, nrow = vcount(g.w), ncol = vcount(g.w))
  
  for(i in 1:length(dif)) {
    if(i %% 100 == 0) {
      message(i, " out of ", length(dif))
    }
    x = strsplit(dif[i], " ")[[1]]
    n1 = as.numeric(x[1])
    for(y in x[2:length(x)]) {
      y = as.numeric(strsplit(y, ":")[[1]])
      n2 = y[1]
      w = y[2]
      weights[n1, n2] = weights[n2, n1] = w
    }
  }
  values = apply(get.edgelist(g.w, names = F), 1, function(x) {
    weights[x[1], x[2]]
  })
  g.w = set.edge.attribute(g.w, attr, E(g.w), values)
  nw = as.character(unlist(read.delim(paste(fout, ".N", sep=""), stringsAsFactor = F, header = F)))
  nw = sapply(nw, function(x) as.numeric(strsplit(x, ":")[[1]]))
  g.w = set.vertex.attribute(g.w, attr, V(g.w)[nw[1,]], nw[2,])
  g.w
}

kwalksHeatmap = function(topnodes, kwalks.tbls, fn = "heatmap.kwalks.pdf", colorByDEG = F, ...) {
  mat = NULL
  for(d in kwalks.treatments) {
    tbl = kwalks.tbls[[d]]
    tbl = tbl[topnodes,]
    
    if(colorByDEG) {
      ## Set sign to DEG direction
      rg = rgraphs[[grep(d, names(rgraphs))]]
      fcattr = grep(paste0("logFC.+", d), list.vertex.attributes(rg), value = T)
      fc = rep(NA, length(topnodes))
      fc[topnodes %in% V(rg)$name] = get.vertex.attribute(rg, fcattr, topnodes[topnodes %in% V(rg)$name])
      tbl = sign(fc) * tbl
    }
    
    if(is.null(mat)) {
      mat = tbl
    } else {
      spacer = rep(0, nrow(tbl))
      mat = cbind(mat, cbind(spacer, spacer))
      mat = cbind(mat, tbl)
      colnames(mat)[colnames(mat) == "spacer"] = ""
    }
  }
  mode(mat) = "numeric"
  rownames(mat) = nodeInfo$symbols[sub("L:", "", rownames(mat))]
  
  annotations = node2mod[topnodes]
  annotations[annotations == ""] = "none"
  annotations = as.data.frame(as.factor(annotations))
  
  annotations$drugInteraction = rep(" ", nrow(annotations))
  for(d in kwalks.treatments) {
    nbs = intersect(drugNbs[[d]], rownames(annotations))
    annotations[nbs, 2] = d
  }
  
  colnames(annotations) = c("module", "drugInteraction")
  rownames(annotations) = rownames(mat)
  
  mod.colors = as.character(unique(annotations[,1]))
  names(mod.colors) = mod.colors
  if("none" %in% names(mod.colors)) mod.colors["none"] = "white"
  dnb.colors = palette()[2:(length(unique(annotations[,2])) + 1)]
  names(dnb.colors) = as.character(unique(annotations[,2]))
  dnb.colors[" "] = "white"
  annot.colors = list(module = mod.colors, drugInteraction = dnb.colors)
  
  library(pheatmap)
  
  nc = 256
  if(colorByDEG) {
    min = -2
    max = 2
    colors = colorRampPalette(c('blue', 'white', 'red'))(nc)
  } else {
    min = 0
    max = 2
    colors = colorRampPalette(c('white', 'blue'))(nc)
  }
  breaks = seq(min, max, (max - min)/(nc-1))
  
  mat.sat = saturateValues(mat, min, max)
  mat.sat[is.na(mat.sat)] = 0
  pheatmap(t(mat.sat),
           color = colors, scale = "none", breaks = breaks,
           legend = T, show_rownames = T, 
           annotation = annotations, annotation_colors = annot.colors,
           cellwidth = 5, cellheight = 5,
           border_color = NA, fontsize = 4.5,
           filename = paste(myOutPath, fn, sep=""),
           cluster_rows = !rankSort,
           ...
  )
}

extractTop = function(tbl, nr = 5) {
  top = list()
  for(cn in colnames(tbl)) {
    mod = c(kwalks.treatments, names(which(node2mod == cn, useNames = T)))
    col = tbl[, cn]
    top[[cn]] = names(col)[order(col, decreasing=T)[1:nr]]
  }
  top
}

kwalksHeatmap.bygroup = function(
  topnodes, kwalks.tbls, fn = "heatmap.kwalks.pdf", colorByDEG = F, 
  firstGroup = c(), distFun = dist, firstGroupCol = "",
  clustMethod = "complete", rankSort = F, ...
) {
  mat = NULL
  
  for(d in names(kwalks.tbls)) {
    tbl = kwalks.tbls[[d]]
    tbl = tbl[topnodes,]
    
    if(colorByDEG) {
      ## Set sign to DEG direction
      rg = rgraphs[[grep(d, names(rgraphs))]]
      fcattr = grep(paste0("logFC.+", d), list.vertex.attributes(rg), value = T)
      fc = rep(NA, length(topnodes))
      fc[topnodes %in% V(rg)$name] = get.vertex.attribute(rg, fcattr, topnodes[topnodes %in% V(rg)$name])
      tbl = sign(fc) * tbl
    }
    
    colnames(tbl) = paste(colnames(tbl), d)
    if(is.null(mat)) {
      mat = tbl
    } else {
      spacer = rep(0, nrow(tbl))
      mat = cbind(cbind(spacer, spacer), mat)
      mat = cbind(tbl, mat)
      colnames(mat)[colnames(mat) == "spacer"] = ""
    }
  }
  mode(mat) = "numeric"
  rownames(mat) = nodeInfo$symbols[sub("L:", "", rownames(mat))]
  firstGroup = nodeInfo$symbols[sub("L:", "", firstGroup)]
  
  mat[is.na(mat)] = 0
  
  docluster = function(m, isFirstGroup = F) {
    if(rankSort) {
      mo = m
      if(isFirstGroup) mo = m[, grep(firstGroupCol, colnames(m))]
      roworder = order(apply(abs(mo), 1, max), decreasing = F)
      m = m[roworder,]
    } else {
      m[hclust(distFun(m), method=clustMethod)$order,]
    }
    m
  }
  
  ## Group and cluster rows
  ## Cluster rest
  mat.add = mat[setdiff(rownames(mat), firstGroup),]
  mat.add = docluster(mat.add)
  ## Cluster group and merge with rest
  mat.grp = mat[firstGroup,]
  if(length(firstGroup) > 0) {
    mat.grp = docluster(mat.grp, T)
    spacer = rep(0, ncol(mat))
    mat = rbind(mat.add, ` ` = spacer, mat.grp)
  }
  
  annotations = node2mod[topnodes]
  annotations[annotations == ""] = "none"
  annotations = as.data.frame(as.factor(annotations))
  
  annotations$drugInteraction = rep(" ", nrow(annotations))
  for(d in kwalks.treatments) {
    nbs = intersect(drugNbs[[d]], rownames(annotations))
    annotations[nbs, 2] = d
  }
  
  colnames(annotations) = c("module", "drugInteraction")
  rownames(annotations) = nodeInfo$symbols[sub("L:", "", rownames(annotations))]
  annotations = rbind(annotations, ` ` = c("none", " "))
  
  annotations = annotations[rownames(mat),]
  mod.colors = as.character(unique(annotations[,1]))
  names(mod.colors) = mod.colors
  if("none" %in% names(mod.colors)) mod.colors["none"] = "white"
  dnb.colors = palette()[2:(length(unique(annotations[,2])) + 1)]
  names(dnb.colors) = as.character(unique(annotations[,2]))
  dnb.colors[" "] = "white"
  annot.colors = list(module = mod.colors, drugInteraction = dnb.colors)
  
  library(pheatmap)
  
  nc = 256
  max = 3
  if(colorByDEG) {
    min = -max
    colors = colorRampPalette(c('blue', 'white', 'red'))(nc)
  } else {
    min = 0
    colors = colorRampPalette(c('white', 'blue'))(nc)
  }
  breaks = seq(min, max, (max - min)/(nc-1))
  
  mat.sat = saturateValues(mat, min, max)
  mat.sat[is.na(mat.sat)] = 0
  pheatmap(t(mat.sat),
           color = colors, scale = "none", breaks = breaks,
           legend = T, show_rownames = T, 
           annotation = annotations, annotation_colors = annot.colors,
           cellwidth = 5, cellheight = 5,
           border_color = NA, fontsize = 4.5,
           filename = paste(myOutPath, fn, sep=""),
           cluster_cols = F, cluster_rows = F,
           ...
  )
}