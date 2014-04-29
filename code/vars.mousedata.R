################################
## Load 10OAD data per mouse  ##
################################
cacheFile = paste0(dataPath, "mouse.data.RData")
if(file.exists(cacheFile)) {
  load(cacheFile) 
} else {
  ## Gene expression (liver)
  expr.liver.file = paste0(dataPath, "liver TX data 10OAD LDLR-KO mice.txt")
  expr.liver.data = read.delim(expr.liver.file, header=F, stringsAsFactor = F)
  expr.liver.data = expr.liver.data[,c(1:3, 5:ncol(expr.liver.data))]
  expr.liver.data.info = expr.liver.data[1:7,]
  expr.liver.data = expr.liver.data[8:nrow(expr.liver.data),]
  rownames(expr.liver.data) = NULL
  
  colnames(expr.liver.data) = as.character(
    c(expr.liver.data.info[7,1:3], expr.liver.data.info[1,4:ncol(expr.liver.data.info)])
  )
  expr.liver.data.genes = expr.liver.data[,1:3]
  ## update the entrez id with more recent annotation from limma analysis
  data.entrez = data[, "entrezID"]
  names(data.entrez) = data[, "probeID"]
  expr.liver.data.genes[,3] = data.entrez[expr.liver.data.genes[,1]]
  
  expr.liver.data = expr.liver.data[, 4:ncol(expr.liver.data)]
  expr.liver.data = as.matrix(expr.liver.data)
  mode(expr.liver.data) = "numeric"
  
  rownames(expr.liver.data) = expr.liver.data.genes[,1]
  
  ## Physiological paramters
  phys.file = paste0(dataPath, "physiological data 10OAD LDLR-KO mice T=16_v5-DAP_20111124.tk.txt")
  ypar.data = read.delim(phys.file, stringsAsFactor = F, skip = 4)
  
  mouse.nrs = ypar.data[,"mouse"]
  mouse.trt = ypar.data[,"treatment"]
  names(mouse.trt) = mouse.nrs
  
  ypar.data = ypar.data[,c(13,15:24,26,28,30,32:42,44,45,48)]
  rownames(ypar.data) = mouse.nrs
  
  ypar.data = t(ypar.data)
  mode(ypar.data) = 'numeric'
  save(mouse.nrs, mouse.trt, ypar.data, expr.liver.data, expr.liver.data.genes, file = cacheFile)
}