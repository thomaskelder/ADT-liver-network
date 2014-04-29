############################
## Load 10OAD limma data  ##
############################
cacheFile = paste0(dataPath, "limma.data.RData")
if(file.exists(cacheFile)) {
  load(cacheFile)  
} else {
  ## Read 10OAD data
  dataFile = paste0(dataPath, "/limma.10OAD.combined.txt")
  data = read.delim(dataFile, sep="\t", stringsAsFactors = F)
  
  ## Select columns that contain p-value
  pvalCols = grep("adj\\.P\\.Val_", colnames(data), value=T)
  pvalCols.unadj = grep("P\\.Value_", colnames(data), value=T)
  tstatCols = grep("t_", colnames(data), value=T)
  logfcCols = grep("logFC_", colnames(data), value=T)
  xrefCol = "xref"
  
  treatmentLabels = sub("adj\\.P\\.Val_LDLR_", "", pvalCols)
  treatmentLabels = sub(".LDLR_Chow", "", treatmentLabels)
  treatmentLabels = sub(".LDLR_16wkHighFat", "", treatmentLabels)
  
  data.xref = as.matrix(data[, 9:ncol(data)])
  mode(data.xref) = 'numeric'
  data.xref = avereps(data.xref, ID=data[,'xref'])
  
  save(
    data, data.xref,
    pvalCols, pvalCols.unadj,
    tstatCols, logfcCols,
    xrefCol, treatmentLabels,
    file = cacheFile
  )
}