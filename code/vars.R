###################################################
## System dependent variables for 10OAD analysis ##
###################################################
projPath = paste0(getwd(), '/')
outPath = paste0(projPath, 'output/')
dataPath = paste0(projPath, 'data/')
codePath = paste0(projPath, 'code/')

########################
## Required libraries ##
########################
require(org.Mm.eg.db)
require(WGCNA)
require(digest)
require(devtools)
require(limma)
require(org.Hs.eg.db)
library(GOstats)

######################
## Common functions ##
######################
source(paste0(codePath, 'functions.R'))

source('code/vars.limmadata.R')
source('code/vars.mousedata.R')

inclPars = c(
  "body.weight",                         
  "liver.total",                            
  "heart.total",                            
  "visceral.WAT",                           
  "gonadal.WAT",                            
  "subcutaneous.WAT",
  "ratio.visc.sub.WAT",                     
  "kidneys",                                
  "liver.triglycerides",                    
  "venes.total.atherosclerotic.lesion.area",
  "urine.glucose",                          
  "cholesterol",                            
  "triglycerides",                          
  "glucose",                                
  "insulin",                               
  "QUICKI"
)                                 

inclPars.labels = gsub("\\.", " ", inclPars)
inclPars.labels = gsub("^(.)", "\\U\\1", inclPars.labels, perl = T)
inclPars.labels = gsub("total", "weight", inclPars.labels)
inclPars.labels[4] = "Omental fat weight"
inclPars.labels[5] = "Epididymal fat weight"
inclPars.labels[6] = "Inguinal fat weight"
inclPars.labels[9] = "Intrahepatic triglycerides"
inclPars.labels[10] = "Atherosclerotic lesion area"

kwalksModules = c("red", "black", "yellow")

p.max = 0.05
treatments = c("16wkHighFat", "LifeStyle", "Fenofibrate", "T0901317")
treatments.p = pvalCols.unadj[unlist(lapply(paste(treatments, "\\.", sep=""), grep, pvalCols.unadj))]
treatments.fdr = pvalCols[unlist(lapply(paste(treatments, "\\.", sep=""), grep, pvalCols))]
treatments.logfc = logfcCols[unlist(lapply(paste(treatments, "\\.", sep=""), grep, logfcCols))]