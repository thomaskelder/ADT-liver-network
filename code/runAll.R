## Uncomment these lines to install the required packages
#install.packages(c("devtools", "igraph", "WGCNA", "pheatmap", "RCurl", "ggplot2", "digest"))
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("limma", "org.Mm.eg.db", "org.Hs.eg.db", "GOstats", "impute", "biomaRt"))

source("code/vars.R")
dir.create(paste0(outPath, "wgcna"))
dir.create(paste0(outPath, "cytoscape"))
dir.create(paste0(outPath, "kwalks"))
dir.create(paste0(outPath, "networks"))
dir.create(paste0(outPath, "roc"))

## WGCNA to find modules
source("code/wgcna.R")

## Extend modules with additional network information
source("code/stitch.graphs.R")
source("code/diseaseAssociations.R")
source("code/intervention.graphs.R")

## Apply kWalks to extract subnetworks that connect intervention
## targets with modules
source("code/intervention.kwalks.R")
