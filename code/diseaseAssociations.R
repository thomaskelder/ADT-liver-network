################################################
## Parse human genetic disease associations   ##
## and write as Cytoscape attributes          ##
################################################
source('code/vars.R')
inPath = paste(dataPath, "disease_association/", sep="")

select = list()
select$`diabetes type 2` = c()
select$`insulin resistance` = c()
select$`metabolic syndrome` = c()

## Read disease associations from geneticassociationdb
gadb = read.delim(paste(inPath, "/geneticassociationdb/all.txt", sep=""), sep="\t", skip=2)

# Remove negative associations
gadb = gadb[gadb[,"Association.Y.N."] != "N", ]

# Find relevant disease terms
pheno = gadb[,"Broad.Phenotype"]
select$`diabetes type 2` = gadb[grep("diabetes, type 2", pheno, ignore.case=T), "Locus.Number"]
select$`insulin resistance` = gadb[grep("insulin resistance", pheno, ignore.case=T), "Locus.Number"]
select$`metabolic syndrome` = gadb[grep("metabolic syndrome", pheno, ignore.case=T), "Locus.Number"]
select$obesity = gadb[grep("obesity", pheno, ignore.case=T), "Locus.Number"]

## Read associations from HUGE
symbol2eg = as.list(org.Hs.egSYMBOL2EG)

getHugeIds = function(file) {
  data = read.delim(file, sep="\t", skip=5)
  symbols = as.character(data[,1])
  unlist(symbol2eg[symbols])
}
select[['insulin resistance']] = 
  union(select[['insulin resistance']], getHugeIds(paste(inPath, "Insulin Resistance_phenopedia_11-10-2011.txt", sep="")))
select[['metabolic syndrome']] = 
  union(select[['metabolic syndrome']], getHugeIds(paste(inPath, "Metabolic Syndrome X_phenopedia_11-10-2011.txt", sep="")))
select[['obesity']] = 
  union(select[['obesity']], getHugeIds(paste(inPath, "Obesity_phenopedia_11-10-2011.txt", sep="")))

names(select) = paste('genetic ', names(select))

## Load PubMed/textmining results (http://fable.chop.edu/)
minArticles = 5
fable = read.delim(paste0(inPath, "fable - type 2 diabetes mellitus or metabolic syndrome.txt"), as.is=T)
select$`fable t2dm or metsyn` = symbol2eg[fable$Symbol[fable$Articles >= minArticles]]

fable = read.delim(paste0(inPath, "fable - (type 2 diabetes mellitus or metabolic syndrome) and liver.txt"), as.is=T)
select$`fable t2dm or metsyn liver` = symbol2eg[fable$Symbol[fable$Articles >= minArticles]]

fable = read.delim(paste0(inPath, "fable - dyslipidemia.txt"), as.is=T)
select$`fable dyslipidemia` = symbol2eg[fable$Symbol[fable$Articles >= minArticles]]

fable = read.delim(paste0(inPath, "fable - type 2 diabetes mellitus and liver.txt"), as.is=T)
select$`fable t2dm and liver` = symbol2eg[fable$Symbol[fable$Articles >= minArticles]]

associated.ids = unique(unlist(select))

associations = sapply(associated.ids, function(id) {
  pts = c()
  for(p in names(select)) {
    if(id %in% select[[p]]) pts = c(pts, p)
  }
  pts = sort(pts)
  paste(pts, collapse = ', ')
})
names(associations) = associated.ids

## Map to mouse homologs (from http://mikedewar.wordpress.com/2010/05/14/generating-homologues-using-biomart/)
library(biomaRt)
gen_hs2mm <- function(){
    ensembl_hs <- useMart(
        "ensembl",
        dataset = "hsapiens_gene_ensembl"
    )
    hs2mm_filters <- c(
        "with_homolog_mmus"
    )
    hs2mm_homo_atts <- c(
        "ensembl_gene_id",
        "mmusculus_homolog_ensembl_gene"
    )
    # the names in these lists are arbitrary
    hs2mm_value = TRUE

    # get the human genes and mouse orthologues
    hs2mm_homo <- getBM(
        attributes = hs2mm_homo_atts,
        filters = hs2mm_filters,
        value = hs2mm_value,
        mart = ensembl_hs
    )
    hs2mm_homo
}
hs2mm = gen_hs2mm()

ids.mouse = c()
associations.mouse = c()
eg2ens = as.list(org.Hs.egENSEMBL)
ens2eg = as.list(org.Mm.egENSEMBL2EG)

for(id in associated.ids) {
  # Get ensembl id
  ens = as.character(eg2ens[[id]])
  ens.mm = hs2mm[hs2mm[,1] %in% ens, 2]
  eg.mm = unlist(ens2eg[ens.mm])
  
  ids.mouse = c(ids.mouse, eg.mm)
  associations.mouse = c(associations.mouse, rep(associations[id], length(eg.mm)))
}
names(associations.mouse) = ids.mouse
associations.mouse = associations.mouse[unique(names(associations.mouse))]

## Write as cytoscape attribute file
attr = "human_genetic_association"
attr.ids = paste("L:", names(associations.mouse), sep="")
attr.values = c(attr, paste(attr.ids, associations.mouse, sep=" = "))
write.table(attr.values, paste(outPath, "human_genetic_association.txt", sep=""), col.names=F, row.names=F, quote=F)

save(associations.mouse, file = paste0(outPath, "associations.RData"))