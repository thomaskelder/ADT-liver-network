##############################################################
## Global effects of treatment                              ##
## creates lists that store which genes have been reverted, ##
## amplified or are an additional effect.                   ##
##############################################################
source("code/vars.R")

getAmplifiedReversed = function(pval.sig, pval.type, pval.sig.hf = pval.sig) {
	sig.hf = data.xref[, paste(pval.type, "_LDLR_16wkHighFat.LDLR_Chow", sep="")] < pval.sig.hf
	up.hf = data.xref[, "t_LDLR_16wkHighFat.LDLR_Chow"] > 0

	reverted = sapply(sub("t_", "", tstatCols[3:length(tstatCols)]), function(trt) {
		sig = data.xref[,paste(pval.type, "_", trt, sep="")] < pval.sig
		up = data.xref[,paste("t_", trt, sep="")] > 0
		sig.hf & sig & (up.hf == !up)
	})

	amplified = sapply(sub("t_", "", tstatCols[3:length(tstatCols)]), function(trt) {
		sig = data.xref[,paste(pval.type, "_", trt, sep="")] < pval.sig
		up = data.xref[,paste("t_", trt, sep="")] > 0
		sig.hf & sig & (up.hf == up)
	})

	additional = sapply(sub("t_", "", tstatCols[3:length(tstatCols)]), function(trt) {
		sig = data.xref[,paste(pval.type, "_", trt, sep="")] < pval.sig
		!sig.hf & sig
	})
	
	list(deg.hf = sig.hf, rev = reverted, amp = amplified, add = additional)
}