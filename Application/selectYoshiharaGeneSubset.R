library(corrplot)
library(affy)
library(ensembleGGM)
library(curatedOvarianData)

data("GSE32062.GPL6480_eset")
data("E.MTAB.386_eset")
data("GSE13876_eset")
data("GSE14764_eset")
data("GSE17260_eset")
data("GSE18520_eset")
data("GSE19829.GPL570_eset")
data("GSE19829.GPL8300_eset")
data("GSE26712_eset")
data("GSE30009_eset")
data("GSE30161_eset")
data("GSE32063_eset")
data("GSE9891_eset")
data("PMID17290060_eset")
data("PMID19318476_eset")
data("TCGA_eset")

identifiers = c("GSE32062.GPL6480",
                "E.MTAB.386",
                "GSE13876",
                "GSE14764",
                "GSE17260",
                "GSE18520",
                "GSE19829.GPL570",
                "GSE19829.GPL8300",
                "GSE26712",
                "GSE30009",
                "GSE30161",
                "GSE32063",
                "GSE9891",
                "PMID17290060",
                "PMID19318476",
                "TCGA")

trainAndValidData = list(exprs(GSE32062.GPL6480_eset),
                 exprs(E.MTAB.386_eset),
                 exprs(GSE13876_eset),
                 exprs(GSE14764_eset),
                 exprs(GSE17260_eset),
                 exprs(GSE18520_eset),
                 exprs(GSE19829.GPL570_eset),
                 exprs(GSE19829.GPL8300_eset),
                 exprs(GSE26712_eset),
                 exprs(GSE30009_eset),
                 exprs(GSE30161_eset),
                 exprs(GSE32063_eset),
                 exprs(GSE9891_eset),
                 exprs(PMID17290060_eset),
                 exprs(PMID19318476_eset),
                 exprs(TCGA_eset))

yoshi = read.table("YoshiharaGeneSet.tsv",sep="\t")
# Figure out which of the 116 genes are in all of the datasets
geneMatch = matrix(rep(NA,124*16),ncol=124)
for(i in 1:length(identifiers))
{
  measuredGenes = colnames(t(trainAndValidData[[i]]))
  geneMatch[i,] = ifelse(yoshi[,1] %in% measuredGenes,1,0)
}

corrplot(geneMatch)
# Omit validation datasets 1,7,9 (indices 2,8,10), which have high missingness
corrplot(geneMatch[-c(2,8,10),])
# Omit any remaining missing genes (~15 genes)
missingGenes = apply(geneMatch[-c(2,8,10),],2,function(x){sum(x==0)>0})
consensusGenes = data.frame(yoshi[which(missingGenes==FALSE),1])
write.table(consensusGenes,file="yoshiharaLimitedGeneSet.tsv",sep="\t",quote=F,row.names=F,col.names = F)
corrplot(geneMatch[-10,(missingGenes==FALSE)])
