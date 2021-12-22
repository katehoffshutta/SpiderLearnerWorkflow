library(affy)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("curatedOvarianData")
library(curatedOvarianData)
library(ensembleGGM)

standardize = function(x){return((x-mean(x))/sd(x))}

data(GSE32062.GPL6480_eset)
lateStage = exprs(GSE32062.GPL6480_eset)

# Extract genes included in the HPO ovarian carcinoma pathway
# See https://hpo.jax.org/app/browse/term/HP:0025318
# Reference: Kohler et al. (2021) The Human Phenotype Ontology in 2021, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D1207–D1217, https://doi.org/10.1093/nar/gkaa1043

lateStageSmall = lateStage[c(1680,1681,3027,4564,8930,12243,12245,13694,13695,13701,13979,16082,16875,17980),]
lateStageSmall = t(lateStageSmall)
names(lateStageSmall) = colnames(lateStageSmall)
lateStageSmall = apply(lateStageSmall,2,standardize)

s1 = SpiderLearner$new()
s2 = SpiderLearner$new()
s3 = SpiderLearner$new()
s4 = SpiderLearner$new()

apple = HugeEBICCandidate$new(gamma = 0)
banana = HugeEBICCandidate$new(gamma = 0.5)
clementine = HugeRICCandidate$new()
date = HGlassoCandidate$new()
elderberry = MLECandidate$new() # don't add MLE forHD case
fraise = HugeStARSCandidate$new(thres = 0.05)
grape = HugeStARSCandidate$new(thres = 0.1)
honeydew = QGraphEBICCandidate$new(gamma = 0)
icewine = QGraphEBICCandidate$new(gamma = 0.5)

candidates_1 = list(apple, 
                  banana, 
                  clementine, 
                  date, 
                  elderberry,
                  fraise,
                  grape,
                  honeydew,
                  icewine)
candidates_2 = list(apple,elderberry,honeydew)
candidates_3 = list(banana, clementine, date, fraise, grape, icewine)
candidates_4 = list(banana, fraise, grape, icewine)

for(candidate in candidates_1)
{
  s1$addCandidate(candidate)
}

for(candidate in candidates_2)
{
  s2$addCandidate(candidate)
}

for(candidate in candidates_3)
{
  s3$addCandidate(candidate)
}

for(candidate in candidates_4)
{
  s4$addCandidate(candidate)
}

set.seed(1222)

slResults1 = s1$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 10)
slResults2 = s2$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 10)
slResults3 = s3$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 10)
slResults4 = s4$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 10)

slResults = list(slResults1, slResults2, slResults3, slResults4)

save(slResults,file="ovarianLibraryExample.rda")


