# This file is for getting the data out for
# the table about the discovery and validation datasets

library(curatedOvarianData)
library(survMisc)

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
data("GSE32062.GPL6480_eset")
data("GSE32063_eset")
data("GSE9891_eset")
data("PMID17290060_eset")
data("PMID19318476_eset")
data("TCGA_eset")

phenoDataValid = list(pData(E.MTAB.386_eset),
                 pData(GSE13876_eset),
                 pData(GSE14764_eset),
                 pData(GSE17260_eset),
                 pData(GSE18520_eset),
                 pData(GSE19829.GPL570_eset),
                 pData(GSE19829.GPL8300_eset),
                 pData(GSE26712_eset),
                 pData(GSE30009_eset),
                 pData(GSE30161_eset),
                 pData(GSE32062.GPL6480_eset),
                 pData(GSE32063_eset),
                 pData(GSE9891_eset),
                 pData(PMID17290060_eset),
                 pData(PMID19318476_eset),
                 pData(TCGA_eset))

platforms = c(annotation(E.MTAB.386_eset),
              annotation(GSE13876_eset),
              annotation(GSE14764_eset),
              annotation(GSE17260_eset),
              annotation(GSE18520_eset),
              annotation(GSE19829.GPL570_eset),
              annotation(GSE19829.GPL8300_eset),
              annotation(GSE26712_eset),
              annotation(GSE30009_eset),
              annotation(GSE30161_eset),
              annotation(GSE32062.GPL6480_eset),
              annotation(GSE32063_eset),
              annotation(GSE9891_eset),
              annotation(PMID17290060_eset),
              annotation(PMID19318476_eset),
              annotation(TCGA_eset))
              

identifiers = c("E.MTAB.386",
                "GSE13876",
                "GSE14764",
                "GSE17260",
                "GSE18520",
                "GSE19829.GPL570",
                "GSE19829.GPL8300",
                "GSE26712",
                "GSE30009",
                "GSE30161",
                "GSE32062.GPL6480_eset",
                "GSE32063",
                "GSE9891",
                "PMID17290060",
                "PMID19318476",
                "TCGA")

sampleSize = sapply(phenoDataValid, nrow)
meanAge = sapply(phenoDataValid, function(x){mean(x$age_at_initial_pathologic_diagnosis, na.rm=T)})
sdAge = sapply(phenoDataValid, function(x){sd(x$age_at_initial_pathologic_diagnosis, na.rm=T)})
missingAge = sapply(phenoDataValid, function(x){sum(is.na((x$age_at_initial_pathologic_diagnosis)))})
       
tumorStage =sapply(phenoDataValid, function(x){mean(x$tumorstage < 4,na.rm=T)})
missingTumorStage = sapply(phenoDataValid, function(x){sum(is.na(x$tumorstage))})

summaryStageEarly = sapply(phenoDataValid, function(x){sum(x$summarystage == "early",na.rm=T)})
summaryStageLate = sapply(phenoDataValid, function(x){sum(x$summarystage == "late",na.rm=T)})
missingSummaryStage = sapply(phenoDataValid, function(x){sum(is.na(x$summarystage))})

summaryGradeLow = sapply(phenoDataValid, function(x){sum(x$summarygrade == "low",na.rm=T)})
summaryGradeHigh = sapply(phenoDataValid, function(x){sum(x$summarygrade == "high",na.rm=T)})
missingSummaryGrade = sapply(phenoDataValid, function(x){sum(is.na(x$summarygrade))})

medianSurvival = sapply(phenoDataValid,function(x){
  mod = survfit(Surv(days_to_death, vital_status == "living") ~ 1, data = x);
  return(summary(mod)$table[7]) })
missingMedianSurvival = sapply(phenoDataValid,function(x){sum(is.na(x$days_to_death))})

censoringRate = sapply(phenoDataValid,function(x){sum(x$vital_status == "living")})/sampleSize

outDF = data.frame("dataset" = identifiers,
                   "platform" = platforms,
                   "sampleSize" = sampleSize,
                   "meanAge" = round(meanAge,2),
                   "sdAge" = round(sdAge,2),
                   "missingAge" = missingAge,
                   "propTumorStageLT4" = round(tumorStage,2),
                   "missingTumorStage" = missingTumorStage,
                   "propSummaryStageEarly" = round(summaryStageEarly/sampleSize,2),
                   "propSummaryStageLate" = round(summaryStageLate/sampleSize,2),
                   "missingSummaryStage" = missingSummaryStage,
                   "propSummaryGradeLow" = round(summaryGradeLow/sampleSize,2),
                   "propSummaryGradeHigh" = round(summaryGradeHigh/sampleSize,2),
                   "missingSummaryGrade" = missingSummaryGrade,
                   "medianSurvival" = medianSurvival,
                   "missingMedianSurvival" = missingMedianSurvival)

write.csv(outDF,"Tables/cohortCharacteristics.csv",row.names=F,quote=F)
