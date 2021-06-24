# SpiderLearner

This code implements an ensemble model for GGM estimation.

A quick start guide is available at https://katehoffshutta.github.io/web/SpiderLearner/SpiderLearnerQuickStart.html. The code to generate the html is in the file SpiderLearnerQuickStart.Rmd.

It should take 5-10 minutes to run. If you would like to see a faster example, you can change from 10-fold to 3-fold cross-validation by changing the line:

```slResults = s$runSpiderLearner(lateStageSmall, K = 10, standardize=T, nCores = 1)```

to

```slResults = s$runSpiderLearner(lateStageSmall, K = 3, standardize=T, nCores = 1)```

or use multiple cores if available. You can also run the demo in the quick start guide by downloading the repository and sourcing the R script SpiderLearnerDemo.R, which is essentially the same as SpiderLearnerQuickStart.Rmd but does not require R Markdown.  

Some R packages you will need to install: affy, BiocManager, curatedOvarianData, doParallel, foreach, hglasso, huge, qgraph, R6, Rsolnp

Feedback encouraged: kshutta@umass.edu
