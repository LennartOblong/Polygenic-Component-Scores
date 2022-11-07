

###Code for a linear mixed effect model, which will identify the need for a linear regression or a logistic regression automatically ###

require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)
require(lattice)

data <- read.csv("projects/0/einf2700/oblongl/logistic_regression/regression_file")

data <- within(data, {
   <- factor(FID, levels = 0:1, labels = c("no", "yes"))
  DID <- factor(DID)
  HID <- factor(HID)
  CancerStage <- factor(CancerStage)
})


glmer(PRS ~ SUBAFF_EVER*Sex + Site + Batch, random=~1 | FID, data = data, na.action = na.omit)
