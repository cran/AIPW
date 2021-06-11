## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 6,
  comment = "#>"#,
  # cache=TRUE
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("remotes")
#  remotes::install_github("yqzhong7/AIPW")

## ---- eval = FALSE------------------------------------------------------------
#  #SuperLearner
#  install.packages("SuperLearner")
#  #sl3
#  remotes::install_github("tlverse/sl3")
#  install.packages("Rsolnp")

## ----example data-------------------------------------------------------------
library(AIPW)
library(SuperLearner)
library(ggplot2)
set.seed(123)
data("eager_sim_obs")
cov = c("eligibility","loss_num","age", "time_try_pregnant","BMI","meanAP")

## ----one_line-----------------------------------------------------------------
AIPW_SL <- AIPW$new(Y= eager_sim_obs$sim_Y,
                    A= eager_sim_obs$sim_A,
                    W= subset(eager_sim_obs,select=cov), 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 10,
                    verbose=FALSE)$
  fit()$
  #Default truncation is set to 0.025; using 0.25 here is for illustrative purposes and not recommended
  summary(g.bound = c(0.25,0.75))$
  plot.p_score()$
  plot.ip_weights()

## ----SuperLearner, message=FALSE,eval=F---------------------------------------
#  library(AIPW)
#  library(SuperLearner)
#  
#  #SuperLearner libraries for outcome (Q) and exposure models (g)
#  sl.lib <- c("SL.mean","SL.glm")
#  
#  #construct an aipw object for later estimations
#  AIPW_SL <- AIPW$new(Y= eager_sim_obs$sim_Y,
#                      A= eager_sim_obs$sim_A,
#                      W= subset(eager_sim_obs,select=cov),
#                      Q.SL.library = sl.lib,
#                      g.SL.library = sl.lib,
#                      k_split = 10,
#                      verbose=FALSE)

## -----------------------------------------------------------------------------
#fit the AIPW_SL object
AIPW_SL$fit()
# or you can use stratified_fit
# AIPW_SL$stratified_fit()

## -----------------------------------------------------------------------------
#estimate the average causal effects from the fitted AIPW_SL object 
AIPW_SL$summary(g.bound = 0.25) #propensity score truncation 

## ----ps_trunc-----------------------------------------------------------------
library(ggplot2)
AIPW_SL$plot.p_score()
AIPW_SL$plot.ip_weights()

## -----------------------------------------------------------------------------
suppressWarnings({
  AIPW_SL$stratified_fit()$summary()
})

## ----parallel, eval=FALSE-----------------------------------------------------
#  # install.packages("future.apply")
#  library(future.apply)
#  plan(multiprocess, workers=2, gc=T)
#  set.seed(888)
#  AIPW_SL <- AIPW$new(Y= eager_sim_obs$sim_Y,
#                      A= eager_sim_obs$sim_A,
#                      W= subset(eager_sim_obs,select=cov),
#                      Q.SL.library = sl3.lib,
#                      g.SL.library = sl3.lib,
#                      k_split = 10,
#                      verbose=FALSE)$fit()$summary()

## ----tmle, eval=F-------------------------------------------------------------
#  # install.packages("tmle")
#  library(tmle)
#  library(SuperLearner)
#  
#  tmle_fit <- tmle(Y=eager_sim_obs$sim_Y,
#                   A=eager_sim_obs$sim_A,
#                   W=eager_sim_obs[,-1:-2],
#                   Q.SL.library=c("SL.mean","SL.glm"),
#                   g.SL.library=c("SL.mean","SL.glm"),
#                   family="binomial",
#                   cvQinit = TRUE)
#  
#  cat("\nEstimates from TMLE\n")
#  unlist(tmle_fit$estimates$ATE)
#  unlist(tmle_fit$estimates$RR)
#  unlist(tmle_fit$estimates$OR)
#  
#  cat("\nEstimates from AIPW\n")
#  a_tmle <- AIPW_tmle$
#    new(A=eager_sim_obs$sim_A,Y=eager_sim_obs$sim_Y,tmle_fit = tmle_fit,verbose = TRUE)$
#    summary(g.bound=0.025)

