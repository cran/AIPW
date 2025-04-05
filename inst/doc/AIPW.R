## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 6,
  comment = "#>"#,
  # cache=TRUE
)

## ---- eval = FALSE------------------------------------------------------------
#  # CRAN version
#  install.packages("AIPW")
#  # github version
#  # install.packages("remotes")
#  # remotes::install_github("yqzhong7/AIPW")

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
  #Default truncation
  summary(g.bound = c(0.025,0.975))$
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

## ----sl3, eval=F--------------------------------------------------------------
#  library(AIPW)
#  library(sl3)
#  
#  ##construct sl3 learners for outcome (Q) and exposure models (g)
#  lrnr_glm <- Lrnr_glm$new()
#  lrnr_mean <- Lrnr_mean$new()
#  #stacking two learner (this will yield estimates for each learner)
#  stacklearner <- Stack$new(lrnr_glm, lrnr_mean)
#  #metalearner is required to combine the estimates from stacklearner
#  metalearner <- Lrnr_nnls$new()
#  sl3.lib <- Lrnr_sl$new(learners = stacklearner,
#                         metalearner = metalearner)
#  
#  #construct an aipw object for later estimations
#  AIPW_sl3 <- AIPW$new(Y= eager_sim_obs$sim_Y,
#                       A= eager_sim_obs$sim_A,
#                       W= subset(eager_sim_obs,select=cov),
#                       Q.SL.library = sl3.lib,
#                       g.SL.library = sl3.lib,
#                       k_split = 10,
#                       verbose=FALSE)

## -----------------------------------------------------------------------------
#fit the AIPW_SL object
AIPW_SL$fit()
# or you can use stratified_fit
# AIPW_SL$stratified_fit()

## -----------------------------------------------------------------------------
#estimate the average causal effects from the fitted AIPW_SL object 
AIPW_SL$summary(g.bound = 0.025) #propensity score truncation 

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

## ----tmle3, eval=F------------------------------------------------------------
#  # remotes::install_github("tlverse/tmle3")
#  library(tmle3,quietly = TRUE)
#  library(sl3,quietly = TRUE)
#  
#  node_list <- list(A = "sim_A",Y = "sim_Y",W = colnames(eager_sim_obs)[-1:-2])
#  or_spec <- tmle_OR(baseline_level = "0",contrast_level = "1")
#  tmle_task <- or_spec$make_tmle_task(eager_sim_obs,node_list)
#  lrnr_glm <- make_learner(Lrnr_glm)
#  lrnr_mean <- make_learner(Lrnr_mean)
#  sl <- Lrnr_sl$new(learners = list(lrnr_glm,lrnr_mean))
#  learner_list <- list(A = sl, Y = sl)
#  tmle3_fit <- tmle3(or_spec, data=eager_sim_obs, node_list, learner_list)
#  
#  cat("\nEstimates from TMLE\n")
#  tmle3_fit$summary
#  
#  # parse tmle3_fit into AIPW_tmle class
#  cat("\nEstimates from AIPW\n")
#  a_tmle3<- AIPW_tmle$
#    new(A=eager_sim_obs$sim_A,Y=eager_sim_obs$sim_Y,tmle_fit = tmle3_fit,verbose = TRUE)$
#    summary(g.bound=0)

