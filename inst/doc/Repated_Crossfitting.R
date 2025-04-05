## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6
)

## ----one_line-----------------------------------------------------------------
library(AIPW)
library(SuperLearner)
library(ggplot2)
set.seed(123)
data("eager_sim_obs")
cov = c("eligibility","loss_num","age", "time_try_pregnant","BMI","meanAP")

AIPW_SL <- AIPW$new(Y= eager_sim_obs$sim_Y,
                    A= eager_sim_obs$sim_A,
                    W= subset(eager_sim_obs,select=cov), 
                    Q.SL.library = c("SL.glm"),
                    g.SL.library = c("SL.glm"),
                    k_split = 2,
                    verbose=TRUE)$
  fit()$
  summary()

## ----refit--------------------------------------------------------------------
# Create a new object from the previous AIPW_SL (Repeated class is an extension of the AIPW class)
repeated_aipw_sl <- Repeated$new(aipw_obj = AIPW_SL)
# Fit repetitively
repeated_aipw_sl$repfit(num_reps = 30, stratified = F)
# Summarise the median estimate, median SE, and the SE of median estimate adjusting for `num_reps` repetitions
repeated_aipw_sl$summary_median()

## ----check refit--------------------------------------------------------------
# Check the distributions of estiamtes from `num_reps` repetitions
s <- repeated_aipw_sl$repeated_estimates
ggplot2::ggplot(ggplot2::aes(x=Estimate),data = s) + ggplot2::geom_histogram(bins = 10) + ggplot2::facet_grid(~Estimand, scales = "free")
ggplot2::ggplot(ggplot2::aes(x=SE),data = s) + ggplot2::geom_histogram(bins = 10) + ggplot2::facet_grid(~Estimand, scales = "free")

