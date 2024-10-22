## --list with allowable gps methods and their arguments-------------------------
.gps_methods <- list(
  "multinom" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = c("mlreg", "mnom"),
    packages_needed = c("nnet"),
    fun.arg.check = list(
      nnet::multinom,
      nnet::nnet.formula
    ),
    link_fun = c("generalized_logit"),
    allowed.treat = c("binary", "multinom"),
    description = c("estimating the GPS using multinomial logistic
                               regression model from nnet package")
  ),
  "polr" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = "propodds",
    packages_needed = "MASS",
    link_fun = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
    description = "estimating gps for ordered treatments using proportional
                odds logistic regression from MASS package"
  ),
  "vglm" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = c("vecGLM"),
    fun.arg.check = list(
      VGAM::vglm,
      VGAM::rrvglm
    ),
    packages_needed = "VGAM",
    link_fun = c("multinomial_logit", "reduced_rank_ml"),
    allowed.treat = c("binary", "multinom"),
    description = "vector generalized linear models for multinomial data"
  ),
  "brglm2" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = "bias_reduced_glm2",
    fun.arg.check = list(brglm2::brmultinom),
    packages_needed = "brglm2",
    link_fun = "baseline_category_logit",
    allowed.treat = c("binary", "multinom"),
    description = "bias reduction for multinomial respones models
                  using the poisson trick"
  ),
  "mblogit" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = "multinomial_logit_model",
    fun.arg.check = list(
      mclogit::mblogit
    ),
    packages_needed = "mclogit",
    link_fun = c("baseline_category_logit"),
    allowed.treat = c("binary", "multinom"),
    description = "baseline-category logit models"
  )
)
