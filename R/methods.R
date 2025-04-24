## --list with allowable gps methods and their arguments------------------------
.gps_methods <- list(
  "multinom" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = c("mlreg", "mnom"),
    packages_needed = c("nnet"),
    fun.arg.check = list(
      "nnet::multinom",
      "nnet::nnet.formula"
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
    fun.arg.check = list(
      "MASS::polr"
    ),
    link_fun = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
    description = "estimating gps for ordered treatments using proportional
                odds logistic regression from MASS package"
  ),
  "vglm" = list(
    missing = c("complete.cases", "mean.imputation"),
    func_used = ".estimate_gps_multinom",
    alias = c("vecGLM"),
    fun.arg.check = list(
      "VGAM::vglm",
      "VGAM::rrvglm"
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
    fun.arg.check = list("brglm2::brmultinom"),
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
      "mclogit::mblogit"
    ),
    packages_needed = "mclogit",
    link_fun = c("baseline_category_logit"),
    allowed.treat = c("binary", "multinom"),
    description = "baseline-category logit models"
  )
)

## --list with allowable significance methods and their arguments---------------
.sig_methods <- list(
  "t_test" = list(
    method_name = "t_test",
    package_used = "rstatix",
    args_check_fun = list(
      rstatix::pairwise_t_test
    )
  ),
  "dunn_test" = list(
    method_name = "dunn_test",
    package_used = "rstatix",
    args_check_fun = list(
      rstatix::dunn_test
    )
  ),
  "tukeyHSD_test" = list(
    method_name = "tukeyHSD",
    package_used = "rstatix",
    args_check_fun = list(
      utils::getS3method("tukey_hsd", "data.frame",
        envir = asNamespace("rstatix")
      )
    )
  ),
  "games_howell_test" = list(
    method_name = "games_howell_test",
    package_used = "rstatix",
    args_check_fun = list(
      rstatix::games_howell_test
    )
  ),
  "wilcoxon_test" = list(
    method_name = "wilcoxon_test",
    package_used = "rstatix",
    args_check_fun = list(
      rstatix::pairwise_wilcox_test
    )
  ),
  "sign_test" = list(
    method_name = "sign_test",
    package_used = "rstatix",
    args_check_fun = list(
      rstatix::pairwise_sign_test,
      rstatix::sign_test
    )
  )
)

## --list with all allowable matching methods and their arguments---------------
.match_methods <- list(
  "nnm" = list(
    allowed_args = c(
      "caliper", "ratio", "replace", "ties", "order"
    ),
    args_check_fun = list(
      Matching::Match,
      Matching::Matchby
    ),
    matching_fun = Matching::Matchby,
    formula_necessary = FALSE,
    data_name = "X",
    treat_var = TRUE
  ),

  # placeholder for Matching::Match() - used only when method = "nnm" and two
  # treatment groups are present
  "nnm_2t" = list(
    allowed_args = c(
      "caliper", "ratio", "replace", "ties", "order"
    ),
    args_check_fun = list(
      Matching::Match
    ),
    matching_fun = Matching::Match,
    formula_necessary = FALSE,
    data_name = "X",
    treat_var = TRUE
  ),
  "fullopt" = list(
    allowed_args = c(
      "caliper", "order", "min_controls", "max_controls"
    ),
    args_check_fun = list(
      optmatch::match_on,
      optmatch::fullmatch
    ),
    matching_fun = optmatch::fullmatch,
    formula_necessary = TRUE,
    data_name = "data",
    treat_var = FALSE
  ),
  "pairopt" = list(
    allowed_args = c(
      "caliper", "order", "ratio"
    ),
    args_check_fun = list(
      optmatch::match_on,
      optmatch::pairmatch
    ),
    matching_fun = optmatch::pairmatch,
    formula_necessary = TRUE,
    data_name = "data",
    treat_var = FALSE
  )
)
