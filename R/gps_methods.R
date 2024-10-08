##--list with allowable gps methods and their arguments-------------------------
.gps_methods <- list(
  'multinom' = list(missing = c('complete.cases', 'mean.imputation'),
                    func_used = '.estimate_gps_multinom',
               alias = c('mlreg', 'mnom'),
               packages_needed = c('nnet'),
               fun.arg.check = list(quote(nnet::multinom),
                                    quote(nnet::nnet.default),
                                    quote(nnet::nnet.formula)),
               link_fun = c('generalized_logit'),
               description = c('estimating the GPS using multinomial logistic
                               regression model from nnet package')),
  'polr' = list(missing = c('complete.cases', 'mean.imputation'),
                func_used = '.estimate_gps_multinom',
                alias = 'propodds',
                packages_needed = 'MASS',
                link_fun = c('logistic', 'probit', 'loglog', 'cloglog', 'cauchit'),
                description = 'estimating gps for ordered treatments using proportional
                odds logistic regression from MASS package'),
  'nnet' = list(missing = c('complete.cases', 'mean.imputation'),
                func_used = '.estimate_gps_multinom',
                alias = 'neuralnet',
                packages_needed = 'nnet',
                link_fun = c('softmax'),
                description = 'estimating gps for ordered treatments using proportional
                odds logistic regression from MASS package')
)
