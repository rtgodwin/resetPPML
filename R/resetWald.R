resetWald<- function(mod, aug.terms = 2, robust = T, fourier = T, sin.link = T) {
  
  etahat <- predict(mod, type = "link")
  original_formula <- mod$formula
  dat <- model.frame(mod)
  dat$etahat <- etahat
  z <- list()

  if(fourier == T) {
    
    if(sin.link == T) {
      dat$v <- 2 * pi * sin(etahat) ^ 2 - pi
    }
    
    if(sin.link == F) {
      dat$v <- pi * (2 * (etahat) - (max(etahat) - min(etahat))) / (max(etahat) - min(etahat))
    }
    
    # Add new variables to the predictors part
    if(aug.terms == 1) {stop("Augmentation terms must be a multiple of 2 for the Fourier transform")}
    if(aug.terms == 2) {new_formula <- update(original_formula, . ~ . + sin(v) + cos(v))}
    if(aug.terms == 3) {stop("Augmentation terms must be a multiple of 2 for the Fourier transform")}
    if(aug.terms == 4) {new_formula <- update(original_formula, . ~ . + sin(v) + cos(v) + I(sin(2*v)) + I(cos(2*v)))}
  }
  
  if(fourier == F) {
    # Add nev variables to the predictors part
    if(aug.terms == 1) {new_formula <- update(original_formula, . ~ . + I(etahat^2))}
    if(aug.terms == 2) {new_formula <- update(original_formula, . ~ . + I(etahat^2) + I(etahat^3))}
    if(aug.terms == 3) {new_formula <- update(original_formula, . ~ . + I(etahat^2) + I(etahat^3) + I(etahat^4))}
    if(aug.terms == 4) {new_formula <- update(original_formula, . ~ . + I(etahat^2) + I(etahat^3) + I(etahat^4) + I(etahat^5))}
  }
  
  aux <- glm(new_formula, data = dat, family = stats::quasipoisson(link = "log"))
  
  if(robust == T) {
    temp_result <- try(lmtest::waldtest(mod, aux, test = "Chisq", vcov = sandwich::vcovHC(aux, type = "HC1")), silent = TRUE)
    if (inherits(temp_result, "try-error")) {
      z$reset <- 1000
      z$pval <- 0
    } else {
      z$reset <- temp_result$Chisq[2]
      z$pval <- temp_result$"Pr(>Chisq)"[2]
    }
  }
  
  if(robust == F) {
    temp_result <- try(lmtest::waldtest(mod, aux, test = "Chisq"), silent = TRUE)
    if (inherits(temp_result, "try-error")) {
      z$reset <- 1000
      z$pval <- 0
    } else {
      z$reset <- temp_result$Chisq[2]
      z$pval <- temp_result$"Pr(>Chisq)"[2]
    }
  }
  return(z)
}
