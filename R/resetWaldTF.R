resetWaldTF<- function(mod, robust = T, sin.link = T, bothlinks = F) {
  
  etahat <- predict(mod, type = "link")
  original_formula <- mod$formula
  dat <- model.frame(mod)
  dat$etahat <- etahat
  z <- list()
  aug.terms = 4
  
  if(bothlinks == F) {
    
    if(sin.link == T) {
      dat$v <- 2 * pi * sin(etahat) ^ 2 - pi
    }
    
    if(sin.link == F) {
      dat$v <- pi * (2 * (etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
    }
    
    new_formula <- update(original_formula, . ~ . + sin(v) + cos(v))
    new_formula <- update(new_formula, . ~ . + I(etahat^2) + I(etahat^3))
    
    aux <- glm(new_formula, data = dat, family = stats::quasipoisson(link = "log"))
    
  }
  
  if(bothlinks == T) {
    
    aug.terms = 6
    
    dat$v <- 2 * pi * sin(etahat) ^ 2 - pi
    dat$z <- pi * (2 * (etahat) - (max(etahat) - min(etahat))) / (max(etahat) - min(etahat))
    
    new_formula <- update(original_formula, . ~ . + sin(v) + cos(v) + sin(z) + cos(z))
    new_formula <- update(new_formula, . ~ . + I(etahat^2) + I(etahat^3))
    
    aux <- glm(new_formula, data = dat, family = stats::quasipoisson(link = "log"))
    
  }
  
  if(robust == T) {
    temp_result <- try(lmtest::waldtest(mod, aux, test = "Chisq", vcov = sandwich::vcovHC(aux, type = "HC0")), silent = TRUE)
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
