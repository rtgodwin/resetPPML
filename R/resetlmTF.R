resetlmTF <- function(mod, robust = T, sin.link = T, bothlinks = F) {
  
  yhat <- mod$fitted.values
  etahat <- predict(mod, type = "link")
  y <- mod$y
  uhat <- y - yhat
  w <- sqrt(yhat)
  util <- uhat / w
  n <- nrow(mod$model)
  z <- list()
  aug.terms = 4
  
  original_formula <- mod$formula
  response_var <- all.vars(original_formula)[1]
  predictor_vars <- all.vars(original_formula)[-1] # Exclude the first variable (response)
  
  if(bothlinks == F) {
    
    if(sin.link == T) {
      v <- 2 * pi * ((sin(etahat)) ^ 2) - pi
    }
    
    if(sin.link == F) {
      v <- pi * (2 * (etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
    }
    
    if(bothlinks == T) {
      aug.terms = 6
      v <- 2 * pi * ((sin(etahat)) ^ 2) - pi
      z <- pi * (2 * (etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
    }
    
    if(robust == T) {
      
      new_response_var1 <- "I(w * sin(v))"
      new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
      new_formula1 <- as.formula(new_formula_string1)
      aux1resids <- lm(new_formula1, data = mod$model)$residuals
      
      new_response_var2 <- "I(w * cos(v))"
      new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
      new_formula2 <- as.formula(new_formula_string2)
      aux2resids <- lm(new_formula2, data = mod$model)$residuals
      
      new_response_var3 <- "I(w * etahat ^ 4)"
      new_formula_string3 <- paste(new_response_var3, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
      new_formula3 <- as.formula(new_formula_string3)
      aux3resids <- lm(new_formula3, data = mod$model)$residuals
      
      new_response_var4 <- "I(w * etahat ^ 5)"
      new_formula_string4 <- paste(new_response_var4, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
      new_formula4 <- as.formula(new_formula_string4)
      aux4resids <- lm(new_formula4, data = mod$model)$residuals
      
      if(bothlinks == T) {
        
        new_response_var5 <- "I(w * sin(z))"
        new_formula_string5 <- paste(new_response_var5, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula5 <- as.formula(new_formula_string5)
        aux5resids <- lm(new_formula5, data = mod$model)$residuals
        
        new_response_var6 <- "I(w * cos(z))"
        new_formula_string6 <- paste(new_response_var6, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula6 <- as.formula(new_formula_string6)
        aux6resids <- lm(new_formula6, data = mod$model)$residuals
      }
      
      ones <- rep(1, n)
      
      if(bothlinks == F) {
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) + I(util * aux3resids) + I(util * aux4resids) -1)
      }
      
      if(bothlinks == T) {
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) + I(util * aux3resids) + I(util * aux4resids) + I(util * aux5resids) + I(util * aux6resids) -1)
      }
      
            z$reset <- n - sum(aux$residuals ^ 2)
      z$pval <- 1 - pchisq(z$reset, aug.terms)
      
    }
    
    if(robust == F) {
      new_response_var <- "util"
      
      if(bothlinks == F) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * sin(v)) + I(w * cos(v)) + I(w * etahat ^ 2) + I(w * etahat ^ 3) -1")
      }
      
      if(bothlinks == T) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), 
                                    "+ I(w * sin(v)) + I(w * cos(v)) + I(w * etahat ^ 2) + I(w * etahat ^ 3) + I(w * sin(z)) + I(w * cos(z)) -1")
      }
      
      new_formula <- as.formula(new_formula_string)
      
      aux <- lm(new_formula, data = mod$model)
      
      z$reset <- n * summary(aux)$r.squared
      z$pval <- 1 - pchisq(z$reset, aug.terms)
    }
  
  }
  
  return(z)
}
