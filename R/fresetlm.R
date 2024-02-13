fresetlm <- function(mod, aug.terms = 2, robust = T) {
  
  yhat <- mod$fitted.values
  y <- mod$y
  uhat <- y - yhat
  w <- sqrt(yhat)
  util <- uhat / w
  n <- nrow(mod$model)
  z <- list()
  
  original_formula <- mod$formula
  response_var <- all.vars(original_formula)[1]
  predictor_vars <- all.vars(original_formula)[-1] # Exclude the first variable (response)
  
  wver <- pi * (2 * (yhat) - (max(yhat) - min(yhat))) / (max(yhat) - min(yhat))
  vver <- 2 * pi * ((sin(yhat)) ^ 2) - pi
  
  if(robust == F) {
    new_response_var <- "util"
    
    # w version
    new_formula_stringw <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * sin(wver)) + I(w * cos(wver)) - 1")
    new_formulaw <- as.formula(new_formula_stringw)
    
    auxw <- lm(new_formulaw, data = mod$model)
    
    z$fresetw <- n * summary(auxw)$r.squared
    z$pvalw <- 1 - pchisq(z$fresetw, aug.terms)
    
    # v version
    new_formula_stringv <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * sin(vver)) + I(w * cos(vver)) - 1")
    new_formulav <- as.formula(new_formula_stringv)
    
    auxv <- lm(new_formulav, data = mod$model)
    
    z$fresetv <- n * summary(auxv)$r.squared
    z$pvalv <- 1 - pchisq(z$fresetv, aug.terms)
    
    return(z)
  }
  
  if(robust == T) {
    # w version
    new_response_var1 <- "I(w * sin(wver))"
    new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula1 <- as.formula(new_formula_string1)
    aux1resids <- lm(new_formula1, data = mod$model)$residuals
    
    new_response_var2 <- "I(w * cos(wver))"
    new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula2 <- as.formula(new_formula_string2)
    aux2resids <- lm(new_formula2, data = mod$model)$residuals
    
    ones <- rep(1, n)
    aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) -1)
    
    z$reset.robust.w <- n - sum(aux$residuals ^ 2)
    z$pval.robust.w <- 1 - pchisq(z$reset.robust.w, aug.terms)
    
    # v version
    new_response_var1 <- "I(w * sin(vver))"
    new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula1 <- as.formula(new_formula_string1)
    aux1resids <- lm(new_formula1, data = mod$model)$residuals
    
    new_response_var2 <- "I(w * cos(vver))"
    new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula2 <- as.formula(new_formula_string2)
    aux2resids <- lm(new_formula2, data = mod$model)$residuals
    
    ones <- rep(1, n)
    aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) -1)
    
    z$reset.robust.v <- n - sum(aux$residuals ^ 2)
    z$pval.robust.v <- 1 - pchisq(z$reset.robust.v, aug.terms)
    
    return(z)
  }
}