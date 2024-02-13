resetlm <- function(mod, aug.terms = 2, robust = T) {
  
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
  
  if(robust == F) {
    new_response_var <- "util"
    # Create the new formula string
    if(aug.terms == 1) {
      new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * log(yhat) ^ 2) - 1")
    } else if(aug.terms == 2) {
      new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * log(yhat) ^ 2) + I(w * log(yhat) ^ 3) - 1")
    } else if(aug.terms == 3) {
      new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * log(yhat) ^ 2) + I(w * log(yhat) ^ 3) + I(w * log(yhat) ^ 4) - 1")
    } else if(aug.terms == 4) {
      new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * log(yhat) ^ 2) + I(w * log(yhat) ^ 3) + I(w * log(yhat) ^ 4) + I(w * log(yhat) ^ 5) - 1")
    }
    new_formula <- as.formula(new_formula_string)
    
    aux <- lm(new_formula, data = mod$model)
    
    z$reset <- n * summary(aux)$r.squared
    z$pval <- 1 - pchisq(z$reset, aug.terms)
    
    return(z)
  }
  
  if(robust == T) {
    new_response_var1 <- "I(w * log(yhat) ^ 2)"
    new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula1 <- as.formula(new_formula_string1)
    aux1resids <- lm(new_formula1, data = mod$model)$residuals
    
    new_response_var2 <- "I(w * log(yhat) ^ 3)"
    new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
    new_formula2 <- as.formula(new_formula_string2)
    aux2resids <- lm(new_formula2, data = mod$model)$residuals
    
    ones <- rep(1, n)
    aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) -1)
    
    z$reset.robust <- n - sum(aux$residuals ^ 2)
    z$pval <- 1 - pchisq(z$reset.robust, 2)
    
    return(z)
  }
}