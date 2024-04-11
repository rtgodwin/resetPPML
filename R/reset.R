reset <- function(mod) {
  z <- list()
  
  #GLM versions
  
  #Wald
  z$Wald.GLM.Taylor.aug1 <- resetWald(mod, aug.terms = 1, robust = F, fourier = F)
  z$Wald.GLM.Taylor.aug2 <- resetWald(mod, aug.terms = 2, robust = F, fourier = F)
  tryCatch(
    resetWald(mod, aug.terms = 3, robust = F, fourier = F),
    warning = function(w) {z$Wald.GLM.Taylor.aug3 <<- "GLM failed to converge"})
  z$Wald.GLM.Taylor.aug4 <- resetWald(mod, aug.terms = 4, robust = F, fourier = F)
  
  z$Wald.GLM.Fourier.sinlink.aug2 <- resetWald(mod, aug.terms = 2, robust = F, fourier = T, sin.link = T)
  z$Wald.GLM.Fourier.sinlink.aug4 <- resetWald(mod, aug.terms = 4, robust = F, fourier = T, sin.link = T)
  
  z$Wald.GLM.Fourier.linlink.aug2 <- resetWald(mod, aug.terms = 2, robust = F, fourier = T, sin.link = F)
  z$Wald.GLM.Fourier.linlink.aug4 <- resetWald(mod, aug.terms = 4, robust = F, fourier = T, sin.link = F)
  
  z$Wald.GLM.combined.sinlink <- resetWaldTF(mod, robust = F, sin.link = T)
  z$Wald.GLM.combined.linlink <- resetWaldTF(mod, robust = F, sin.link = F)
  
  #LM
  z$lm.GLM.Taylor.aug1 <- resetlm(mod, aug.terms = 1, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug2 <- resetlm(mod, aug.terms = 2, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug3 <- resetlm(mod, aug.terms = 3, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug4 <- resetlm(mod, aug.terms = 4, robust = F, fourier = F)
  
  z$lm.GLM.Fourier.sinlink.aug2 <- resetlm(mod, aug.terms = 2, robust = F, fourier = T, sin.link = T)
  z$lm.GLM.Fourier.sinlink.aug4 <- resetlm(mod, aug.terms = 4, robust = F, fourier = T, sin.link = T)
  
  z$lm.GLM.Fourier.linlink.aug2 <- resetlm(mod, aug.terms = 2, robust = F, fourier = T, sin.link = F)
  z$lm.GLM.Fourier.linlink.aug4 <- resetlm(mod, aug.terms = 4, robust = F, fourier = T, sin.link = F)
  
  z$lm.GLM.combined.sinlink <- resetlmTF(mod, robust = F, sin.link = T)
  z$lm.GLM.combined.linlink <- resetlmTF(mod, robust = F, sin.link = F)
  
  #robust versions
  
  #Wald
  z$Wald.robust.Taylor.aug1 <- resetWald(mod, aug.terms = 1, robust = T, fourier = F)
  z$Wald.robust.Taylor.aug2 <- resetWald(mod, aug.terms = 2, robust = T, fourier = F)
  z$Wald.robust.Taylor.aug3 <- tryCatch(
    resetWald(mod, aug.terms = 3, robust = T, fourier = F),
    warning = function(w) {z$Wald.robust.Taylor.aug3 <<- "GLM failed to converge"})
  z$Wald.robust.Taylor.aug4 <- resetWald(mod, aug.terms = 4, robust = T, fourier = F)
  
  z$Wald.robust.Fourier.sinlink.aug2 <- resetWald(mod, aug.terms = 2, robust = T, fourier = T, sin.link = T)
  z$Wald.robust.Fourier.sinlink.aug4 <- resetWald(mod, aug.terms = 4, robust = T, fourier = T, sin.link = T)
  
  z$Wald.robust.Fourier.linlink.aug2 <- resetWald(mod, aug.terms = 2, robust = T, fourier = T, sin.link = F)
  z$Wald.robust.Fourier.linlink.aug4 <- resetWald(mod, aug.terms = 4, robust = T, fourier = T, sin.link = F)
  
  z$Wald.robust.combined.sinlink <- resetWaldTF(mod, robust = T, sin.link = T)
  z$Wald.robust.combined.linlink <- resetWaldTF(mod, robust = T, sin.link = F)
  
  #LM
  z$lm.robust.Taylor.aug1 <- resetlm(mod, aug.terms = 1, robust = T, fourier = F)
  z$lm.robust.Taylor.aug2 <- resetlm(mod, aug.terms = 2, robust = T, fourier = F)
  z$lm.robust.Taylor.aug3 <- resetlm(mod, aug.terms = 3, robust = T, fourier = F)
  z$lm.robust.Taylor.aug4 <- resetlm(mod, aug.terms = 4, robust = T, fourier = F)
  
  z$lm.robust.Fourier.sinlink.aug2 <- resetlm(mod, aug.terms = 2, robust = T, fourier = T, sin.link = T)
  z$lm.robust.Fourier.sinlink.aug4 <- resetlm(mod, aug.terms = 4, robust = T, fourier = T, sin.link = T)
  
  z$lm.robust.Fourier.linlink.aug2 <- resetlm(mod, aug.terms = 2, robust = T, fourier = T, sin.link = F)
  z$lm.robust.Fourier.linlink.aug4 <- resetlm(mod, aug.terms = 4, robust = T, fourier = T, sin.link = F)
  
  z$lm.robust.combined.sinlink <- resetlmTF(mod, robust = T, sin.link = T)
  z$lm.robust.combined.linlink <- resetlmTF(mod, robust = T, sin.link = F)
  
  return(z)
}