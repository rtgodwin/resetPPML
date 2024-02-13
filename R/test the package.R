library(wooldridge)
data("crime1")
mod1 <- glm(narr86 ~ pcnv + avgsen + tottime + ptime86 + qemp86 + inc86 
            + durat + black + hispan + born60, family=poisson, data = crime1)
resetlm(mod1)
resetlm(mod1, aug.terms = 3)
fresetlm(mod1)

mod2 <- glm(narr86 ~ pcnv + I(pcnv^2) + avgsen + tottime + ptime86 
            + I(ptime86^2) + qemp86 + inc86 + I(inc86^2) + durat + black 
            + hispan + born60, family=poisson, data = crime1)

resetlm(mod2)
resetlm(mod2, robust = F)
fresetlm(mod2)