feeplot <- function(model, data, maxpred, ylimit, ccex) {
  
  plotpp <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=25, col="darkmagenta",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="darkmagenta",lwd=ccex)
    leg <<- c(leg, "counterfactual distribution")
    cols <<- c(cols, "darkmagenta")
    pchs <<- c(pchs, 25)
  }
  
  plotztnb <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=18, col="red",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="red",lwd=ccex)
    leg <<- c(leg, "counterfactual distribution")
    cols <<- c(cols, "red")
    pchs <<- c(pchs, 18)
  }
  
  plotoipp <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=17, col="green",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="green",lwd=ccex)
    leg <<- c(leg, "factual distribution")
    cols <<- c(cols, "green")
    pchs <<- c(pchs, 17)
  }
  
  plotoiztnb <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=16, col="blue",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="blue",lwd=ccex)
    leg <<- c(leg, "factual distribution")
    cols <<- c(cols, "blue")
    pchs <<- c(pchs, 16)
  }
  
  #b <- model$beta
  #g <- model$gamma
  #if (model$dist == "negbin") {a <- model$alpha}
  
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  #X <- cleandata$X
  #Z <- cleandata$Z
  y <- cleandata$y
  
  if(missing(maxpred)) {
    maxpred = max(y)
  }
  
  if(missing(ylimit)) {
    ylimit = max(tabulate(y)) * 1.1
  }
  
  if(missing(ccex)) {
    ccex = 1.5
  }
  
  df.bar <- barplot(tabulate(y)[1:maxpred], names=1:maxpred, xlab="count", ylab="frequency", col="gray", ylim = c(0, ylimit))
  leg <- "actual data"
  cols <- "gray"
  pchs <- 15
  
  if(class(model) == "oneinflmodel" & model$dist == "Poisson") {
    plotoipp(model, data, maxpred)
    modelcounter <- model
    class(modelcounter) <- "truncmodel"
    plotpp(modelcounter, data, maxpred)
  } else if(class(model) == "oneinflmodel" & model$dist == "negbin") {
    plotoiztnb(model, data, maxpred)
    modelcounter <- model
    class(modelcounter) <- "truncmodel"
    plotztnb(modelcounter, data, maxpred)
  }
  
  legend("topright", legend=leg, col=cols, pch=pchs, cex = 1)

}
