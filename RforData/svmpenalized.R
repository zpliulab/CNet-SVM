# function ----------------------------------------------------------------
## svmpenalized function 
## input:
svmpenalized <- function(data, methodlist){
  # data <- data
  # methodlist <- name
  
  condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)
  
  lambda <- c(seq(0.01 ,0.05, .01), seq(0.1,0.5, 0.2), 1 )
  lambda <-lambda[2:3]

  bounds <- t(data.frame(log2lambda1 = c(-10, 10)))
  colnames(bounds)<-c("lower", "upper")
  
  ## first
  set.seed(123*1)
  training.samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
  train.data  <- as.matrix(data[training.samples, ])
  test.data <- as.matrix(data[-training.samples, ])  
  
  ## train
  x <- train.data[,-1]    # sample*gene
  y <-  2*as.numeric(as.factor(train.data[,1]))-3 
  x_test <- test.data[,-1]    # sample*gene
  y_test = 2*as.numeric(as.factor(test.data[,1]))-3 
  
  
  if (methodlist == "scad+L2") {
    print("start Fixed grid -- scad+L2 method")
    fit <- svmfs(x, y, fs.method = methodlist,
                        cross.outer = 0, 
                        grid.search = "discrete",
                        lambda1.set = lambda,
                        parms.coding = "none",
                        show = "none",
                        maxIter = 10, 
                        inner.val.method = "cv", 
                        cross.inner = 5,
                        maxevals = 500,
                        seed = seed, 
                        verbose = FALSE)
  } else {
    print("start interval search -- XXXXmethod")
    fit <- svmfs(x, y, fs.method = methodlist,
                 cross.outer = 0, 
                 bounds = bounds,
                 gird.search = "interval",
                 parms.coding = "log2",  
                 show = "none",
                 maxIter = 10,    # 700
                 inner.val.method = "cv", 
                 cross.inner = 5,
                 maxevals = 500,
                 seed = seed, 
                 verbose = FALSE) 
}
  # SCAD result
  fea <- c()
  fea <- c(fea, list(colnames(x)[fit$model$xind]))
  
  # predict
  test <- predict(fit, x_test, y_test)
  pred <- cbind(y_test, as.numeric(test$fitted))
  p <- plot.roc(y_test, as.numeric(test$fitted), print.auc=T)
  auclist <- as.numeric(p$auc)

  ## second
  for (i in condset) {
    print(paste("******* i= ********",i))
    # i = 7
    set.seed(123*i)
    training.samples <- data$outcome %>% createDataPartition(p=0.7, list = FALSE)
    train.data  <- as.matrix(data[training.samples, ])
    test.data <- as.matrix(data[-training.samples, ])  
    
    ## train
    x <- train.data[,-1]    # sample*gene
    y <-  2*as.numeric(as.factor(train.data[,1]))-3 
    x_test <- test.data[,-1]    # sample*gene
    y_test = 2*as.numeric(as.factor(test.data[,1]))-3 

    ## four method
    if (methodlist == "scad+L2") {
      print("start Fixed grid -- scad+L2 -- method")
      fit <- svmfs(x, y, fs.method = methodlist,
                   cross.outer = 0, 
                   grid.search = "discrete",
                   lambda1.set = lambda,
                   parms.coding = "none",
                   show = "none",
                   maxIter = 10, 
                   inner.val.method = "cv", 
                   cross.inner = 5,
                   maxevals = 500,
                   seed = seed, 
                   verbose = FALSE)
    } else {
      print("start interval search -- XXXX -- method")
	    fit <- svmfs(x, y, fs.method = methodlist, 
	                 cross.outer = 0, 
	                 bounds = bounds,
	                 gird.search = "interval",
	                 parms.coding = "log2",  
	                 show = "none",
	                 maxIter = 10,    # 700
	                 inner.val.method = "cv", 
	                 cross.inner = 5,
	                 maxevals = 500,
	                 seed = seed, 
	                 verbose = FALSE) 
    }
    
    ## Integration
    fea <- c(fea, list(colnames(x)[fit$model$xind]))
    # predict
    test <- predict(fit, x_test, y_test)
    pred <- cbind(pred, cbind(y_test, as.numeric(test$fitted)))
    p <- plot.roc(y_test, as.numeric(test$fitted), print.auc=T)
    auclist <- c(auclist, as.numeric(p$auc))
    print(paste("*************** i= *************",i)) 
  }
  return(list(fea, auclist, pred))
}