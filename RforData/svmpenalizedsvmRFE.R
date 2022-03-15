## svmpenalized mRMR£¬a kind of filter method, set feature num = 30


svmpenalizedsvmRFE <- function(data){
  condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)
  ## the length of solution -- the number of solution
  solength <- 30
  
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
  
  
  # SVM-RFE -----------------------------------------------------------------
  # result <- svmRFE(train.data, k=1, halve.above=100)    # standard SVM-RFE, you can use k=1.
  featureRankedList = svmrfeFeatureRanking(x,y)
  rank <- featureRankedList[1:solength]
  feature <- colnames(x)[rank]
  
  
  # fit
  svmfit = svm(x[ , rank], y, cost = 10, kernel="linear")  # linear svm
  summary(svmfit)
  
  # result
  fea <- c()
  fea <- c(fea, list(feature))
  
  # predict
  test <- predict(svmfit, x_test[ , rank])
  pred <- cbind(y_test, as.numeric(test))
  p <- plot.roc(y_test, as.numeric(test), print.auc=T)
  auclist <- as.numeric(p$auc)

  ## second
  for (i in condset) {
    print(paste("******* i= ********",i))
    set.seed(123*i)
    training.samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
    train.data  <- as.matrix(data[training.samples, ])
    test.data <- as.matrix(data[-training.samples, ]) 
    
    ## train
    x <- train.data[,-1]    # sample*gene
    y <-  2*as.numeric(as.factor(train.data[,1]))-3 
    x_test <- test.data[,-1]    # sample*gene
    y_test = 2*as.numeric(as.factor(test.data[,1]))-3 
    
    
    # SVM-RFE -----------------------------------------------------------------
    # result <- svmRFE(train.data, k=1, halve.above=100)    # standard SVM-RFE, you can use k=1.
    featureRankedList = svmrfeFeatureRanking(x,y)
    rank <- featureRankedList[1:solength]
    feature <- colnames(x)[rank]
    
    # fit
    svmfit = svm(x[ , rank], y, cost = 10, kernel="linear")  # linear svm
    summary(svmfit)
    
    ## Integration
    fea <- c(fea, list(feature))
    
    # predict
    test <- predict(svmfit, x_test[ , rank])
    pred <- cbind(pred, cbind(y_test, as.numeric(test)))
    p <- plot.roc(y_test, as.numeric(test), print.auc=T)
    auclist <- c(auclist, as.numeric(p$auc))
    print(paste("*************** i= *************",i)) 
  }
  return(list(fea, auclist, pred))
}
