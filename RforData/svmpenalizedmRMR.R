## svmpenalized mRMR£¬a kind of filter method, set feature num = 30


svmpenalizedmRMR <- function(data){
  condset <- c(3,13,14,17:20,23:25,28,32,33,36,38:41,44,46)
  ## the length of solution -- the number of solution
  solength <- 30
  
  ## first
  set.seed(123*1)
  training_samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
  ##  # sample*gene(col 1 -label)
  train_data  <- as.matrix(data[training_samples, ])  
  test_data <- as.matrix(data[-training_samples, ])  
  ## data.fram  --  train and test
  train_data <- data.frame(train_data) 
  test_data <- data.frame(test_data) 
  
  mrmre_train <- mRMR.data(data=train_data[,-1], strata = factor(train_data[,1]))
  
  ## train -- classical mRMR feature selection
  mdata <- mrmre_train
  classic.m <- mRMR.classic(mdata, feature_count=solength, target_indices=1)
  ## unable to find an inherited method for function ¡®featureNames¡¯ for signature 
  sefs_rmre <- mdata@feature_names[solutions(classic.m)[[1]]]
  
  # fit ---------------------------------------------------------------------
  ens1_fs <- as.matrix(sefs_rmre)
  ## ens1_fs -- must have dim (matrix) 
  model <- apply(ens1_fs, 2, function(x, y) {
    ff <- as.formula(sprintf("%s ~ %s", colnames(y)[1], paste(x, collapse=" + ")))
    mm <- lm(formula=ff, data=y, model=FALSE)
    return(mm)
  }, y=train_data)
  
  
  # mRMR result -------------------------------------------------------------
  fea <- c()
  fea <- c(fea, list(sefs_rmre))
  
  # predict ---------------------------------------------------------------------
  test <- predict(object=model[[1]], newdata=test_data, type="response")
  pred <- cbind(test_data[,1], as.numeric(test))
  p <- plot.roc(test_data[,1], as.numeric(test), print.auc=T)
  auclist <- as.numeric(p$auc)
  
  ## second
  for (i in condset) {
    print(paste("******* i= ********",i))
    set.seed(123*i)
    training_samples <- data$outcome %>% createDataPartition(p = 0.7, list = FALSE)
    ##  # sample*gene(col 1 -label)
    train_data  <- as.matrix(data[training_samples, ])  
    test_data <- as.matrix(data[-training_samples, ])  
    ## data.fram  --  train and test
    train_data <- data.frame(train_data) 
    test_data <- data.frame(test_data) 
    
    mrmre_train <- mRMR.data(data=train_data[,-1], strata = factor(train_data[,1]))
    
    ## train -- classical mRMR feature selection
    mdata <- mrmre_train
    classic.m <- mRMR.classic(mdata, feature_count=solength, target_indices=1)
    ## unable to find an inherited method for function ¡®featureNames¡¯ for signature 
    sefs_rmre <- mdata@feature_names[solutions(classic.m)[[1]]]
    
    # fit ---------------------------------------------------------------------
    ens1_fs <- as.matrix(sefs_rmre)
    ## ens1_fs -- ±ØÐëÓÐÎ¬¶È£¨¾ØÕó£©
    model <- apply(ens1_fs, 2, function(x, y) {
      ff <- as.formula(sprintf("%s ~ %s", colnames(y)[1], paste(x, collapse=" + ")))
      mm <- lm(formula=ff, data=y, model=FALSE)
      return(mm)
    }, y=train_data)
    
    ## Integration
    fea <- c(fea, list(sefs_rmre))
    # predict
    test <- predict(object=model[[1]], newdata=test_data, type="response")
    pred <- cbind(pred, cbind(test_data[,1], as.numeric(test)))
    p <- plot.roc(test_data[,1], as.numeric(test), print.auc=T)
    auclist <- c(auclist, as.numeric(p$auc))
    print(paste("*************** i= *************",i)) 
  }
  return(list(fea, auclist, pred))
}



