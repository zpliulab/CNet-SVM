#########################################################################
## Name:        svmrfeFeatureRanking.R
## Author:      Guyon et al.
##
##
## Change Log:
##   * April 4, 2018
##       - split sigFeature.R into several files for easier readability
##         and maintenance
#########################################################################

#########################################################################
## svmrfeFeatureRanking():  This function was proposed by Guyon, Isabelle,
## et al. named SVM-RFE to solve the classification problem with the help
## of ranking the features. In this algorithm, the dataset has been
## trained with SVM linear kernel model and removed the feature with the
## smallest ranking criterion. This criterion is the w value of the
## decision hyperplane given by the SVM. The function is included in this
## package to make a comparison with sigFeature() function.
##
##
#########################################################################


svmrfeFeatureRanking = function(x,y){

    #Checking for the variables
    stopifnot(!is.null(x) == TRUE, !is.null(y) == TRUE)

    n = ncol(x)
    survivingFeaturesIndexes = seq_len(n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n

    while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,
                scale=FALSE, type="C-classification", kernel="linear" )

    #compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV

    #compute ranking criteria
    rankingCriteria = w * w

    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix

    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1

    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    }

    return (featureRankedList)
}
