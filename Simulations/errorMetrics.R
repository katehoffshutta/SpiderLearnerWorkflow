library(MatrixCorrelation)

relativeFrobNormAfter = function(estimatedMatrix, trueMatrix)
{
  Delta = as.matrix(estimatedMatrix - trueMatrix)
  return(norm(Delta,"F")/norm(as.matrix(trueMatrix),"F"))
}


ellOneNormLowerTri = function(estimatedMatrix, trueMatrix)
{
  Delta = as.matrix(estimatedMatrix - trueMatrix)
  
  return(sum(abs(Delta[lower.tri(Delta)])))
}


squaredDiffLowerTri = function(estimatedMatrix, trueMatrix)
{
  Delta = as.matrix(estimatedMatrix - trueMatrix)
  
  return(sum((Delta[lower.tri(Delta)])^2))
}

loglikLossfunction = function(thetaEst,dataTest)
{
  n = nrow(dataTest)
  p = ncol(dataTest)
  firstTerm = -n/2*log(det(solve(thetaEst)))
  secondTerm= 0
  for(i in 1:n)
    secondTerm = secondTerm - 1/2*t(dataTest[i,]) %*% thetaEst %*% dataTest[i,]
  return(firstTerm + secondTerm)
}

matrixRV = function(estimatedMatrix, trueMatrix)
{
  return(RV(estimatedMatrix,trueMatrix))
}

deltaMatrix = function(estimatedMatrix, trueMatrix)
{
  # returns a matrix, entries are estimated - true
  return(estimatedMatrix - trueMatrix) 
}

thresholdMatrix = function(thisMat, threshold)
{
  thisMatThres = ifelse(abs(thisMat) <= threshold, 0, 1)
  return(thisMatThres)
}

binaryMatrixDiff = function(estimatedMatrix, trueMatrix, threshold)
{
  return(sum(abs(thresholdMatrix(estimatedMatrix,threshold) - thresholdMatrix(trueMatrix,0))))
}

truePositiveCount= function(estimatedMatrix, trueMatrix, threshold)
{
  return(sum(thresholdMatrix(estimatedMatrix,threshold) == 1 & thresholdMatrix(trueMatrix,0) == 1))
}

falsePositiveCount= function(estimatedMatrix, trueMatrix, threshold)
{
  return(sum(thresholdMatrix(estimatedMatrix,threshold) == 1 & thresholdMatrix(trueMatrix,0) == 0))
}

trueNegativeCount= function(estimatedMatrix, trueMatrix, threshold)
{
  return(sum(thresholdMatrix(estimatedMatrix,threshold) == 0 & thresholdMatrix(trueMatrix,0) == 0))
}

falseNegativeCount= function(estimatedMatrix, trueMatrix, threshold)
{
  return(sum(thresholdMatrix(estimatedMatrix,threshold) == 0 & thresholdMatrix(trueMatrix,0) == 1))
}

logLik = function(thetaEst,dataTest)
{
  n = nrow(dataTest)
  p = ncol(dataTest)
  firstTerm = -n/2*log(det(solve(thetaEst)))
  secondTerm= 0
  for(i in 1:n)
    secondTerm = secondTerm - 1/2*t(dataTest[i,]) %*% thetaEst %*% dataTest[i,]
  return(firstTerm + secondTerm)
}

# make a matrix that specifies whether the entry of the input matrix
# is in each of the following groups:
# 0: a zero entry
# 1: first quartile of off-diagonal elements
# 2: second or third quartile of off-diagonal elements
# 3: fourth quartile of off-diagonal elements
# 4: diagonal elements

splitToFive = function(M)
{
  MAbs = as.matrix(abs(M))
  offDiagonal = lower.tri(MAbs) | upper.tri(MAbs)
  category0 = (MAbs == 0)
  MBlank = MAbs-diag(diag(MAbs))
  if(sum(diag(MBlank))!=0) print("!ERROR!")
  offDiagonalCutoffs = quantile(MBlank[MBlank!=0], c(0.25,0.75))
  category1 = (0 < MBlank & MBlank <= offDiagonalCutoffs[1])
  category2 = (offDiagonalCutoffs[1] < MBlank & MBlank <= offDiagonalCutoffs[2])
  category3 = (offDiagonalCutoffs[2] < MBlank)
  category4 = !(lower.tri(MAbs) | upper.tri(MAbs))
  outMat = ifelse(category0,0,
                  ifelse(category1,1,
                         ifelse(category2,2,
                                ifelse(category3,3,
                                       ifelse(category4,4,NA)))))
  return(outMat)
}
