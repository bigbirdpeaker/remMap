\name{remMap.df}
\alias{remMap.df}
\title{A function to calculate the degrees of freedom of the remMap estimator.}
\description{A function to calculate the degrees of freedom of the remMap estimator. 
It assumes orthogonality of the design matrix to estimate the degrees of freedom 
and tends to overestimate the actual degrees of freedom when the design matrix is far from orthogonal.}
   
\usage{
remMap.df(X.m, Y.m, lamL1.v, lamL2.v, C.m=NULL)
}

\arguments{
     \item{X.m}{numeric matrix (n by p): columns correspond to predictor variables and rows correspond to samples. 
                        Missing values are not allowed.}  
       \item{Y.m}{numeric matrix (n by q): columns correspond to response variables and rows correspond to samples. 
                     Missing values are not allowed.}                   
    \item{lamL1.v}{numeric vector: a set of \eqn{l_1} norm penalty parameters.} 
    \item{lamL2.v}{numeric vector: a set of \eqn{l_2} norm penalty parameters.}    
   \item{C.m}{numeric matrix (p by q).  \eqn{C_m[i,j]=0} means the corresponding coefficient beta[i,j] is set to be zero in the model;
                        \eqn{C_m[i,j]=1} means the corresponding beta[i,j] is included in the MAP penalty; 
                        \eqn{C_m[i,j]=2} means the corresponding beta[i,j] is not included in the MAP penalty; 
                  default(=NULL): \eqn{C_m[i,j]} are all set to be 1.}
          }        
\details{
    \code{remMap.df} calculates the estimated degrees of freedom of the remMap estimator for each pair of tuning parameters (lamL1, lamL2) (Peng and et.al., 2008).}

\value{
    A numeric matrix. Each element corresponds to the estimated  degrees of freedom of the resulting remMap estimator under a pair of lamL1 (rows) and lamL2 (columns). 
       }

\references{
 J. Peng, J. Zhu, A. Bergamaschi, W. Han, D.-Y. Noh, J. R. Pollack, P. Wang, 
       Regularized Multivariate Regression for Identifying Master Predictors with Application to Integrative Genomics Study of Breast Cancer.
  (http://arxiv.org/abs/0812.3671)
}

\author{Jie Peng, Pei Wang, Ji Zhu}

\keyword{methods}

\examples{


############################################
############# Generate an example data set
############################################
n=50
p=30
q=30
set.seed(1)

## generate X matrix
X.m<-matrix(rnorm(n*p),n,p)

## generate coefficient
coef.m<-matrix(0,p,q)
hub.n=10
hub.index=sample(1:p, hub.n)
for(i in 1:q){
  cur=sample(1:3,1)
  temp=sample(hub.index, cur)
  coef.m[temp,i]<-runif(length(temp), min=2, max=3)
 }

## generate responses
E.m<-matrix(rnorm(n*q),n,q)
Y.m<-X.m%*%coef.m+E.m

##############################################
# 1. ## fit model for one pair of (lamL1, lamL2)
##############################################

try1=remMap(X.m, Y.m,lamL1=10, lamL2=5, phi0=NULL, C.m=NULL)

################################################################################################
# 2. ## Select tuning parameters with BIC:
##   ## computationally easy; but the BIC procedure assumes orthogonality of the design matrix
##   ## to estimate the degrees of freedom;
##   ## thus it tends to select too small models when the actual design matrix (X.m) is far
##   ##  from orthogonal
################################################################################################

lamL1.v=exp(seq(log(10),log(20), length=3))
lamL2.v=seq(0,5, length=3)
df.m=remMap.df(X.m, Y.m, lamL1.v, lamL2.v, C.m=NULL)
##  The estimated degrees of freedom can be used to select the ranges of tuning parameters.

try2=remMap.BIC(X.m, Y.m,lamL1.v, lamL2.v, C.m=NULL)
pick=which.min(as.vector(t(try2$BIC)))
result=try2$phi[[pick]]
FP=sum(result$phi!=0 & coef.m==0) ## number of false positives
FN=sum(result$phi==0 & coef.m!=0) ## number of false negatives

#BIC selected tuning parameters
print(paste("lamL1=", round(result$lam1,3), "; lamL2=", round(result$lam2,3), sep=""))
print(paste("FP=", FP, "; FN=", FN, sep=""))

#################################################################################################
# 3. ## Select tuning parameters with v-fold cross-validation;
##   ## computationally demanding;
##   ## but cross-validation assumes less assumptions than BIC and thus is recommended unless
##   ## computation is a concern;
#    ## alos cv based on unshrinked estimator (ols.cv) is recommended over cv based on shrinked
##   ## estimator (rss.cv);
##   ## the latter tends to select too large models.
################################################################################################

lamL1.v=exp(seq(log(10),log(20), length=3))
lamL2.v=seq(0,5, length=3)
try3=remMap.CV(X=X.m, Y=Y.m,lamL1.v, lamL2.v, C.m=NULL, fold=5, seed=1)

############ use CV based on unshrinked estimator (ols.cv)
pick=which.min(as.vector(try3$ols.cv))
#pick=which.min(as.vector(try3$rss.cv))
lamL1.pick=try3$l.index[1,pick]    ##find the optimal (LamL1,LamL2) based on the cv score
lamL2.pick=try3$l.index[2,pick]
##fit the remMap model under the optimal (LamL1,LamL2).
result=remMap(X.m, Y.m,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
FP=sum(result$phi!=0 & coef.m==0) ## number of false positives
FN=sum(result$phi==0 & coef.m!=0) ## number of false negatives
##CV (unshrinked) selected tuning parameters
print(paste("lamL1=", round(lamL1.pick,3), "; lamL2=", round(lamL2.pick,3), sep=""))
print(paste("FP=", FP, "; FN=", FN, sep=""))

}
