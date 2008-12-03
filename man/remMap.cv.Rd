\name{remMap.CV}
\alias{remMap.CV}
\title{Fit remMap models for a series of tuning parameters and return the corresponding v-fold cross-validation scores.}
\description{Fit remMap models for a series of tuning parameters and return the corresponding v-fold cross-validation scores.
Two types of cross-validation scores are computed:  cv based on unshrinked estimator (ols.cv); and cv based on shrinked estimator (rss.cv); 
ols.cv is  recommended. rss.cv tends to select very large models and thus is not recommended in general (especially for very sparse models).
v-fold CV is computationally demanding. But it assumes less 
assumptions than BIC and thus ols.cv is recommended over BIC unless computation is a concern.
}
   
\usage{
remMap.CV(X, Y,lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)
}

\arguments{
    \item{X}{numeric matrix (n by p): columns correspond to predictor variables and rows correspond to samples. 
                     Missing values are not allowed.}  
    \item{Y}{numeric matrix (n by q): columns correspond to response variables and rows correspond to samples. 
                     Missing values are not allowed.}                   
    \item{lamL1.v}{numeric vector: a set of \eqn{l_1} norm penalty parameters.} 
    \item{lamL2.v}{numeric vector: a set of \eqn{l_2} norm penalty parameters.}    
    \item{C.m}{numeric matrix (p by q).  \eqn{C_m[i,j]=0} means the corresponding coefficient beta[i,j] is set to be zero in the model;
                      \eqn{C_m[i,j]=1} means the corresponding beta[i,j] is included in the MAP penalty; 
                      \eqn{C_m[i,j]=2} means the corresponding beta[i,j] is not included in the MAP penalty; 
                  default(=NULL): \eqn{C_m[i,j]} are all set to be 1.}
    \item{fold}{positive integer. It specifies the cross validation fold number; default=10.}   
    \item{seed}{numeric scalar. It specifies the seed of the random number generator in R for generating cross validation subsets; default=1.}
          }        
\details{
    \code{remMap.CV} is used to perform two-dimensional grid search of the tuning parameters (lamL1.v, lamL2.v) based on v-fold
    cross-validation scores. (Peng and et.al., 2008).}

\value{
     A list with four components
    \item{ols.cv}{a numeric matrix recording the cross validation scores based on unshrinked OLS estimators for each pair of (lamL1, lamL2).}
    \item{rss.cv}{a numeric matrix recording the cross validation scores based on shrinked remMap estimators for each pair of (lamL1, lamL2).}
    \item{phi.cv}{a list recording the fitted remMap coefficients on cross validation training subsets. Each component corresponds to one
    cv fold and it is again a list with each component corresponding to the estimated remMap coefficients for one pair of (lamL1, lamL2) on that training subset. }
    \item{l.index}{numeric matrix with two rows: each column is a pair of (lamL1, lamL2) and the kth column corresponds to the kth cv score in as.vector(ols.cv) and  as.vector(rss.cv).}           
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
