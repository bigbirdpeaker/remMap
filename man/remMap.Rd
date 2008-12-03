\name{remMap}
\alias{remMap}
\title{A function to fit REgularized Multivariate regression model using the MAP penalty (remMap).}
\description{A function to fit regularized multivariate regression model using the MAP penalty.}
   
\usage{
remMap(X.m, Y.m,lamL1, lamL2, phi0=NULL, C.m=NULL)
}

\arguments{
    \item{X.m}{numeric matrix (n by p): columns correspond to predictor variables and rows correspond to samples. 
                     Missing values are not allowed.}  
    \item{Y.m}{numeric matrix (n by q): columns correspond to response variables and rows correspond to samples. 
                     Missing values are not allowed.}                    
    \item{lamL1}{numeric value: \eqn{l_1} norm penalty parameter.} 
    \item{lamL2}{numeric value: \eqn{l_2} norm penalty parameter.}    
    \item{phi0}{numeric matrix (p by q): an initial estimate of the coefficient matrix of the multivariate 
                regression model; default(=NULL): univariate estimates are used as initial estimates.}
    \item{C.m}{numeric matrix (p by q): \eqn{C_m[i,j]=0} means the corresponding coefficient beta[i,j] is set to be zero in the model;
                  \eqn{C_m[i,j]=1} means the corresponding beta[i,j] is included in the MAP penalty; 
                  \eqn{C_m[i,j]=2} means the corresponding beta[i,j] is not included in the MAP penalty; 
                  default(=NULL): \eqn{C_m[i,j]} are all set to be 1.}
          }        
\details{
    \code{remMap} uses a computationally efficient approach for performing multivariate regression 
         under the high-dimension-low-sample-size setting (Peng and et. al., 2008).}

\value{
    A list with two components
    \item{phi}{the estimated coefficient matrix (p by q) of the regularized multivariate regression model.}
    \item{rss.v}{a vector of length q recording the RSS values of the q regressions.}
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
n=100 
p=300 
q=300
set.seed(1)

### generate X matrix
rho=0.5; Sig<-matrix(0,p,p)
for(i in 2:p){ for(j in 1: (i-1)){
    Sig[i,j]<-rho^(abs(i-j))
    Sig[j,i]<-Sig[i,j]
}}
diag(Sig)<-1
R<-chol(Sig)
X.m<-matrix(rnorm(n*p),n,p)
X.m<-X.m\%*\%R

### generate coefficient
coef.m<-matrix(0,p,q)
hub.n=20
hub.index=sample(1:p, hub.n)
for(i in 1:q){
  cur=sample(1:3,1)
  temp=sample(hub.index, cur)
  coef.m[temp,i]<-runif(length(temp), min=2, max=3)
 }

### generate responses
E.m<-matrix(rnorm(n*q),n,q)
Y.m<-X.m\%*\%coef.m+E.m


##############################################################################################
############ perform analysis
##############################################################################################

###############################################
## 1. ## fit model for one pair of (lamL1, lamL2)
###############################################

try1=remMap(X.m, Y.m,lamL1=100, lamL2=50, phi0=NULL, C.m=NULL)

#################################################################################################
## 2. ## Select tuning parameters with BIC: 
###   ## computationally easy; but the BIC procedure assumes orthogonality of the design matrix to estimate the degrees of freedom; 
###   ## thus it tends to select too small models when the actual design matrix (X.m) is far from orthogonal
#################################################################################################

lamL1.v=exp(seq(log(51),log(150), length=5))
lamL2.v=seq(0,100, length=5)
df.m=remMap.df(X.m, Y.m, lamL1.v, lamL2.v, C.m=NULL) 
###  The estimated degrees of freedom can be used to select the ranges of tuning parameters.

try2=remMap.BIC(X.m, Y.m,lamL1.v, lamL2.v, C.m=NULL)
pick=which.min(as.vector(t(try2$BIC)))
result=try2$phi[[pick]]
FP=sum(result$phi!=0 & coef.m==0) ## number of false positives
FN=sum(result$phi==0 & coef.m!=0) ## number of false negatives


print(paste("lamL1=", round(result$lam1,3), "; lamL2=", round(result$lam2,3), sep=""))  ##BIC selected tuning parameters
print(paste("FP=", FP, "; FN=", FN, sep="")) 


################################################################################################################
## 3. ## Select tuning parameters with v-fold cross-validation;
###   ## computationally demanding; 
###   ## but cross-validation assumes less assumptions than BIC and thus is recommended unless computation is a concern; 
##    ## alos cv based on unshrinked estimator (ols.cv) is recommended over cv based on shrinked estimator (rss.cv); 
###   ## the latter tends to select too large models.
################################################################################################################

lamL1.v=exp(seq(log(51),log(150), length=5))
lamL2.v=seq(0,100, length=5)
try3=remMap.CV(X=X.m, Y=Y.m,lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)

############# use CV based on unshrinked estimator (ols.cv)
pick=which.min(as.vector(try3$ols.cv))
lamL1.pick=try3$l.index[1,pick]    ##find the optimal (LamL1,LamL2) based on the cv score
lamL2.pick=try3$l.index[2,pick]
result=remMap(X.m, Y.m,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)  ##fit the remMap model under the optimal (LamL1,LamL2).
FP=sum(result$phi!=0 & coef.m==0) ## number of false positives
FN=sum(result$phi==0 & coef.m!=0) ## number of false negatives
print(paste("lamL1=", round(lamL1.pick,3), "; lamL2=", round(lamL2.pick,3), sep="")) ##CV (unshrinked) selected tuning parameters
print(paste("FP=", FP, "; FN=", FN, sep=""))

############# use CV based on shrinked estimator (rss.cv); it tends to select very large models: thus is not recommended in general
pick=which.min(as.vector(try3$rss.cv))
lamL1.pick=try3$l.index[1,pick]    ##find the optimal (LamL1,LamL2) based on the cv score
lamL2.pick=try3$l.index[2,pick]
result=remMap(X.m, Y.m,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
FP=sum(result$phi!=0 & coef.m==0) ## number of false positives
FN=sum(result$phi==0 & coef.m!=0)  ## number of false negatives
print(paste("lamL1=", round(lamL1.pick,3), "; lamL2=", round(lamL2.pick,3), sep="")) ##CV (shrinked) selected tuning parameters
print(paste("FP=", FP, "; FN=", FN, sep=""))


}
