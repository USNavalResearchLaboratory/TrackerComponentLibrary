%%DEMOCONJUGATEPRIORESTIMATION Conjugate prior distributions provide a
%              simple method of performing recursive Bayesian estimation of
%              unknown parameters of probability density functions based on
%              samples. This file demonstrates how the statistical
%              functions in the Tracker Component Library to perform such
%              estimation.
%
%Each of the examples initialize the conjugate prior distirbutions with
%unifnromative priors. Then, they are updated given a set of measurements
%and the estimates are shown. The estimates can be compared to the original
%variables being estimated in the code. 
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of samples to use for the estimation problems.
numSamp=5000;
display('DISTRIBUTION I: GAMMA')
display('Estimating the scale parameter of a central gamma distribution having known shape parameter.')
display('The conjugate prior distribution is an inverse gamma distribution.')

kMeas=5;%The shape parameter of the gamma distribution.
thetaMeas=10;%The scale parameter of the gamma distribution.

display('Generating random samples; the shape parameter is k=5, the scale parameter to be estimated is theta=10.')
%Generate the random samples.
xMeas=GammaD.rand([numSamp,1],kMeas,thetaMeas);

display('An uninformative conjugate prior is used; it has k=0 and beta=0.')
%The prior should be uninformative. In this case, that means k=0 and
%beta=0.
kEst=0;
betaEst=0;
display('Updating the uninformative conjugate prior using the measurements.')
[kEst,betaEst]=gammaScaleLikeInvGammaConjUpdate(xMeas,kMeas,kEst,betaEst);

display('The mean of the estimated scale parameter is')
thetaEst=InverseGammaD.mean(kEst,betaEst)
display('The variance of the estimate is')
thetaVar=InverseGammaD.var(kEst,betaEst)

display('Estimating the shape parameter of a central gamma distribution having known scale parameter.')
display('The conjugate prior distribution is a gamma conjugate type I distribution.')

display('An uninformative conjugate prior is used; it has delta=0 and mu=1.')
deltaEst=0;
muEst=1;

display('Evaluating the mean of the gamma conjugate type I distribution with a lot of samples')
display('is problematic due to finite precision problems. However, since the estimates')
display('obtained using disjoint sets of samples are independent, we can compute')
display('independent estimates from small groups of samples and then average them, lessening the')
display('precision issues.')
kEst=0;

%We are using 40 groups, so for 5000 samples, these are 125 measurements
%apiece. We want to use the largest groups without overflow, because all of
%the numerical integration also does not guarantee that the results are
%unbiased.
numGroups=40;
numMeasPerGroup=numSamp/numGroups;
minIdx=1;
for curMeas=1:numGroups
    measSet=minIdx:(minIdx+numMeasPerGroup-1);
    
    %Each time in the loop, an uninformative prior is used.
    [muCur,deltaCur]=gammaShapeLikeGammaConjIConjUpdate(xMeas(measSet),thetaMeas,muEst,deltaEst);
    kEst=kEst+GammaConjugateID.mean(muCur,deltaCur,[],'RelTol',1e-12,'AbsTol',1e-15);
    minIdx=minIdx+numMeasPerGroup;
end
display('The mean of the estimated shape parameter is')
kEst=kEst/numGroups

display('DISTRIBUTION II: POISSON')
display('Estimating the rate parameter of a Poisson distribution.')
display('The conjugate prior distribution is a gamma distribution.')

lambdaMeas=50;
xMeas=PoissonD.rand([numSamp,1],lambdaMeas);

display('An uninformative conjugate prior is used; it has k=0 and theta=Inf.')
%The prior should be uninformative. In this instance, this means that k=0
%and theta=Inf.
kEst=0;
thetaEst=Inf;

display('Updating the uninformative conjugate prior using the measurements.')
[kEst,thetaEst]=PoissonRateLikeGammaConjUpdate(xMeas,kEst,thetaEst);
%Obtain the mean and variance of the estimate
display('The mean of the estimated rate parameter is')
lambdaEst=GammaD.mean(kEst,thetaEst)
display('The variance of the estimated rate parameter is')
lambdaVar=GammaD.var(kEst,thetaEst)


display('DISTRIBUTION III: BINOMIAL')

nMeas=30;
pMeas=0.3;
xMeas=BinomialD.rand([numSamp,1],nMeas,pMeas);

%The prior should be uninformative. This means that a=0 and b=0.
aEst=0;
bEst=0;

[aEst,bEst]=binomialProbLikeBetaConjUpdate(xMeas,nMeas,aEst,bEst);
pEst=BetaD.mean(aEst,bEst)
pEstVar=BetaD.var(aEst,bEst)

display('DISTRIBUTION IV: EXPONENTIAL')

lambdaMeas=10;
xMeas=ExponentialD.rand([numSamp,1],lambdaMeas);

%The prior should be uninformative. In this instance, this means that k=0
%and theta=Inf.
kEst=0;
thetaEst=Inf;

[kEst,thetaEst]=expRateLikeGammaConjUpdate(xMeas,kEst,thetaEst);
%Obtain the mean and variance of the estimate
lambdaEst=GammaD.mean(kEst,thetaEst)
lambdaVar=GammaD.var(kEst,thetaEst)

display('DISTRIBUTION V: PARETO TYPE I')

xMinMeas=5;
aMeas=4;
xMeas=ParetoTypeID.rand([numSamp,1],xMinMeas,aMeas);

%The prior should be uninformative. In this instance, this means that k=0
%and theta=Inf.
kEst=0;
thetaEst=Inf;

[kEst,thetaEst]=ParetoIShapeLikeGammaConjUpdate(xMeas,xMinMeas,kEst,thetaEst);

%Obtain the mean and variance of the estimate
xMinMeasEst=GammaD.mean(kEst,thetaEst)
xMinMeasVar=GammaD.var(kEst,thetaEst)

display('DISTRIBUTION VI: SCALAR UNIFORM')

UBMeas=12;
xMeas=UniformD.rand(numSamp,[0;UBMeas]);

%The prior should be uninformative. In this instance, this means that
%xMin=0 and a=0;
xMinEst=0;
aEst=0;

[xMinEst,aEst]=scalUniformBoundLikeParetoIConjUpdate(xMeas,xMinEst,aEst);
UBEst=ParetoTypeID.mean(xMinEst,aEst)
UBVar=ParetoTypeID.var(xMinEst,aEst)

display('DISTRIBUTION VII: MULTIVARIATE GAUSSIAN')

muMeas=[1;2;3;4];
SigmaMeas=[62, 7,  12,  17;
            7, 52, 17, 22;
           12, 17, 42, 27;
           17, 22, 27, 32];
SigmaMeasInv=inv(SigmaMeas);
xMeas=GaussianD.rand(numSamp,muMeas,SigmaMeas);

%Estimating the mean for a given covariance. The prior should be
%uninformative.
muEst=zeros(4,1);
SigmaEstInv=zeros(4,4);

[muEst,SigmaEstInv]=GaussMeanLikeGaussConjUpdate(xMeas,SigmaMeasInv,muEst,SigmaEstInv);
%The mean estimate is
muEst=GaussianD.mean(muEst)
muEstCov=GaussianD.cov(inv(SigmaEstInv))

%Estimating the covariance matrix for a given mean. The prior should be
%uninformative.
PsiEst=zeros(4,4);
nuEst=0;

[PsiEst,nuEst]=GaussCovLikeInvWishartConjUpdate(xMeas,muMeas,PsiEst,nuEst);
InverseWishartD.mean(PsiEst,nuEst)

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
