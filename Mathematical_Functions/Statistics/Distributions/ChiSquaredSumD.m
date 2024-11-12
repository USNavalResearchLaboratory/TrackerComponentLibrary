classdef ChiSquaredSumD
%%CHISQUAREDSUMD These are function to handle the generalized chi squared
%      distribution. This is equivalent to the sum of N weighted noncentral
%      chi squared random variables plus an additional weighted standard
%      normal random variable plus a constant. Basically,
%      z=sum_{i=1}^N w(i)*chi2(k(i),lambda(i))+s*z+m
%      where w(i) are weights, k(i) and lambda(i) are the number of degrees
%      of freedom and the noncentrality parameter of the ith distribution,
%      and z is a standard normal random variable.
%
%Implemented methods are: normalQuadForm2ChiSquaredSum,
%                         chiSquaredSum2NormalQuadForm, mean, var, PDF,
%                         CDF, invCDF, and rand.
%
%Such distributions are discussed in [1] and [2]. The integration method
%used for the PDF and the CDF is that of [3], which is also slightly
%generalized in [2]. If one needs higher precision at certain extreme
%values, then the other algorithms in [1] might be better, but are not
%implemented here.
%
%REFERENCES:
%[1] A. Das, "New methods to compute the generalized chi-square
%    distribution," ArXiv, 31 Jul. 2024. [Online]. Available:
%    https://arxiv.org/pdf/2404.05062
%[2] A. Das and W. S. Geisler, "Methods to integrate multinormals and
%    compute classification measures," A Method to Integrate and Classify
%    Normal Distributions, vol. 21, no. 1, Sep. 2021, a corrected version
%    of on ArXiv. [Online]. Available: https://arxiv.org/pdf/2012.14331v11
%[3] J. P. Imhof, "Computing the distribution of quadratic forms in normal
%    variables," Biometrika, vol. 48, no. 3-4, pp. 419-426, Dec. 1961.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function [w,k,lambda,s,m]=normalQuadForm2ChiSquaredSum(mu,Sigma,Q2,q1,q0)
%%NORMALQUADFORM2CHISQUAREDSUM Consider a random variable qx defined in
%       terms of the quadratic form of a Gaussian vector x with mean mu and
%       covariance matrix Sigma: qx=x'*Q2*x+q1'*x+q0. This function takes
%       this form and rewrites it in terms of a weighted sum of noncentral
%       chi squared random variables plus an extra standard normal random
%       variable and also an additive constant
%       qx=sum_{i=1}^N chi2_{k(i),lambda(i)}+s*z+m
%       where chi2_{a,b} represents a noncentral chi squared distribution
%       with a degrees of freedom and a noncentrality parameter b and z is
%       a scalar standard normal (mean 0 variance 1) random variable.
%       This is the opposite of the chiSquaredSum2NormalQuadForm method.
%
%INPUTS: mu The xDimX1 mean value of x.
%     Sigma The xDimXxDim covariance matrix of x.
%        Q2 A real xDimXxDim matrix.
%        q1 A real xDimX1 vector.
%        q0 A real scalar quantity.
%
%OUTPUTS: w The 1XwDim set of weights of the chi-squared terms in the sum
%           of chi squared random variables.
%         k The 1XwDim set of degrees of freedom of each value in the sum
%           of chi squared random variables.
%    lambda The 1XwDim set of noncentrality parameters of the chi squared
%           random variables.
%         s The scalar coefficient of the normal random variable.
%         m The scalar additive term.
%
%The formulae from Section 2.1 of [1] are implemented for the conversion.
%
%EXAMPLE:
%A quadratic normal form is converted into a sum of chi squared random
%variables. Then, random samples of the normal distribution are generated
%and transformed according to the quadratic formula and random samples of
%the chi-squared form are generated and it is shown that the histograms of
%the samples agree.
% numMCRuns=1e4;
% Sigma=[132,   7,  12,  17;
%          7, 122,  17,  22;
%         12,  17, 112,  27;
%         17,  22,  27, 102];
% mu=[1;2;3;4];
% Q2=[-5, -12,  -3,   8;
%     18,  -9,   7,  -7;
%    -10,   3,  -5,  -1;
%      6,  -3,  17,   3];
% q1=[-6;6;8;-5];
% q0=12;
% [w,k,lambda,s,m]=ChiSquaredSumD.normalQuadForm2ChiSquaredSum(mu,Sigma,Q2,q1,q0);
% 
% qx=zeros(1,numMCRuns);
% qChi2=zeros(1,numMCRuns);
% for i=1:numMCRuns
%     x=GaussianD.rand(1,mu,Sigma);
%     qx(i)=x'*Q2*x+x'*q1+q0;
%     qChi2(i)=ChiSquaredSumD.rand(1,w,k,lambda,s,m);
% end
% 
% figure(1)
% clf
% histogram(qx,'normalization','pdf');
% hold on
% histogram(qChi2,'normalization','pdf');
%
%REFERENCES:
%[1] A. Das and W. S. Geisler, "Methods to integrate multinormals and
%    compute classification measures," A Method to Integrate and Classify
%    Normal Distributions, vol. 21, no. 1, Sep. 2021, a corrected version
%    of on ArXiv. [Online]. Available: https://arxiv.org/pdf/2012.14331v11
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(q0))
    q0=0;
end

%Force Q2 to be symmetric. This doesn't change the equation.
Q2s=(1/2)*(Q2+Q2');
%Next, as in Equation 2 of [1], get new Q2, q1 and q0 values that are with
%respect to a standard normal distribution.
S=sqrtm(Sigma);
Q2=S*Q2s*S;
q1New=S*(2*Q2s*mu+q1);
q0=mu'*Q2s*mu+q1'*mu+q0;
q1=q1New;

%Get the parameters shown in the equations after Equation 2.
[R,D]=eig(Q2);
d=diag(D)';
b=(R'*q1)';

%W is a 1XnumUniqueEig vector holding the unique eigenvalues in d.
[w,~,ic]=uniquetol(nonzeros(d)',1e-12); 
numUniqueEig=length(w);
%k is a numUniqueEigX1 vector indicating how many times each unique
%eigenvalue of Q2 is repeated.
k=accumarray(ic,1)';
lambda=zeros(1,numUniqueEig);
for i=1:numUniqueEig
    %The last equation on page 2 of [1]. The noncentrality parameters.
    lambda(i)=sum((b(d==w(i))).^2)/(4*w(i).^2);
end
m=q0-sum(w(:).*lambda(:));
s=norm(b(~d));

end

function [Q2,q1,q0]=chiSquaredSum2NormalQuadForm(w,k,lambda,s,m)
%%CHISQUAREDSUM2NORMALQUADFORM Given a weighted sum of noncentral chi
%       squared random variables plus an extra standard normal random
%       variable and also an additive constant
%       qx=sum_{i=1}^N chi2_{k(i),lambda(i)}+s*z+m 
%       This function find an equivalently distributed expression in terms
%       of the quadratic form of a standard normal  Gaussian vector x (mean
%       0 and the covariance matrix is the identity matrix). This is the
%       opposite of the normalQuadForm2ChiSquaredSum method.
%
%INPUTS: w The 1XwDim or wDimX1 set of weights of the chi-squared terms in
%          the sum of chi squared random variables.
%        k The 1XwDim or wDimX1 set of degrees of freedom of each value in
%          the sum of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
%
%INPUTS: Q2 A real xDimXxDim matrix.
%        q1 A real xDimX1 vector.
%        q0 A real scalar quantity.
%
%The transformation scheme is described in Section 2 of [1].
%
%EXAMPLE:
%A sum of chi squared random variables with an additional additive normal
%random variable is converted into a quadratic form of normal random
%varibales. Then, random samples of the normal distribution are generated
%and transformed according to the quadratic formula and random samples of
%the chi-squared form are generated and it is shown that the histograms of
%the samples agree.
% numMCRuns=1e4;
% w=[1,3,-2,4];
% k=[3,12,1,1];
% lambda=[0,12,8,3];
% s=1;
% m=14;
% [Q2,q1,q0]=ChiSquaredSumD.chiSquaredSum2NormalQuadForm(w,k,lambda,s,m);
% xDim=length(q1);
% mu=zeros(xDim,1);
% Sigma=eye(xDim,xDim);
% 
% qx=zeros(1,numMCRuns);
% qChi2=zeros(1,numMCRuns);
% for i=1:numMCRuns
%     x=GaussianD.rand(1,mu,Sigma);
%     qx(i)=x'*Q2*x+x'*q1+q0;
%     qChi2(i)=ChiSquaredSumD.rand(1,w,k,lambda,s,m);
% end
% 
% figure(1)
% clf
% histogram(qx,'normalization','pdf');
% hold on
% histogram(qChi2,'normalization','pdf');
%
%REFERENCES:
%[1] A. Das, "New methods to compute the generalized chi-square
%    distribution," ArXiv, 31 Jul. 2024. [Online]. Available:
%    https://arxiv.org/pdf/2404.05062
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<5||isempty(m))
        m=0;
    end

    if(nargin<4||isempty(s))
        s=0;
    end

    if(nargin<3||isempty(lambda))
        numW=length(w);
        lambda=zeros(numW,1);
    end

    q0=sum(w(:).*lambda(:))+m;

    hasS=nargin>3&&~isempty(s)&&(s~=0);
    numW=length(w);
    kSum=sum(k);
    numTerms=kSum+hasS;
    q2Diag=zeros(numTerms,1);
    q1=zeros(numTerms,1);
    startIdx=1;
    for i=1:numW
        %w(i) is repeated k(i) times.
        q2Diag(startIdx:(startIdx+k(i)-1))=w(i);
        q1(startIdx)=-2*w(i)*sqrt(lambda(i));
        startIdx=startIdx+k(i);
    end
    Q2=diag(q2Diag);

    if(hasS)
        q1(end)=s;
    end
end

function val=mean(w,k,lambda,m)
%MEAN Find the mean of a sum of weighted chi squared random variables plus
%     a normal random variable.
%
%INPUTS: w The 1XwDim set of weights of the chi-squared terms in the sum
%           of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
%
%OUTPUTS: val The scalar mean of the distribution.
%
%Note that this does not depend on the scale coefficient of the normal
%random variable and thus that coefficient is not requested as an input.
%
%EXAMPLE:
%Here, we show that the algorithm is consistent with the sample mean.
%Typically here to more than 3 digits.
% numMCRuns=1e6;
% w=[1,3,2,4];
% k=[3,12,1,1];
% lambda=[0,12,8,3];
% s=1;
% m=14;
% vals=ChiSquaredSumD.rand([1,numMCRuns],w,k,lambda,s,m);
% meanVal=ChiSquaredSumD.mean(w,k,lambda,m);
% sampleMean=mean(vals);
% RelDiff=(sampleMean-meanVal)/meanVal
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(m))
    m=0;
end

if(nargin<3||isempty(lambda))
    numW=length(w);
    lambda=zeros(numW,1);
end

val=sum(w(:).*(k(:)+lambda(:)))+m;
end

function val=var(w,k,lambda,s)
%%VAR Find the variance of a sum of weighted chi squared random variables
%     plus a normal random variable.
%
%INPUTS: w The 1XwDim set of weights of the chi-squared terms in the sum
%           of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%
%OUTPUTS: val The scalar variance of the distribution.
%
%The variance doesn't depend on the scalar additive term m and thus m is
%not requested as an input.
%
%EXAMPLE:
%Here, we show that the algorithm is consistent with the sample mean.
%Typically here to more than 3 digits.
% numMCRuns=1e6;
% w=[1,3,2,4];
% k=[3,12,1,1];
% lambda=[0,12,8,3];
% s=3;
% m=14;
% vals=ChiSquaredSumD.rand([1,numMCRuns],w,k,lambda,s,m);
% varVal=ChiSquaredSumD.var(w,k,lambda,s);
% sampleVar=var(vals);
% RelDiff=(sampleVar-varVal)/varVal
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(s))
    s=0;
end

if(nargin<3||isempty(lambda))
    numW=length(w);
    lambda=zeros(numW,1);
end

val=2*sum((w(:).^2).*(k(:)+2*lambda(:)))+s^2;

end

function vals=PDF(x,w,k,lambda,s,m,AbsTol,RelTol)
%%PDF Evaluate the probability distribution function (PDF) of the
% distribution of a sum of weighted chi squared random variables plus a
% normal random variable (a generalized chi squared distribution).
%
%INPUTS: x A matrix of one or more points at which the PDF should be
%          evaluated.
%        w The 1XwDim or wDimX1 set of weights of the chi-squared terms in
%          the sum of chi squared random variables.
%        k The 1XwDim or wDimX1 set of degrees of freedom of each value in
%          the sum of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
% AbsTol, RelTol Absolute and relative tolerances to use in the integral
%          function. The defaults if omitted or empty matrices are passed
%          are 1e-11 and 1e-8.
%
%OUTPUTS: vals A matrix the same size as x holding values of the PDF
%              evaluated at x.
%
%The general idea is that given in [1]. The form of the solution used here
%is that described in Section 4.1. of [2]. In [1], the additive constant
%and the term for the normal PDF/CDF are missing.
%
%EXAMPLE:
%Here, the PDF is validate by showing that it agrees with a histogram
%generated from random samples.
% numMCRuns=1e5;
% numPts=500;
% w=[1,3,2,4];
% k=[3,12,1,1];
% lambda=[0,12,8,3];
% s=3;
% m=14;
% vals=ChiSquaredSumD.rand([1,numMCRuns],w,k,lambda,s,m);
% x=linspace(0,300,numPts);
% y=ChiSquaredSumD.PDF(x,w,k,lambda,s,m);
% figure(1)
% clf
% histogram(vals,'normalization','pdf');
% hold on
% plot(x,y,'linewidth',2);
%
%REFERENCES:
%[1] J. P. Imhof, "Computing the distribution of quadratic forms in normal
%    variables," Biometrika, vol. 48, no. 3-4, pp. 419-426, Dec. 1961.
%[2] A. Das, "New methods to compute the generalized chi-square
%    distribution," ArXiv, 31 Jul. 2024. [Online]. Available:
%    https://arxiv.org/pdf/2404.05062
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<8||isempty(RelTol))
    RelTol=1e-8;
end

if(nargin<7||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<6||isempty(m))
    m=0;
end

if(nargin<5||isempty(s))
    s=0;
end

if(nargin<4||isempty(lambda))
    numW=length(w);
    lambda=zeros(numW,1);
end

w=w(:);
k=k(:);
lambda=lambda(:);

%Get rid of this warning.
warnStruct=warning();
warning('off','MATLAB:integral:MaxIntervalCountReached')

vals=zeros(size(x));
numVals=numel(x);
for i=1:numVals
    f=@(u)costFun(u,x(i),w,k,lambda,s,m);
    %The max operation keeps the PDF valid given finite precision
    %limitiations.
    vals(i)=max(0,integral(f,0,Inf,'AbsTol',AbsTol,'RelTol',RelTol)/(2*pi));
end
warning(warnStruct);

function val=costFun(u,x,w,k,lambda,s,m)
    %The equations for theta are rho are given before Equation 3.3 in [1].
    %Note the offset due to m shifting the mean.
    theta=sum(k.*atan(w*u)+(lambda.*(w*u))./(1+w.^2*u.^2),1)/2-u*(x-m)/2;
    %Note that this has an extra product term at the end as compared to
    %[1].The extra term represents the s*z term  (extra additive normal
    %distribution). The characteristic function of the sum of two random
    %variables is the product of the characteristic functions of the
    %variables, hence an extra term.
    rho=prod(((1+w.^2*u.^2).^(k/4)).*exp(((w.^2*u.^2).*lambda)./(2*(1+w.^2*u.^2))),1) .* exp(u.^2*s^2/8);
    %The combination of theta and rho for the integral are given in 3-2 in
    %[1] for the CDF. Section 4.2 of [2] notes that a derivative of the CDF
    %integral expression is needed for the PDF.
    val=cos(theta)./rho;
end
end

function vals=CDF(x,w,k,lambda,s,m,AbsTol,RelTol)
%%CDF Evaluate the cumulative distribution function (CDF) of the
% distribution of a sum of weighted chi squared random variables plus a
% normal random variable (a generalized chi squared distribution).
%
%INPUTS: x A matrix of one or more points at which the CDF should be
%          evaluated.
%        w The 1XwDim or wDimX1 set of weights of the chi-squared terms in
%          the sum of chi squared random variables.
%        k The 1XwDim or wDimX1 set of degrees of freedom of each value in
%          the sum of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
% AbsTol, RelTol Absolute and relative tolerances to use in the integral
%          function. The defaults if omitted or empty matrices are passed
%          are 1e-11 and 1e-8;
%
%OUTPUTS: vals A matrix the same size as x holding values of the PDF
%              evaluated at x.
%
%The general idea is that given in [1]. The form of the solution used here
%is that described in Section 4.1. of [2]. In [1], the additive constant
%and the term for the normal PDF/CDF are missing.
%
%EXAMPLE:
%Here, the CDF is validate by showing that it agrees with the empirical
%distribution of a bunch of random samples.
% numMCRuns=1e5;
% numPts=500;
% w=[1,3,2,4];
% k=[3,12,1,1];
% lambda=[0,12,8,3];
% s=3;
% m=14;
% vals=ChiSquaredSumD.rand([1,numMCRuns],w,k,lambda,s,m);
% x=linspace(0,300,numPts);
% y=ChiSquaredSumD.CDF(x,w,k,lambda,s,m);
% yEmp=EmpiricalD.CDF(x,vals);
% 
% figure(1)
% clf
% plot(x,yEmp,'linewidth',4)
% hold on
% plot(x,y,'linewidth',2)
%
%REFERENCES:
%[1] J. P. Imhof, "Computing the distribution of quadratic forms in normal
%    variables," Biometrika, vol. 48, no. 3-4, pp. 419-426, Dec. 1961.
%[2] A. Das, "New methods to compute the generalized chi-square
%    distribution," ArXiv, 31 Jul. 2024. [Online]. Available:
%    https://arxiv.org/pdf/2404.05062
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<8||isempty(RelTol))
    RelTol=1e-8;
end

if(nargin<7||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<6||isempty(m))
    m=0;
end

if(nargin<5||isempty(s))
    s=0;
end

if(nargin<4||isempty(lambda))
    numW=length(w);
    lambda=zeros(numW,1);
end

w=w(:);
k=k(:);
lambda=lambda(:);

vals=zeros(size(x));
numVals=numel(x);
%Get rid of this warning.
warnStruct=warning();
warning('off','MATLAB:integral:MaxIntervalCountReached')
for i=1:numVals
    f=@(u)costFun(u,x(i),w,k,lambda,s,m);
    vals(i)=(1/2)-integral(f,0,Inf,'AbsTol',AbsTol,'RelTol',RelTol)/pi;
end
warning(warnStruct);

%Clip to the limits, as finite precision can affect things.
vals=min(1,max(0,vals));

function val=costFun(u,x,w,k,lambda,s,m)
    %The equations for theta are rho are given before Equation 3.3 in [1].
    %Note the offset due to m shifting the mean.
    theta=sum(k.*atan(w*u)+(lambda.*(w*u))./(1+w.^2*u.^2),1)/2-u*(x-m)/2;
    %Note that this has an extra product term at the end as compared to
    %[1].The extra term represents the s*z term  (extra additive normal
    %distribution). The characteristic function of the sum of two random
    %variables is the product of the characteristic functions of the
    %variables, hence an extra term.
    rho=prod(((1+w.^2*u.^2).^(k/4)).*exp(((w.^2*u.^2).*lambda)./(2*(1+w.^2*u.^2))),1) .* exp(u.^2*s^2/8);
    %The combination of theta and rho for the integral are given in 3-2 in
    %[1].
    val=sin(theta)./(u.*rho);
end
end

function [x,exitCode]=invCDF(prob,w,k,lambda,s,m,numStdInt,BrentVals,intVals)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
% of the distribution of a sum of weighted chi squared random variables
% plus a normal random variable (a generalized chi squared distribution).
%
%INPUTS: prob A matrix of one or more probabiltiies at which the inverse of
%          the CDF should be evaluated. These are between 0 and 1.
%        w The 1XwDim or wDimX1 set of weights of the chi-squared terms in
%          the sum of chi squared random variables.
%        k The 1XwDim or wDimX1 set of degrees of freedom of each value in
%          the sum of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
% numStdInt The number of standard deviations about the mean of the
%          distribution (in both directions) for which the inverse value is
%          searched. The default if omitted or an empty matrix is passed is
%          7.
% BrentVals An optional structure that can have fields AbsTol and maxIter,
%          which set those values in the BrentRootFind function. The
%          defaults if omitted or an empty matrix is passed are 1e-7 and
%          500.
%  intVals An optional structure that can have field AbsTol and RelTol,
%          which are the absolute and relative tolerance to use in the
%          integral function. The defaults if omitted or empty matrices are
%          passed are 1e-11 and 1e-8.
%
%OUTPUTS: x A matrix the same size as prob holding values of the inverse
%           CDF evaluated at prob.
%  exitCode A matrix having the same size as x indicating how the
%           BrentRootFind terminated (or 0 if an exact solution on the
%           boundary is found. The values are the same as the BrentRootFind
%           function, and defined in the comments to that function, except
%           if the probability fell outside of the initial search region
%           and is neither 0 nor 1, then -1 is set for the exit code.
%
%This function just tries to bracket the desired probability within
%+/-numStdInt standard deviation of the mean. If it is outside of those
%bounds, then those bounds are returned, unless 0 or 1 are passed in which
%case 0 or Inf are returned as they are exact solutions. Otherwise, the
%BrentRootFind function is used to numerically search for the solution. If
%numStdInt is too big, then BrentRootFind will not be able to bracket the
%roots. The default value of numStdInt is probably fine most of the time.
%
%EXAMPLE:
%Here, for a given distribution, we get the probability at a particular
%point and then show that when using invCDF, one can get the same point
%back. Also, we start with a given probability and use invCDF to get a
%point.
% w=[1,3,2,4];
% k=[3,12,1,1];
% lambda=[0,12,80,3];
% s=3;
% m=14;
% xTrue=250;
% prob=ChiSquaredSumD.CDF(xTrue,w,k,lambda,s,m);
% xBack=ChiSquaredSumD.invCDF(prob,w,k,lambda,s,m);
% RelxErr=(xBack-xTrue)/xTrue
% 
% probTrue=0.999;
% x=ChiSquaredSumD.invCDF(probTrue,w,k,lambda,s,m);
% probBack=ChiSquaredSumD.CDF(x,w,k,lambda,s,m);
% RelProbErr=(probBack-probTrue)/probTrue
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

AbsTolInt=1e-11;
RelTolInt=1e-8;
if(nargin>8&&~isempty(intVals))
    if(isfield(intVals,'AbsTol'))
        AbsTolInt=intVals.AbsTol;
    end
    if(isfield(intVals,'RelTol'))
        RelTolInt=intVals.RelTol;
    end
end

AbsTol=1e-7;
maxIter=500;
if(nargin>7&&~isempty(BrentVals))
    if(isfield(BrentVals,'AbsTol'))
        AbsTol=BrentVals.AbsTol;
    end
    if(isfield(BrentVals,'maxIter'))
        AbsTol=BrentVals.maxIter;
    end
end

if(nargin<7||isempty(numStdInt))
    numStdInt=7;
end

if(nargin<6||isempty(m))
    m=0;
end

if(nargin<5||isempty(s))
    s=0;
end

numW=length(w);

if(nargin<4||isempty(lambda))
    lambda=zeros(numW,1);
end

if(numW==1&&lambda==0&&s==0)
    %The special case of a single chi squared random varibale and no normal
    %random variable.
    x=zeros(size(prob));
    exitCode=zeros(size(prob));
    numEl=numel(prob);
    for i=1:numEl
        x(i)=w*2*gammaincinv(prob(i),k/2)+m;
    end
    return
end

meanVal=ChiSquaredSumD.mean(w,k,lambda,m);
sigmaVal=sqrt(ChiSquaredSumD.var(w,k,lambda,s));
xSpan=[max(0,meanVal-numStdInt*sigmaVal);meanVal+numStdInt*sigmaVal];

%These are precomputed and used to check if the probabilities are outside
%of the span given.
minVal=ChiSquaredSumD.CDF(xSpan(1),w,k,lambda,s,m,AbsTolInt,RelTolInt);
maxVal=ChiSquaredSumD.CDF(xSpan(2),w,k,lambda,s,m,AbsTolInt,RelTolInt);

x=zeros(size(prob));
exitCode=zeros(size(prob));
numEl=numel(prob);
for i=1:numEl
    %Check the upper and lower bounds based on the span. If the probability
    %is not bounded, then return an extremum.
    if(prob(i)<minVal)
        if(prob(i)==0)
            exitCode(i)=0;
            x(i)=0;
        else
            exitCode(i)=-1;
            x(i)=minVal;
        end
        return
    elseif(prob(i)>maxVal)
        if(prob(i)==1)
            exitCode(i)=0;
            x(i)=Inf;
        else
            exitCode(i)=-1;
            x(i)=maxVal;
        end
        return
    end

    f=@(x)(prob-ChiSquaredSumD.CDF(x,w,k,lambda,s,m,AbsTolInt,RelTolInt));

    [x(i),exitCode(i)]=BrentRootFind(f,xSpan,AbsTol,maxIter);
end
end

function vals=rand(N,w,k,lambda,s,m)
%%RAND Generate random variables from the the distribution of a sum of
% weighted chi squared random variables plus a normal random variable (a
% generalized chi squared distribution).
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        w The 1XwDim or wDimX1 set of weights of the chi-squared terms in
%          the sum of chi squared random variables.
%        k The 1XwDim or wDimX1 set of degrees of freedom of each value in
%          the sum of chi squared random variables.
%   lambda The 1XwDim or wDimX1 set of noncentrality parameters of the chi
%          squared random variables. If omitted or an empty matrix is
%          passed, a vector of all zeros is used.
%        s The scalar coefficient of the normal random variable. If omitted
%          or an empty matrix is passed, 0 is used.
%        m The scalar additive term. If omitted or an empty matrix is
%          passed, 0 is used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generalized chi-squared random variables.
%
%This function just generates the specified chi squared and normal
%variables and sums them.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<6||isempty(m))
    m=0;
end

if(nargin<5||isempty(s))
    s=0;
end

if(nargin<4||isempty(lambda))
    numW=length(w);
    lambda=zeros(numW,1);
end

if(isscalar(N))
    dims=[N, N];
else
    dims=N;
end

numW=length(w);

vals=zeros(dims);
numEl=prod(dims);
for curVar=1:numEl
    sumVal=m+s*randn(1);
    for i=1:numW
        sumVal=sumVal+w(i)*ChiSquareD.rand(1,k(i),lambda(i));
    end
    vals(curVar)=sumVal;
end
end

end
end

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
