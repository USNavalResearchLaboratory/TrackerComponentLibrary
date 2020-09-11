function [H,AMISE,exitCode]=plugIn2DGaussKernelBW(xi,numStages,preProcAlg,HInit,epsilon,deltaTestDist,delta,lineSearchParams,cholSemiVal,maxIter)
%%PLUGIN2DGAUSSKERNELBW Estimating the optimal bandwidth when approximating
%             a PDF given sampled 2D points using a Gaussian kernel
%             estimator. A kernel estimator puts some type of a kernel
%             around each point to estimate a continuous PDF from the
%             discrete points. This function assumes that a N(0,H*H')
%             Gaussian kernel is used. In this case, the "bandwidth" to be
%             estimated is the optimal lower-triangular covariance matrix H
%             to use. The estimation algorithm used is a so-called
%             "plug-in" estimator.
%
%INPUTS: xi A 2XN vector of N samples of the PDF.
%   numStages The number of stages to use in the estimator. This must be
%           >=1, and is normally not very large. If omitted or an empty
%           matrix is passed, the default value of 2 is used.
%  preProcAlg This specifies how the data is preprocessed before running
%           the algorithm. The data must be preprocessed as the algorithm
%           assumes that the standard deviation in each dimension is the
%           same. After the algorithm is run, the data is transformed back.
%           Possible values are:
%           0 (The default if omitted or an empty matrix is passed) The
%             data is sphered. This means that it is pre-multiplied by
%             the inverse of the sample covariance matrix. The scaled
%             data has the identity matrix as its covariance matrix.
%           1 The data is scaled. This does not eliminate the
%             cross-correlation terms from the data. In some scenarios,
%             this can produce a better estimate.
%     HInit An initial estimate of the bandwidth matrix. If this
%           parameter is omitted or an empty matrix is passed, a standard
%           Gaussian approximation as given in Chapter 3 of [2] is used.
%   epsilon The parameter determining the accuracy of the desired
%           solution in terms of the gradient. The function terminates
%           when norm(g) < epsilon*max([1, norm(x)]) where g is the
%           gradient. The default if omitted or an empty matrix is as
%           given in NewtonsMethod.
% deltaTestDist The number of iterations back to use to compute the
%           decrease of the objective function if a delta-based convergence
%           test is performed. If zero, then no delta-based convergence
%           testing is done. The default is as given in NewtonsMethod.
%     delta The delta for the delta convergence test. This determines
%           the minimum rate of decrease of the objective function.
%           Convergence is determined if (f'-f)<=delta*f, where f' is the
%           value of the objective function f deltaTestDist iterations ago,
%           and f is the current objective function value. The default if
%           this parameter is omitted or an empty matrix is passed is 0.
% lineSearchParams An optional structure whose members specify tolerances
%           for the line search. The parameters are described as in the
%           lineSearch function, except if -1 is passed instead of a
%           parameter structure, then no line search is performed. The
%           default if omitted or an empty matrix is passed is -1.
% cholSemiVal A value indicating whether a method similar to that suggested
%           in Chapter 1.4 of [1] for using a modified Cholesky
%           decomposition of the Hessian matrix should be used so as to
%           assure that one always obtains a descent direction. If a
%           positive value <1 of cholSemiVal is given, then the Hessian
%           matrix is decomposed with cholSemiDef with epsVal of
%           cholSemiVal. If a negative value is provided, then no
%           decomposition is done and no effort is made to assure that the
%           direction traveled is a descent direction. The default if this
%           parameter is omitted or an empty matrix is passed is 1e-6, This
%           parameter assumes the "Hessian" matrix is symmetric. Thus, this
%           option should be set to -1 if this function is used to zero a
%           vector rather than minimize a function.
%   maxIter The maximum number of iterations to perform. If omitted or an
%           empty matrix is passed, the default value of 50 is used.
%
%OUTPUTS: H The 2X2 lower-triangular estimate of the kernel bandwidth.
%     AMISE The estimated AMISE of the given bandwidth.
%  exitCode A value indicating how the function terminated. Possible values
%           are those returned by the NewtonsMethod function. A 0 or
%           positive number means a successful termination.
%
%This file implements the algorithm of [1]. Unlike in [1], where a
%quasi-Newton method is used, the Hessian matrix was explcitly derived
%and Newton's method is used directly. Also, the algorithm of [1] defined
%the "bandwidth" to be the square of what is used as the bandwidth as the
%output of this function, so care must be taken when interpreting the code.
%
%EXAMPLE 1:
%Consider density D in [1]. It is a two-component Gaussian mixture.
% prob1=0.5;%First component is 50% likely.
% mu1=[1;-1];
% Sigma1=[4/9,  14/45;
%         14/45,4/9];
% S1=chol(Sigma1,'lower');
% mu2=[-1;1];
% Sigma2=[4/9, 0;
%         0,   4/9];
% S2=chol(Sigma2,'lower');
% %Generate the random samples.
% numSamples=1000;
% xi=zeros(2,numSamples);
% for curSample=1:numSamples
%     if(rand(1)<prob1)
%         xi(:,curSample)=mu1+S1*randn(2,1);
%     else
%         xi(:,curSample)=mu2+S2*randn(2,1);
%     end
% end
% %Given the random samples, find a bandwidth to use for a kernel-based
% %estimator
% [H,AMISE,exitCode]=plugIn2DGaussKernelBW(xi,[],0);
% %Plot the true PDF.
% numPoints=250;
% vals=linspace(-3,3,numPoints);
% [X,Y]=meshgrid(vals,vals);
% points=[X(:)';Y(:)'];
% PDFVals=prob1*GaussianD.PDF(points,mu1,Sigma1)+(1-prob1)*GaussianD.PDF(points,mu2,Sigma2);
% PDFVals=reshape(PDFVals,numPoints,numPoints);
% figure()
% surf(X,Y,PDFVals,'EdgeColor','none')
% title('True PDF')
% %Plot the kernel estimation of the PDF.
% PDFEst=kernelApprox(points,xi,H);
% PDFEst=reshape(PDFEst,numPoints,numPoints);
% figure()
% surf(X,Y,PDFEst,'EdgeColor','none')
% title('Approximated PDF')
% 
%EXAMPLE 2:
%Consider density F in [1]. It is a single Gaussian.
% mu=[0;0];
% Sigma=[1,   0.9;
%        0.9, 1];
% S=chol(Sigma,'lower');
% %Generate the random samples.
% numSamples=500;
% xi=S*randn(2,numSamples);
% %Given the random samples, find a bandwidth to use for a kernel-based
% %estimator
% [H,AMISE,exitCode]=plugIn2DGaussKernelBW(xi,[],0);
% %Plot the true PDF.
% numPoints=250;
% vals=linspace(-3,3,numPoints);
% [X,Y]=meshgrid(vals,vals);
% points=[X(:)';Y(:)'];
% PDFVals=GaussianD.PDF(points,mu,Sigma);
% PDFVals=reshape(PDFVals,numPoints,numPoints);
% figure()
% surf(X,Y,PDFVals,'EdgeColor','none')
% title('True PDF')
% %Plot the kernel estimation of the PDF.
% PDFEst=kernelApprox(points,xi,H);
% PDFEst=reshape(PDFEst,numPoints,numPoints);
% figure()
% surf(X,Y,PDFEst,'EdgeColor','none')
% title('Approximated PDF')
% 
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.
%[2] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(numStages))
    numStages=2;
end

if(nargin<3||isempty(preProcAlg))
    preProcAlg=0;
end

if(nargin<5||isempty(epsilon))
    epsilon=[];
end

if(nargin<6||isempty(deltaTestDist))
    deltaTestDist=[];
end

if(nargin<7||isempty(delta))
    delta=[];
end

%Take any line search parameters given.
if(nargin<8||isempty(lineSearchParams))
    lineSearchParams=-1;%No line search.
end

if(nargin<9||isempty(cholSemiVal))
    cholSemiVal=1e-6;
elseif(cholSemiVal>1)
    error('cholSemiVal should be less than 1.')
end

if(nargin<10||isempty(maxIter))
    maxIter=50; 
end

n=size(xi,2);
jMax=2*numStages+4;
%This holds all of the values of psi_r that are needed for Equation 5.
psiMat=zeros(jMax+1,jMax+1);

[~,R]=calcMixtureMoments(xi);
switch(preProcAlg)
    case 0%Sphere the data.
        SFull=chol(R,'lower');
        xi=SFull\xi;
        SScale=SFull;
        inv2SampCov=inv(2*eye(2));
    case 1%Scale the data
        SScale=diag(sqrt(diag(R)));
        xi=SScale\xi;
        inv2SampCov=inv(2*SScale\R/SScale);
    otherwise
        error('Unknown algorithm for the data preprocessing step given')
end

%If no initial estimate is given, then use the standard Gaussian
%approximation based on the standard deviation.
if(nargin<4||isempty(HInit))
    %Due to the preprocessing step, the standard deviation value to use is
    %always 1 in each dimension.
    HInitSqrt=(4/3)^(1/5)*eye(2)*n^(-1/5);
    HInit=HInitSqrt*HInitSqrt';
else
    HInit=SScale\HInit;%Account for pre-scaling.
    HInit=HInit*HInit';
end

%Here the steps of Section 4.1 of [1] are implemented.
%%STEP 1
%We must obtain the normal reference estimates.
rVals=genAllTCompositions(jMax+2,2)-1;
numR=size(rVals,2);

for k=1:numR
    r=rVals(:,k);
    psiMat(r(1)+1,r(2)+1)=psiRND(r,inv2SampCov);
end

%Obtain g from Equation 8.
[A2,A3,A4]=getA2A3A4(jMax-2,psiMat);
g=gjSAMSE(jMax-2,n,A2,A3,A4);

%%STEP 2
for j=(jMax-2):-2:4
    %STEP 2, part a
    %Use g and Equation 4 to estimate values of Psi for the current j
    %value.
    rVals=genAllTCompositions(j+2,2)-1;
    numR=size(rVals,2);
    for k=1:numR
        r=rVals(:,k);
        psiMat(r(1)+1,r(2)+1)=psiR(rVals(:,k),xi,g);
    end
    
    %STEP 2, part b
    %Get g from Equation 8 using the psi values for j-2.
    if(j>4)
        [A2,A3,A4]=getA2A3A4(j-2,psiMat);
        g=gjSAMSE(j-2,n,A2,A3,A4);
    end
end

%%STEP 3
%Get the big Psi4 matrix from Equation 3.
Psi4Val=Psi4(psiMat);

%%STEP 4
%Now, numerically minimize the AMISE in Equation 2 for the given Psi4Val
%matrix using Newton's method.

%f returns the function value, the gradient, and the Hessian.
f=@(h)calcAMISE(vech2Mat(h),n,Psi4Val);
fHess=[];
h0=vech(HInit);
[h,AMISE,exitCode]=NewtonsMethod(f,fHess,h0,epsilon,deltaTestDist,delta,lineSearchParams,cholSemiVal,maxIter);
H=vech2Mat(h);

%Finally undo the prescale/sphering step.
H=SScale*H*SScale';

%And, take the square root
H=chol(H,'lower');

end

function [AMISE,gradAMISE,HessAMISE]=calcAMISE(H,n,Psi4Val)
%%CALCAMISE This implements Equation 2 in [1] along with the graidient and
%           Hessian with respect to [H(1,1);H(2,1);H(2,2)]. The gradient
%           and Hessian are not given in the paper and are explained below.
%           The paper must have used the gradient for the quasi-Newton
%           method it used for minimization. We are using the Hessian too
%           so that we can use Newton's method directly.
%
%Equation 2 can be written as
% AMISE=c1/sqrt(det(H))+vec(H)'*C2*vec(H)
%where
% vech(H)=[H(1,1);H(2,1);H(2,2)]
% c1=n^(-2)*R(K);
% C2=(1/4)*mu2(K)^2*Psi4;
%and one could write sqrt(det(H))=sqrt(vec(H)'*A*vec(H))
%where  A=[0,   0, 1/2;
%          0,  -1, 0;
%          1/2, 0, 0];
%The gradient with respect to [H(1,1);H(2,1);H(2,2)] is thus
%grad(AMISE)=-c1*A*vec(H)/sqrt(det(H))^3+2*C2*vec(H)
%To be able to use Newton's method to zero this function, we need the
%Hessian matrix. The Hessian matrix is
%Hess(AMISE)=3*c1*(A*vec(h))*(A*vec(h))'/sqrt(det(H))^5-c1*A/sqrt(det(H))^3+2*C2
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

%These values are for a standard Gaussian kernel.
RK=1/(4*pi);
mu2K=1;

vechH=vech(H);
sqrtDetH=sqrt(det(H));
c1=RK/n;
C2=(1/4)*mu2K^2*Psi4Val;

AMISE=c1/sqrtDetH+vechH'*C2*vechH;

A=[0,   0, 1/2;
   0,  -1, 0;
   1/2, 0, 0];
gradAMISE=-c1*A*vechH/sqrtDetH^3+2*C2*vechH;
HessAMISE=3*c1*(A*vechH)*(A*vechH)'/sqrtDetH^5-c1*A/sqrtDetH^3+2*C2;
end

function val=psiRND(r,inv2SampCov)
%%PSIRND Compute \hat{\psi}_r^{ND} using the formula at the beginning of
%        Section 4 of [1]. inv2S is the inverse of two times the sample
%        covariance matrix of the data.
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

[~,coeffs]=GaussianD.PDFIDerivs(zeros(2,1),inv2SampCov,r);
val=(-1)^(sum(r))*polyValMultiDim(coeffs,zeros(2,1))/(2*pi);
end

function [A2,A3,A4]=getA2A3A4(j,psiMat)
%%GETA2A3A4 Compute the values of A2, A3, and A4 that go into Equation 8 of
%        [1]. j is the order of the sums. psiMat(r(1)+1,r(2)+1) holds the
%        value for the given r of the psi values already found. This means
%        that memory is wated as the sums only use a particular order, but
%        indexation is simpler. The memory wasted does not matter as the
%        maximum number of stages used in the algorithm will generally be
%        small, so j will not be very large.
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

%All of the r values for the sums in the Equations for A2,A3, and A4. The
%values start from 0.
rVals=genAllTCompositions(j+2,2)-1;
numRVals=size(rVals,2);

%The value of mu_2(K) when using the standard bivariate normal distribution
%as the kernel.
mu2=1;

A2=0;
A3=0;
A4=0;
for k=1:numRVals
    r=rVals(:,k);
    %Evaluate K^r(0)
    [~,coeffs]=GaussianD.PDFIDerivs(zeros(2,1),eye(2),r);
    KTerm=polyValMultiDim(coeffs,zeros(2,1))/(2*pi);
   
    %Update the A2 sum.
    A2=A2+KTerm^2;
    
    PsiSum=psiMat(r(1)+2+1,r(2)+1)+psiMat(r(1)+1,r(2)+2+1);
    A3=A3+KTerm*PsiSum;
    A4=A4+PsiSum^2;
end

A3=mu2*A3;
A4=mu2^2*A4;
end

function val=gjSAMSE(j,n,A2,A3,A4)
%This implements Equation 8 of [1].
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

innerTerm=(4*j+8)*A2/(n*(-j*A3+sqrt(j^2*A3^2+(8*j+16)*A2*A4)));
val=nthroot(innerTerm,j+4);
end

function val=psiR(r,xi,g)
%%PSIR Compute \hat{\psi}_r(G) from Equation 4 of [1] given the sample
%      points xi, g, as G=g^2*eye(2) (only a bivariate distribution
%      considered), and r, which specifies the number of derivatives of the
%      kernel (a bivariate normal distribution with covariance matrix G) to
%      use.
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

n=size(xi,2);

%Get the coefficients for the polynomial term of the kernel derivative.
[~,coeffs]=GaussianD.PDFIDerivs(zeros(2,1),1/g^2*eye(2),r);

%Sum the off-diagonal terms.
sumVal=0;
for i=1:(n-1)
    diff=bsxfun(@minus,xi(:,i),xi(:,(i+1):n));
    expVal=exp((-1/2)*sum(diff.*diff,1)/g^2);
    sumVal=sumVal+(2*pi)^(-1)*(1/g^2)*sum(expVal.*polyValMultiDim(coeffs,diff));
end
%The off-diagonal terms are all the same.
sumVal=sumVal+n*polyValMultiDim(coeffs,zeros(2,1));
%The common term.
val=sumVal/n^2;
end

function Psi4Val=Psi4(psiMat)
%%PSI4VAL Given a matrix of the psi_{a,b} terms, create the big Psi4 matrix
%         that is in Equation 3 in [1].
%
%REFERENCES:
%[1] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.

Psi4Val=[psiMat(4+1,0+1),   2*psiMat(3+1,1+1),  psiMat(2+1,2+1);
        2*psiMat(3+1,1+1), 4*psiMat(2+1,2+1),  2*psiMat(1+1,3+1);
        psiMat(2+1,2+1),   2*psiMat(1+1,3+1),  psiMat(0+1,4+1)];

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
