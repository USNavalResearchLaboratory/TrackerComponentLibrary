function [H,costVal,exitCode]=UCVKernelBW(xi,preProcAlg,HInit,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter)
%%UCVKERNELBW Estimate the optimal bandwidth when approximating a PDF given
%             sampled numDim-dimensional points using a Gaussian kernel
%             estimator.  A kernel estimator puts some type of a kernel
%             around each point to estimate a continuous PDF from the
%             discrete points. This function assumes that a N(0,H*H')
%             Gaussian kernel is used. In this case, the "bandwidth" to be
%             estimated is the optimal lower-triangular covariance matrix H
%             to use. The bandwidth estimate is obtained by minimizing the
%             unbiased cross-validation (UCV) criterion using a
%             quasi-Newton method.
%
%INPUTS: xi A numDimXN vector of N samples of the PDF.
%  preProcAlg This specifies how the data is preprocessed before running
%           the algorithm. After the algorithm is run, the data is
%           transformed back. Possible values are:
%           0 The data is sphered. This means that it is pre-multiplied by
%             the inverse of the sample covariance matrix. The scaled data
%             has the identity matrix as its covariance matrix.
%           1 The data is scaled. This does not eliminate the cross-
%             correlation terms from the data. In some scenarios, this can
%             produce a better estimate.
%           2 (The default if omitted or an empty matrix is passed) No
%             preprocessing is performed.
%     HInit An initial estimate of the bandwidth matrix. If this parameter
%           is omitted or an empty matrix is passed, a standard Gaussian
%           approximation as given in Chapter 3 of [2] is used.
%   epsilon The parameter determining the accuracy of the desired solution
%           in terms of the gradient. The function terminates when
%           norm(g) < epsilon*max([1, norm(x)]) where g is the gradient.
%           The default if omitted or an empty matrix is as given in
%           quasiNetwonBFGS.
% deltaTestDist The number of iterations back to use to compute the
%           decrease of the objective function if a delta-based convergence
%           test is performed. If zero, then no delta-based convergence
%           testing is done. The default is as given in quasiNetwonBFGS.
%     delta The delta for the delta convergence test. This determines the
%           minimum rate of decrease of the objective function. Convergence
%           is determined if (f'-f)<=delta*f, where f' is the value of the
%           objective function f deltaTestDist iterations ago, and f is the
%           current objective function value. The default if this parameter
%           is omitted or an empty matrix is passed is  as given in
%           quasiNetwonBFGS.
% lineSearchParams An optional structure whose members specify tolerances
%           for the line search. The parameters are described as in the
%           lineSearch function.
%    scaleD A boolean parameter indicating whether the inverse Hessian
%           estimate used in the quasi-Newton method should be scaled as in
%           Equation 1.201 of Section 1.7 of [3]. The default if omitted or
%           an empty matrix is passed is as given in quasiNetwonBFGS.
%   maxIter The maximum number of iterations to perform. If omitted or an
%           empty matrix is passed, the default value of 50 is used.
%
%OUTPUTS: H The numDimXnumDim lower-triangular estimate of the kernel
%           bandwidth.
%   costVal The value of the unbiased cross-validation (UCV) criterion
%           obtained at the optimal point.
%  exitCode A value indicating the termination condition of the algorithm.
%           Nonnegative values indicate success; negative values indicate
%           some type of failure. Possible values are:
%           0 The algorithm termiated successfully based on the gradient
%             criterion.
%           1 The algorithm terminated successfully based on the accuracy
%             criterion.
%          -1 The maximum number of overall iterations was reached.
%           Other negative values correspond to a failure in lineSearch and
%           correspond to the exitCode returned by the lineSearch function.
%
%This function implements the multivariate unbiased cross validation
%algorithm of [1] using a quasi-Newton method.
%
%The derivation of the gradient of the cost functions requires being able
%to take the derivative of 1/sqrt(det(C*C')) with respect to vech(C). This
%is done for each element of vech(C) using a number of identities. First,
%is the identity for taking the derivative of the determinant of a matrix Y
%with respect to a scalar x.
%d(det(Y))/dx=det(Y)*trace(inv(Y)*dY/dx)
%The second identity, needed to use the chain rule, is the derivative of
%the product of two matrices A and B with respect to a scalar x.
%d(A*B)/dx=(dA/dx)*B+A*dB/dx
%Finally, the standard scalar identity (n~=1) that
%d(x^n)/dt=n*x^(n-1)dx/dt
%Putting it all together, one gets one of the terms used in the
%implementation.
%
%EXAMPLE 1:
%Consider bivariate density D in [2]. It is a two-component Gaussian
%mixture.
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
% numSamples=100;
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
% [H,costVal,exitCode]=UCVKernelBW(xi);
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
%Consider bivariate density F in [2]. It is a single Gaussian.
% mu=[0;0];
% Sigma=[1,   0.9;
%        0.9, 1];
% S=chol(Sigma,'lower');
% %Generate the random samples.
% numSamples=100;
% xi=S*randn(2,numSamples);
% %Given the random samples, find a bandwidth to use for a kernel-based
% %estimator
% [H,AMISE,exitCode]=UCVKernelBW(xi);
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
%[1] T. Duong and M. L. Hazelton, "Cross-validation bandwidth matrices for
%    multivariate kernel density estimation," Scandinavian Journal of
%    Statistics, vol. 32, no. 3, pp. 485-506, Sep. 2005.
%[2] T. Duong and M. L. Hazelton, "Plug-in bandwidth matrices for bivariate
%    kernel density estimation," pp. 17-30, 2003.
%[3] D. P. Bertsekas, Nonlinear Programming, 2nd ed. Belmont, MA: Athena
%    Science, 1999.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(preProcAlg))
    preProcAlg=2;
end

if(nargin<4||isempty(epsilon))
    epsilon=[];
end

if(nargin<5||isempty(deltaTestDist))
    deltaTestDist=[];
end

if(nargin<6||isempty(delta))
    delta=[];
end

if(nargin<7||isempty(lineSearchParams))
    lineSearchParams=[];
end

if(nargin<8||isempty(scaleD))
    scaleD=[];
end

if(nargin<9||isempty(maxIter))
    maxIter=50; 
end

numDim=size(xi,1);
n=size(xi,2);

[~,R]=calcMixtureMoments(xi);
switch(preProcAlg)
    case 0%Sphere the data.
        SFull=chol(R,'lower');
        xi=SFull\xi;
        SScale=SFull;
        
        %The lower-triangular square root of the sample covariance matrix
        %after scaling.
        SMod=eye(numDim,numDim);
    case 1%Scale the data
        SScale=diag(sqrt(diag(R)));
        xi=SScale\xi;
        
        %The lower-triangular square root of the sample covariance matrix
        %after scaling.
        SMod=SScale\chol(R,'lower');
    case 2%No scaling at all.
        SScale=eye(numDim,numDim);
        
        %The lower-triangular square root of the sample covariance matrix
        %after scaling.
        SMod=chol(R,'lower');
    otherwise
        error('Unknown algorithm for the data preprocessing step given')
end

%If no initial estimate is given, then use the standard Gaussian
%approximation based on the standard deviation.
if(nargin<3||isempty(HInit))
    HInit=(4/3)^(1/5)*SMod*n^(-1/5);
else
    HInit=SScale\HInit;%Account for pre-scaling.
end

f=@(h)UCV(h,xi);
h0=vech(HInit);

%Use the true Hessian as the initial Hessian estimate. Do a modified
%Cholesky decomposition to get an inverse Hessian matrix that is guaranteed
%positive definite (In other words, if HessInit is not positive definite,
%HessInvInit will not be a true inverse of it.)
HessInit=UCVHessian(HInit,xi);
numVechEls=size(HessInit,1);
L=cholSemiDef(HessInit,'lower',0,1e-6);
optsLT.LT=true;
optsLT.UT=false;
optsUT.UT=true;
optsUT.LT=false;
HessInvInit=linsolve(L',linsolve(L,eye(numVechEls,numVechEls),optsLT),optsUT);

[h,costVal,exitCode]=quasiNetwonBFGS(f,h0,HessInvInit,epsilon,deltaTestDist,delta,lineSearchParams,scaleD,maxIter);
if(isempty(h))
   H=[];
   return;
end

H=vech2Mat(h,0);

%Finally undo the prescale/sphering step.
H=SScale*H;
end

function [val,grad]=UCV(h,xi)
%This function implements Equation 10 and its gradient in [1]. HSqrt must
%be a lower-triangular square root matrix and h is vech(HSqrt);

numDim=size(xi,1);
n=size(xi,2);

HSqrt=vech2Mat(h,0);

numVechEls=(numDim*(numDim+1)/2);

%The first term in Equation 10. Due to the symmetry of the normal
%distribution, terms for xi(:,i)-xi(:,j) are the same as those for 
%xi(:,j)-xi(:,i), so the sums can be simplified.
val=0;
for i=1:(n-1)
    diff=bsxfun(@minus,xi(:,i),xi(:,(i+1):n));
    val=val+sum(GaussianD.PDFS(diff,zeros(numDim,1),sqrt(2)*HSqrt))-2*sum(GaussianD.PDFS(diff,zeros(numDim,1),HSqrt));
end
val=(2/(n*(n-1)))*val;

%Add in the second term in Equation 10.
val=val+(1/n)*(4*pi)^(-numDim/2)*1/abs(det(HSqrt));

%For the gradient terms, we need the gradient of the multivariate Gaussian
%distributions with covariance matrices H and 2*H with respect to vech(HSqrt).
grad=zeros(numVechEls,1);
for i=1:(n-1)
    diff=bsxfun(@minus,xi(:,i),xi(:,(i+1):n));
    
    [gradG,HDetGrad]=GaussianD.PDFSGradHessVechS(diff,zeros(numDim,1),HSqrt);
    gradG2=GaussianD.PDFSGradHessVechS(diff,zeros(numDim,1),sqrt(2)*HSqrt);
    %gradG2 is the gradient with respect to the elements of sqrt(2)*HSqrt.
    %We need to transform it to be with respect to the elements of HSqrt.
    %This is just an application of the chain rule.
    gradG2=sqrt(2)*sum(gradG2,2);

    grad=grad+gradG2-2*sum(gradG,2);
end
grad=(2/(n*(n-1)))*grad;

%Add in the effects of the second term in Equation 10.
grad=grad+(1/n)*(4*pi)^(-numDim/2)*HDetGrad;
end

function Hess=UCVHessian(HSqrt,xi)
%This function implements the Hessian of Equation 10 in [1]. HSqrt must
%be a lower-triangular square root matrix.

numDim=size(xi,1);
n=size(xi,2);

numVechEls=(numDim*(numDim+1)/2);

%For the gradient terms, we need the gradient of the multivariate Gaussian
%distributions with covariance matrices H and 2*H with respect to vech(HSqrt).
Hess=zeros(numVechEls,numVechEls,1);
for i=1:(n-1)
    diff=bsxfun(@minus,xi(:,i),xi(:,(i+1):n));
    
    [~,~,HessG,CDetHess]=GaussianD.PDFSGradHessVechS(diff,zeros(numDim,1),HSqrt);
    [~,~,HessG2]=GaussianD.PDFSGradHessVechS(diff,zeros(numDim,1),sqrt(2)*HSqrt);
    %gradG2 is the gradient with respect to the elements of sqrt(2)*HSqrt.
    %We need to transform it to be with respect to the elements of HSqrt.
    %This is just an application of the chain rule.
    HessG2=2*sum(HessG2,3);

    Hess=Hess+HessG2-2*sum(HessG,3);
end
Hess=(2/(n*(n-1)))*Hess;

%Add in the effects of the second term in Equation 10.
Hess=Hess+(1/n)*(4*pi)^(-numDim/2)*CDetHess;
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
