function [h,exitCode]=plugIn1DGaussKernalBW(xi,method,hInit,RelTol,AbsTol,maxIter)
%%PLUGIN1DGAUSSKERNELBW Estimate the optimal bandwidth when approximating a
%             PDF given sampled scalar points using a Gaussian kernel
%             estimator. A kernel estimator puts some type of a kernel
%             around each point to estimate a continuous PDF from the
%             discrete points. This function assumes that a N(0,h^2)
%             Gaussian kernel is used. In this case, the "bandwidth" to be
%             estimated is the optimal standard deviation h to use. The
%             estimation algorithms here are so-called "plug-in"
%             estimators.
%
%INPUTS: xi An NX1 or 1XN vector of N scalar samples of the PDF.
%    method A parameter indicating which algorithm should be used to
%           estimate the bandwidth. Possible values are
%           0 (The default if omitted or an empty matrix is passed) Use hS2
%             from Section 5 of [2].
%           1 Use the description from Chapter 3.4.6 of [1]. This
%             algorithm uses the continuous interpolating kernel function
%             to estimate the second derivative of the "true" PDF directly.
%             This method is biased.
%     hInit An initial estimate of h. If this parameter is omitted or an
%           empty matrix is passed, then the default value is obtained
%           using the standardGaussKernelBW function with its default
%           options.
%    RelTol The maximum relative error tolerance to declare convergence,
%           a positive scalar. If omitted or an empty matrix is passed, the
%           default value of 1e-3 is used.
%    AbsTol The absolute error tolerance allowed, a positive scalar. If
%           omitted or an empty matrix is passed, the default value of 1e-6
%           is used.
%   maxIter The maximum number of steps of Newton's method to use. If
%           omitted or an empty matrix is passed, the default of 50 is
%           used.
%
%OUTPUT: h The estimate of h or an empty matrix is a non-finite number or a
%          negative value of h arose during any of the iterations.
% exitCode A value indicating how the function terminated. Possible values
%          are
%          0 The relative or absolute error tolerances were achieved.
%          1 A non-finite number arose during an iteration.
%          2 The maximum number of iterations was reached.
%            
%The basic idea behind the procedure used is described in Chapter 3 of [1],
%specifically Chapter 3.4.6, where the use of Newton's method to find the
%solution is suggested. The algorithm of [2] addresses some bias issues
%that are present in the expression given in [1].
%
%EXAMPLE:
% %A scalar bimodal Gaussian mixture distribution with well-separated
% %components.
% prob1=0.4;%First component is 40% likely.
% mu1=-12;
% mu2=11;
% sigma1=2;
% sigma2=1;
% %Generate the random samples.
% numSamples=1000;
% xi=zeros(1,numSamples);
% for curSample=1:numSamples
%     if(rand(1)<prob1)
%         xi(curSample)=mu1+sigma1*randn();
%     else
%         xi(curSample)=mu2+sigma2*randn();
%     end
% end
% %Given the random samples, find a bandwidth to use for a kernel-based
% %estimator using method 0
% [h0,exitCode]=plugIn1DGaussKernalBW(xi,0);
% %And using method 1
% [h1,exitCode]=plugIn1DGaussKernalBW(xi,1);
% %Plot a histogram of the points, the interpolated PDF using the
% %kernel-based estimators, and the true PDF.
% numPoints=500;
% xVals=linspace(-20,20,numPoints);
% PDFEst0=kernelApprox(xVals,xi,h0);
% PDFEst1=kernelApprox(xVals,xi,h1);
% figure()
% hold on
% histogram(xi,'Normalization','probability','BinWidth',2,'FaceColor',[1,1,0]);
% %Plot the true PDF
% PDFVals=prob1*GaussianD.PDF(xVals,mu1,sigma1^2)+(1-prob1)*GaussianD.PDF(xVals,mu2,sigma2^2);
% plot(xVals,PDFVals,'-b','linewidth',4)
% %Plot the approximations.
% plot(xVals,PDFEst0,'--r','linewidth',4)
% plot(xVals,PDFEst1,'--g','linewidth',4)
% legend('Histogram','True PDF','Kernel-Based Estimate 0', 'Kernel-Based Estimate 1','Location','NorthWest')
%One can see that method 0 tends to have a smooth fit, but undershoots the
%second maximum a bit, whereas method one has a number of local minima and
%maxima near the first maximum and undershoots the second maximum more.
%
%REFERENCES:
%[1] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.
%[2] S. J. Sheather and M. C. Jones, "A reliable data-based bandwidth
%    selection method for kernel density estimation," Journal of the Royal
%    Statistical Society. Series B (Methodological), vol. 53, no. 3, pp.
%    683-690, 1991.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(maxIter))
    maxIter=50;
end

if(nargin<5||isempty(AbsTol))
   AbsTol=1e-6;
end

if(nargin<4||isempty(RelTol))
    RelTol=1e-3;
end

%If no initial bandwidth estimate is given.
if(nargin<3||isempty(hInit))
    hInit=standardGaussKernelBW(xi(:)');
end

if(nargin<2||isempty(method))
    method=0;
end

h=hInit;
switch(method)
    case 0
        %Take one step to get the constant term so that it does not have to
        %be recomputed during the other steps.
        [hEst,dhEstdh,alphaHat2ConstTerm]=h2S(h,xi);
        optFun=@(h)h2S(h,xi,alphaHat2ConstTerm);
        f=hEst-h;
        df=dhEstdh-1;
        hNext=h-f/df;
        if(~isfinite(hNext)||hNext<=0)
            h=[];
            exitCode=1;
            return;
        end
        
        diffMag=abs(hNext-h);
        h=hNext;
        %Test relative error.
        if(diffMag/hNext<RelTol)
            exitCode=0;
            return;
        end
    
        %Test absolute error
        if(diffMag<AbsTol)
            exitCode=0;
            return;
        end
        
        startIter=2;
    case 1
        optFun=@(h)hOpt(h,xi);
        startIter=1;
    otherwise
        error('Unknonw method specified');
end

%Continue the iterations.
for curIter=startIter:maxIter
    [hEst,dhEstdh]=optFun(h);
    f=hEst-h;
    df=dhEstdh-1;
    hNext=h-f/df;
    
    if(~isfinite(hNext)||hNext<=0)
        h=[];
        exitCode=1;
        return;
    end
    
    diffMag=abs(hNext-h);
    h=hNext;
    
    %Test relative error.
    if(diffMag/hNext<RelTol)
        exitCode=0;
        return;
    end
    
    %Test absolute error
    if(diffMag<AbsTol)
        exitCode=0;
        return;
    end
end

%Maximum number of iterations reached without convergence.
exitCode=2;

end

function [hEst,dhEstdh]=hOpt(h,xi)
%%HOPT This function find hOpt based on Equation 3.45 in Chapter 3.4.6 of
%      [1] for a Gaussian kernel. The integral over the second derivative
%      is approximated using the second derivative of the scalar Gaussian
%      kernel estimator itself for a given value of h.
%
%INPUTS: h The value of h going into the second derivative approximation
%          for beta (from after Equation 3.45 in [1].)
%       xi The scalar sample points.
%
%OUTPUTS:      hEst The value of the left-most term of Equation 12 in [1].
%           dhEstdh The derivative of hEst woith respect to h.
%
%REFERENCES:
%[1] B. W. Silverman, Density Estimation for Statistics and Data Analysis.
%    Chapman and Hall, 1986.

n=length(xi);

%Add all the terms where j~=k in the expression in [1].
sumVal=0;
derivSumVal=0;
for k=1:(n-1)
    u=(xi(k)-xi((k+1):end))/h;
    u2=u.*u;
    u4=u2.*u2;
    dudh=-u/h;
    
    expVal=exp(-(1/4)*u2);
    Psi=(1-u2+u4/12).*expVal;
    sumVal=sumVal+sum(Psi);
    
    PsiDeriv=-(1/24)*expVal.*u.*(60-20*u2+u4).*dudh;
    derivSumVal=derivSumVal+sum(PsiDeriv);
end

%Due to the symmetry of Psi and PsiDeriv, each of those terms appears
%twice.
sumVal=2*sumVal;
derivSumVal=2*derivSumVal;

%The initial sum value is Phi(0) times n -the values of all instances
%where j=k in the two sums in [1]. Add in those terms.
sumVal=sumVal+n;
%The terms for xi==xj in PsiDeriv are all zero.

c=(3/8)/(n^2*sqrt(pi));

rootn5cPhi=nthroot(c*sumVal,5);

beta=h/rootn5cPhi;
dbetadh=rootn5cPhi*(1-(1/5)*(h*derivSumVal)/sumVal);

%The value of alpha for a Gaussian kernel. alpha is defined before equation
%3.45 in [1].
alpha=(1/(2*sqrt(pi)))^(1/5);
hEst=alpha*beta*n^(-1/5);
dhEstdh=alpha*dbetadh*n^(-1/5);

end

function [hEst,dhEstdh,alphaHat2ConstTerm]=h2S(h,xi,alphaHat2ConstTerm)
%%H2S This function finds the value of h2S from a given value of h as
%     described in Section 5 of [1]. This function does not implement
%     Newton's method; it just returns h2S given h as well as the
%     derivative of h2S with respect to the given h. This function can
%     provide the values to be used in Newton's method to solve Equation
%     12 in [1].
%
%INPUTS: h The value of h going into the left-most term of Equation 12 in
%          [1].
%       xi The scalar sample points.
% alphaHat2ConstTerm A term that depends on xi and not h. Passing this term
%          avoid repreated computations. It is the expression for
%          \hat{\alpha}}_2 in Section 5 of [1] without the h term.
%
%OUTPUTS:      hEst The value of the left-most term of Equation 12 in [1].
%           dhEstdh The derivative of hEst woith respect to h.
%alphaHat2ConstTerm The term alphaHat2ConstTerm computed from xi that can
%                   be passed on a subsequent call to this function to make
%                   it faster.
%
%REFERENCES
%[1] S. J. Sheather and M. C. Jones, "A reliable data-based bandwidth
%    selection method for kernel density estimation," Journal of the Royal
%    Statistical Society. Series B (Methodological), vol. 53, no. 3, pp.
%    683-690, 1991.

n=length(xi);

if(nargin<3||isempty(alphaHat2ConstTerm))
    %The scale factor is the sample interquartile range.
    lambdaHat=interquartileRange(xi);

    a=0.920*lambdaHat*n^(-1/7);
    b=0.912*lambdaHat*n^(-1/9);
    TDb=0;
    SDa=0;
    for k=1:(n-1)
        xDiff=(xi(k)-xi((k+1):end));
        xDiff2=xDiff.*xDiff;
        TDb=TDb+sum(exp(-(1/2)*xDiff2/b^2))/(sqrt(2*pi));
        SDa=SDa+sum(exp(-(1/2)*xDiff2/a^2))/(sqrt(2*pi));
    end
    %Deal with symmetry of the cross terms.
    TDb=2*TDb;
    SDa=2*SDa;
    %Add in the diagonal terms.
    TDb=TDb+n;
    SDa=SDa+n;
    %The sign of TDb has been changed from the paper, since TDb should be
    %positive.
    TDb=TDb/(b^7*n*(n-1));
    SDa=SDa/(a^5*n*(n-1));

    alphaHat2ConstTerm=1.357*nthroot(SDa/TDb,7);
end

alphaHat2=alphaHat2ConstTerm*h^(5/7);

%The derivative of alphaHat2 with respect to h.
dalphaHat2dh=(5/7)*alphaHat2ConstTerm/nthroot(h^2,7);

SDSum=0;
dSDdalphaHat2Sum=0;%Fro the derivative with respect to alpha.
for k=1:(n-1)
    xDiff=(xi(k)-xi((k+1):end));
    xDiff2=xDiff.*xDiff;
    expVal=exp(-(1/2)*xDiff2/alphaHat2^2);
    SDSum=SDSum+sum(expVal)/(sqrt(2*pi));
    
    dSDdalphaHat2Sum=dSDdalphaHat2Sum+sum((expVal.*xDiff2)/(sqrt(2*pi)*alphaHat2^3));
end
SDSum=2*SDSum;
dSDdalphaHat2Sum=2*dSDdalphaHat2Sum;
SDSum=SDSum+n;%Diagonal terms.
%The derivative with respect to alphaHat2 has no diagonal terms in the sum
%part.
SD=SDSum/(alphaHat2^5*n*(n-1));

%The derivative of SD with respect to alphaHat2 (chain rule used).
dSDdalphaHat2=(1/(n*(n-1)))*(-5*SDSum/alphaHat2^6+dSDdalphaHat2Sum/alphaHat2^5);

%Using the chain rule to get the derivative of SD with respect to h.
dSDdh=dSDdalphaHat2*dalphaHat2dh;

beta=1/nthroot(SD,5);
%The chain rule to get the derivative of beta with respect to h.
dbetadh=-(1/5)*dSDdh/nthroot(SD^6,5);

%The value of alpha for a Gaussian kernel. alpha is defined before equation
%3.45 in [1].
alpha=(1/(2*sqrt(pi)))^(1/5);
hEst=alpha*beta*n^(-1/5);

%Derivative of hEst with respect to h.
dhEstdh=alpha*dbetadh*n^(-1/5);

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
