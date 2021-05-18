function [coeffs,exitCode]=multilinMRegress(x,y,hasIntercept,scalFactAlg,MFunction,RelTol,AbsTol,maxIter)
%%MULTILINMREGRESS Perform multilinear regression with outliers using an M
%           estimator with a chosen weight function. The algorithm finds
%           coeffs such that y=coeffs'*x for an input x and y. If
%           hasIntercept, then the equation is y=coeffs'*[x;1] as the last
%           element of coeffs will be an additive constant.
%
%INPUTS: x A pXn set of n p-dimensional samples of the p-dimensional
%          independent parameter. It is required that n>p.
%        y A 1Xn or nX1 vector of the dependant parameter.
% hasIntercept A boolean value indicating whether the fit has a nonzero y
%          intercept. The default if omitted or an empty matrix is passed
%          is true.
% scalFactAlg This input specifies how when given a set of error residuals,
%          a robust estiamte of the standard deviation is obtained.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            scaled median absolute deviation via medAbsDev(e,[],2). 
%          1 Use the scaled Sn estimator SnEst(e,2);
%          2 Use the scaled Qn estimator QnEst(e,2);
%          A function handle taking a vector and returning the scale factor
%          estimate can also be passed.
% MFunction The specifies the weights used. Note that if an M function
%          is used where so many of the w values are zero that there are no
%          longer enough nonzero values of w for the coefficient to be
%          observable, then as opposed to using u as in Chapter 25.2 of
%          [1], the error residuals are centered as e=e-median(e) and the
%          MFunction value is recomputed. Possible values are:
%          A numeric value: This corresponds to the input funcType in the 
%          MEstFunction function with the default other parameters; only
%          the weight output is used. A function handle: The function
%          handle is called. It should be done in the manner of the
%          functions in MEstFunction, but only one output (the weight) is
%          expected. The default if omitted or an empty matrix is passed is
%          the default of the MFunction function.
% RelTol, AbsTol Parameters determining when the algorithm has converged.
%          If diff is the absolute value fo the difference difference
%          between the coefficient estimate from the current and the
%          previous iteration, then convergence is declared if
%          all(diff<=AbsTol)||all(diff<=coeffs*RelTol). The defaults if
%          omitted or empty matrices are passed are respectively 1e-10 and
%          1e-12.
%  maxIter The maximum number of iterations to allow without convergence.
%          The default if omitted or an empty matrix is passed is 10*n.
%
%OUTPUTS: coeffs The regression coefficients as described above.
%       exitCode A value indicating how the algorithm termianted. Possible
%                values are:
%                0 The algorithm converged.
%                1 Termiantion was due to reching maxIter iterations.
%
%This implements the algorithm described in Chapter 25.2 of [1], with the
%modifications described above for when too many values in the weights from
%the MFunction are zero.
%
%EXAMPLE:
%Here, we show that given noise-free samples of a line polluted 10% by
%samples from another line, one is able to find the dominant line exactly
%even though standard least squares does poorly.
% numTrue=1000;
% betaTrue=[14;-8];
% yIntTrue=42;
% xTrue=rand(2,numTrue);
% yTrue=sum(bsxfun(@times,betaTrue,xTrue),1)+yIntTrue;
% numBad=100;
% betaBad=[-80;32];
% yIntBad=10;
% xBad=rand(2,numBad);
% yBad=sum(bsxfun(@times,betaBad,xBad),1)+yIntBad;
% x=[xTrue,xBad];
% y=[yTrue,yBad];
% coeffsStd=multilinRegress(x,y)
% coeffs=multilinMRegress(x,y)
%One will see that coeffs matches [betaTrue;yIntTrue], but coeffsStd are
%off by quite a bit.
%
%REFERENCES:
%[1] N. R. Draper and H. Smith, Applied Regression Analysis, 3rd ed. New
%    York: John Wiley and Sons, Inc., 1998.
%
%September 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(maxIter))
    n=size(x,2);
    maxIter=10*n;
end

if(nargin<7||isempty(AbsTol))
    AbsTol=1e-12; 
end

if(nargin<6||isempty(RelTol))
    RelTol=1e-10; 
end

if(nargin<5||isempty(MFunction))
    MFunction=[];
end

if(nargin<4||isempty(scalFactAlg))
    scalFactAlg=0;
end

if(nargin<3||isempty(hasIntercept))
   hasIntercept=true; 
end

if(hasIntercept)
    numPoints=size(x,2);
    x=[x;ones(1,numPoints)];
end

numUnknowns=size(x,1);

%Estimate the regression coefficients with unweighted ordinary least
%squares.
coeffs=multilinRegress(x,y,[],false);

if(~isa(scalFactAlg,'function_handle'))
    switch(scalFactAlg)
        case 0
            scalFactAlg=@(e)medAbsDev(e,[],2);
        case 1
            scalFactAlg=@(e)SnEst(e,2);
        case 2
            scalFactAlg=@(e)QnEst(e,2);
        otherwise
            error('Unknown scale factor algorithm chosen.')
    end
end

MFunIsHandle=isa(MFunction,'function_handle');

exitCode=1;
for curIter=1:maxIter
    %error residuals
    e=y(:)-x'*coeffs;

    s=scalFactAlg(e);

    %e=e-median(e);
    if(MFunIsHandle)
        w=MFunction(e/s);
    else
        [~,~,w]=MEstFunction(e/s,MFunction);
    end
    
    %Definition  after Equation 25.2.7 in [1] for any possible zero-error
    %terms.
    sel=(e==0);
    w(sel)=1;
    
    %All zero weights can occur in some instances
    if(sum(w~=0)<numUnknowns)
        e=e-median(e);
        if(MFunIsHandle)
            w=MFunction(e/s);
        else
            [~,~,w]=MEstFunction(e/s-median(e/s),MFunction);
        end
    end

    coeffsNew=multilinRegress(x,y,w,false);

    diff=abs(coeffsNew-coeffs);
    
    coeffs=coeffsNew;
    
    if(all(diff<=AbsTol)||all(diff<=coeffs*RelTol))
        exitCode=0;
        break;%The algorithm converged.
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
