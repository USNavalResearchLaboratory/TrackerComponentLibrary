classdef SinPowD
%%SINPOWD Functions to hand the scalar sine power distribution. This is a
%         distribution where the probability density function (PDF) is
%         given by sin(x)^n for 0<=x<=pi. This distribution arises when
%         generating uniformly distributed samples within a hyperellipsoid
%         as in [1].
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand
%
%REFERENCES:
%[1] H. Sun and M. Farooq, "Note on the generation of random points
%    uniformly distributed in hyper-ellipsoids," in Proceedings of the
%    Fifth International Conference on Information Fusion, Annapolis, MD,
%    8-11 Jul. 2002, pp. 489-496.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(n)
%%MEAN  Obtain the mean of the scalar sine power probability distribution.
%
%INPUTS: n The real integer power of the sine term in the PDF (n>=0).
%
%OUTPUTS: val The mean of the PDF.
%
%This function simply implments an explicit solution, which was found by
%integration.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    normConst=intSinPow([0;pi],n);
    val=(1/normConst)*(pi^(3/2)*gamma((1+n)/2))/(2*gamma(1+n/2));
end

function val=var(n)
%%VAR  Obtain the variance of the scalar sine power probability
%      distribution.
%
%INPUTS: n The real integer power of the sine term in the PDF (n>=0).
%
%OUTPUTS: val The variance of the PDF.
%
%This function simply implments an explicit solution, which was found by
%integration.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    normConst=intSinPow([0;pi],n);
    E2Val=(1/normConst)*intPowSinPow([0;pi],2,n);%Second noncentral moment.
    val=E2Val-SinPowD.mean(n)^2;
end

function val=PDF(x,n)
%%PDF Evaluate the scalar sine power probability distribution function
%     (PDF) at one or more desired points.
%
%INPUTS: x The vector or matrix of points at which the PDF should be
%          evaluated.
%        n The real integer power of the sine term in the PDF (n>=0).
%
%OUTPUTS: val The value(s) of sine power PDF evaluated at x.
%
%This function simply goes from the definition of the PDF, using the
%intSinPow function to obtain the normalizing constant.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

normConst=intSinPow([0;pi],n);
val=sin(x).^n/normConst;

val(x<0||x>pi)=0;

end

function val=CDF(x,n)
%%CDF Evaluate the cumulative distribution function (CDF) of the
%     scalar sine power distribution at desired points.
%
%INPUTS: x The vector or matrix of points at which the CDF should be
%          evaluated.
%        n The real integer power of the sine term in the PDF (n>=0).
%
%OUTPUTS: val The value(s) of sine power CDF evaluated at x.
%
%The function intSinPow provides the integral of the non-normalized PDF.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

normConst=intSinPow([0;pi],n);
val=intSinPow([0;x(:)'],n)/normConst;
 
val=reshape(val,size(x));
val(x<0)=0;
val(x>=pi)=1;

end

function val=invCDF(prob,n,options)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the scalar sine power distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%           n The real integer power of the sine term in the PDF (n>=0).
%     options The function utilizes the fminbnd function. This is the
%             optional set of option parameters for the fminbnd function.
%             If this parameter is omitted or an empty matrix is passed,
%             the defaults for the fminbnd function are used except TolX is
%             set to 1e-12.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The sine distribution is complicated and not easily inverted, so the
%fminbnd function is used to invert the CDF.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(options))
    options=optimset('fminbnd');
    options=optimset(options,'TolX',1e-12);
elseif(isempty(optimget(options,'TolX')))
    options=optimset(options,'TolX',1e-12);
end

normConst=intSinPow([0;pi],n);
prob=prob(:);
N=length(prob);

val=zeros(size(prob));

for curEl=1:N
    pCur=prob(curEl);
    if(pCur<=0)
        val(curEl)=0;
        continue;
    elseif(pCur>=1)
        val(curEl)=pi;
        continue;
    end
    
    f=@(u)costFun(u,pCur);
    val(curEl)=fminbnd(f,0,pi,options);
end

function x=costFun(u,p)
    x=(intSinPow([0;u(:)'],n)/normConst-p)^2;
end
end

function vals=rand(N,n,algorithm,options)
%%RAND Generate sine power distributed random variables with the given
%      parameter.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element  row vector, then rand
%          returns an MXN1 matrix of random variables.
%        n The real integer power of the sine term in the PDF (n>=0).
% algorithm A parameter specifying the algorithm to use to generate the
%          random variables. Possible values are
%          0 (The default if this parameter is omitted or an empty matrix
%            is passed) Use the rejection method of Chapter 5.2 of [1].
%          1 Use the inverse transformation method of Chapter 5.1 of [1].
%  options If algorithm=1, then this function utilizes the fminbnd
%          function. This input is the optional set of option parameters
%          for the fminbnd function. If this parameter is omitted or an
%          empty matrix is passed, the defaults for the fminbnd function
%          are used except TolX is set to 1e-12. This parameter is not used
%          if algorithm=0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated sine power random variables.
%
%EXAMPLE:
%We can verify that the sample mean and variance are close to the true mean
%and variance.
% N=[1e5,1];
% n=4;
% vals=SinPowD.rand(N,n);
% mean(vals)
% SinPowD.mean(n)
% var(vals)
% SinPowD.var(n)
%One will see that the sample mean and variance are close to the true mean
%and variance.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(algorithm))
       algorithm=0; 
    end

    if(nargin<4||isempty(options))
       options=[]; 
    end

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end
    
    switch(algorithm)
        case 0%The rejection method of Chapter 5.2 of [1].
            %We choose the sample density to be the uniform distribution.
            %If g(x) is the uniform distribution and f(x) is the desired
            %distribution, we need to find c such that f(x)/g(x)=c if we
            %want to minimize the number of rejections. The solution is pi
            %divided by the normalizing constant for the sine power
            %distribution. Since g(x)=pi, the pi cancels out of
            %f(x)/(c*g(x) as does the normalizing constant.
            
            vals=zeros(dims);
            numVals=prod(dims);
            
            for curVal=1:numVals
                while(1)
                    Y=pi*rand(1);
                    U=rand();
                    if(U<=sin(Y)^n)
                        vals(curVal)=Y;
                        break;
                    end
                end
            end
        case 1%The inverse transform method of Chapter 5.1 of [1].    
            U=rand(dims);
            vals=SinPowD.invCDF(U,n,options);
        otherwise
            error('Unknown algorithm specified');
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
