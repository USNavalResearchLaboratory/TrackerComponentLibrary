function H=numDiffHess(x,f,fDim,epsilon,algorithm)
%%NUMDIFFHESS Numerically estimate a matrix of second derivatives,
%         including cross terms, of a given function f at the point x.
%         Formulae are derived by repeatedly applying forward differencing,
%         central differencing, or backward differencing. For functions f
%         that return vector values, the matrix H is such that H(:,:,i) is
%         the (symmetric) Hessian matrix for the ith component of the
%         output of f. The algorithm is made to alert (and not crash) in
%         the event that f fails. Specifically, if f returns an empty
%         matrix or it any of the components of f are NaNs, then the
%         numDiff function will terminate early, returning an empty matrix
%         to indicate failure.      
%
%INPUTS: x The xDimX1 vector or scalar point at which the derivative of the
%          (possibly vector) function is desired.
%        f The scalar or vector function that is to be differentiated. The
%          function f must take x as its parameter and its output should be
%          a scalar or a column vector.
%     fDim The dimensionality of the output of f.
%  epsilon A scalar or xDimX1 vector quantity specifying the finite step
%          size used for numerical differentiation. If a scalar value is
%          given, that value is used for differentiating with respect to
%          elements of xDim. If an xDimX1 value is given, then the
%          corresponding element of epsilon is used to differentiate each
%          element of x. If epsilon is omitted or an empty matrix is
%          passed, then epsilon=max(1e-5*x,1e-7); is used.
% algorithm A value indicating the finite difference method used. Possible
%          values are:
%          0 Forward differencing.
%          1 (The default if omitted or an empty matrix is used) central
%            differencing.
%          2 Backward differencing.
%
%OUTPUTS: H An xDimXxDimXfDim set of Hessian matrices. If at any point the
%           function f returned a NaN or an empty matrix, H will be an
%           empty matrix.
%
%Forward differencing to estimate the first derivatives of f with respect
%to x(i) is:
%dfdx=(f(x+ei)-f(x))/epsilon(i)
%where ei is a vector of zeros with an epsilon(i) in the ith entry.
%Backward differencing is 
%dfdx=(f(x)-f(x-ei))/epsilon(i)
%and central differencing is 
%dfdx=(f(x+ei)-f(x-ei))/(2*epsilon(i))
%To get second derivatives, we do finite differencing on dfdx. The finite
%differencing above are just first order estimates. It can be shown that
%the error for the forward and backward differences scales O(epsilon), but
%the error for central differencing scales as O(epsilon^2), so central
%differences are generally preferable.
%
%EXAMPLE 1:
%Here we compare the numeric Hessian of a bivariate function with its
%analytic solution.
% x=[pi;exp(1)];
% eps=1e-3;
% fDim=1;
% H=numDiffHess(x,f,fDim,eps);
% abs(ddf(x)-H)./abs(ddf(x))
%The realtive error will be less than 1e-8, indicating good agreement.
%
%EXAMPLE 2:
%This is an example where f returns a vector.
% f=@(x)[-18+x(1)-2*x(2)+5*x(2)^2-x(2)^3;
%        -34+x(1)*x(2)^3-14*x(2)+x(2)^2+x(2)^3];
% ddf1=@(x)[0,0;
%           0,10-6*x(2)];
% ddf2=@(x)[0,3*x(2)^2;
%           3*x(2)^2,2+6*x(2)+6*x(1)*x(2)];
% x=[12;3];
% HExact=zeros(2,2,2);
% HExact(:,:,1)=ddf1(x);
% HExact(:,:,2)=ddf2(x)
% fDim=2;
% H=numDiffHess(x,f,fDim,epsVal)
%One will see that H and HExact are very close.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algorithm))
    algorithm=1;
end

xDim=size(x,1);

if(nargin<4||isempty(epsilon))
    %If epsilon is not specified, then use some ad-hoc default value
    epsilon=max(1e-5*x,1e-7);
end
if(isscalar(epsilon))
   epsilon=repmat(epsilon,[xDim,1]); 
end

if(algorithm==2)
    epsilon=-epsilon;
    algorithm=0;
end

H=zeros(xDim,xDim,fDim);
switch(algorithm)
    case 0%Order 1,forward difference
        fxDelta=zeros(fDim,xDim);
        e1=zeros(xDim,1);
        for curDim=1:xDim
            e1(curDim)=1;
            fVal=f(x+e1*epsilon(curDim));
            if(isempty(fVal)||any(isnan(fVal)))
                H=[];
                return;
            end
            fxDelta(:,curDim)=fVal;
            e1(curDim)=0;
        end

        e2=zeros(xDim,2);
        fx=f(x);
        if(isempty(fx)||any(isnan(fx)))
            H=[];
            return;
        end
        for curDim1=1:xDim
            e1(curDim1)=1;
            epsH=epsilon(curDim1);
            h=epsH*e1;
            fxh=fxDelta(:,curDim1);

            H(curDim1,curDim1,:)=reshape((fx-2*fxh+f(x+2*h))./(epsH^2),[1,1,fDim]);
            for curDim2=(curDim1+1):xDim
                e2(curDim2)=1;
                epsK=epsilon(curDim2);
                k=epsK*e2;
                fxk=fxDelta(:,curDim2);
                fxkh=f(x+k+h);
                if(isempty(fxkh)||any(isnan(fxkh)))
                    H=[];
                    return;
                end

                H(curDim1,curDim2,:)=reshape((fx-fxh-fxk+fxkh)/(epsH*epsK),[1,1,fDim]);
                H(curDim2,curDim1,:)=H(curDim1,curDim2,:);

                e2(curDim2)=0;
            end
            e1(curDim1)=0;
        end
    case 1%Order 1, central difference.
        e1=zeros(xDim,1);
        e2=zeros(xDim,2);
        for curDim1=1:xDim
            e1(curDim1)=1;
            epsH=epsilon(curDim1);
            h=epsH*e1;
            for curDim2=1:xDim
                e2(curDim2)=1;
                epsK=epsilon(curDim2);
                k=epsK*e2;
                f1=f(x-h-k);
                f2=f(x+h-k);
                f3=f(x-h+k);
                f4=f(x+h+k);

                if(isempty(f1)||any(isnan(f1))||isempty(f2)||any(isnan(f2))||isempty(f3)||any(isnan(f3))||isempty(f4)||any(isnan(f4)))
                    H=[];
                    return;
                end

                H(curDim1,curDim2,:)=reshape((f1-f2-f3+f4)/(4*epsH*epsK),[1,1,fDim]);
                H(curDim2,curDim1,:)=H(curDim1,curDim2,:);
                e2(curDim2)=0;
            end
            e1(curDim1)=0;
        end
    otherwise
        error('Unknown differencing algorithm specified.')
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
