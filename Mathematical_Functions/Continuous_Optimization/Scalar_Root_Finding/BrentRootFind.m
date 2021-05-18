function [zeroVal,exitFlag]=BrentRootFind(f,xSpan,AbsTol,maxIter)
%%BRENTROOTFIND Given two points that are of different signs
%          surrounding the root of a scalar funding, find the root of the
%          function using Brent's method. Brent's method is a combination
%          of bisection, linear interpolation and inverse quadratic
%          interpolation.
%
%INPUTS: f The handle to the function whose root is desired. The function
%          is called as f(x), where x is a scalar value.
%    xSpan A 2X1 or 1X2 vector such that
%          sign(f(xSpan(1)))~=sign(f(xSpan(2)))
%   AbsTol The absolute tolerance on the value of f for converence, as
%          defined in [1]. The default if omitted or an empty matrix is
%          passed is eps().
%  maxIter The maximum number of ierations to perform. The default if
%          omitted or an empty matrix is passed is 100.
%
%OUTPUTS: zeroVal This is the value of x at the zero point.
%        exitCode A value indicating how the function terminated. Possible
%                values are:
%                0 Either AbsTol was satisfied or the bounds on the value
%                  zeroVal went below finite precision constraints.
%                1 The maximum number of iterations occurred.
%
%This function implements the first algorithm described by Brent in [1].
%That is the zero algorithm.
%
%This function usually outperforms bisectionRootFind. However, if the
%function being zeroed is very flat over much of the bounded region, then
%this function might get stuck and return a bad value.
%
%EXAMPLE:
% f=@(x)(x.^2-4);
% xSpan=[0;5];
% [zeroVal,exitFlag]=BrentRootFind(f,xSpan)
%In this example convergence is to the exact root, zeroVal=2.
%
%REFERENCES:
%[1] R. P. Brent, "An algorithm with guaranteed convergence for finding a
%    zero of a function," The Computer Journal, vol. 14, no. 4, pp. 422-
%    425, 1971.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(maxIter))
   maxIter=100; 
end

if(nargin<3||isempty(AbsTol))
   AbsTol=eps(); 
end

xa=xSpan(1);
xb=xSpan(2);

fa=f(xa);
fb=f(xb);

if((fa>0)==(fb>0))
   error('There is no sign change in the initial span.');
end

xc=xa;
fc=fa;
diffBa=xb-xa;
d=diffBa;

exitFlag=1;
for curIter=1:maxIter
    if((fb>0)==(fc>0))
        xc=xa;
        fc=fa;
        diffBa=xb-xa;
        d=diffBa;
    end
    %At this point, we have fb*fc<=0
    
    if(abs(fc)<abs(fb))
        xa=xb;
        xb=xc;
        xc=xa;
        fa=fb;
        fb=fc;
        fc=fa;
    end
    %At this point, we have fb*fc<0 and abs(fb)<=abs(fc),
    %which meets the assumption in Section 3.
    
    tol=2*eps()*abs(xb)+AbsTol;
    xMid=(xc-xb)/2;

    %The extra two conditions deal with x going below finite precision
    %constraints.
    if(~(abs(xMid)>tol&&fb~=0)||xMid==xc||xMid==xb)
        exitFlag=0;
        break;
    end
    
    if(abs(diffBa)<tol||abs(fa)<=abs(fb))
        %Perform a bisection.
        d=xMid;
        diffBa=xMid;
    else
        s=fb/fa;
        if(xa==xc)%Linear interpolation
            p=2*xMid*s;
            q=1-s;
        else%Inverse quadratic interpolation (Section 4)
            %This is done if a, b, and c are distinct.
            q=fa/fc;
            r=fb/fc;
        
            p=s*(s*xMid*q*(q-r)-(xb-xa)*(r-1));
            q=(q-1)*(r-1)*(s-1);
        end
        
        if(p>0)
            q=-q;
        else
            p=-p;
        end
        
        s=diffBa;
        diffBa=d;
        
        if((2*p<3*xMid*q-abs(tol*q))&&p<abs(s*q/2))
            d=p/q;
        else
            d=diffBa;
            diffBa=xMid;
        end
        xa=xb;
        fa=fb;
        if(abs(d)>tol)
            xb=xb+d;
        elseif(xMid>0)
            xb=xb+tol;
        else
            xb=xb-tol;
        end
        fb=f(xb);
    end
end

zeroVal=xb;
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
