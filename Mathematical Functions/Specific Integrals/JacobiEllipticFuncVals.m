function [snVals,cnVals,dnVals]=JacobiEllipticFuncVals(xVals,m)
%%JACOBIELLIPTICFUNCVALS This evaluate the Jacobi elliptic functions sn, cn
%       and dn for real inputs x and m. The functions sn solves the
%       equation
%       x=integral_0^sn 1/sqrt((1-t^2)*(1-m*t^2) dt
%       The functions cn and dn stem from the relations sn^2+cn^2=1 and
%       m*sn^2+dn^2=1. However, for values of x being complex or with
%       abs(x)>1, a definition given in terms of differential equations as
%       in [2] is more appropriate.
%
%INPUTS: x A real or complex scalar value or matrix of values at which the
%          Jacobi elliptic function is to be evaluated. These values must
%          be such that abs(x(k))<=1 for all k.
%        m The real, finite parameter of the integral.
%
%OUTPUTS: sn The upper bound of integration in the above integral and/or
%            the solution to a particular set of differential equations
%            given in [2].
%         cn A value such that sn^2+cn^2=1.
%         dn A value such that m*sn^2+dn^2=1.
%
%This function implements the algorithm for sncndn(x,m) that is given as
%Algorithm 5 of [1]. Note that the definition of m is different from the
%definition of mc in [1] and must be transformed.
%
%REFERENCES:
%[1] R. Burlisch, "Numerical calculation of elliptic integrals and elliptic
%    functions," Numerische Mathematik, vol. 7, no. 1, pp. 78-90, 1965.
%[2] Weisstein, Eric W. "Jacobi Elliptic Functions." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/JacobiEllipticFunctions.html
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=numel(xVals);

%Change to the definition of mc used in [1].
kc=1-m;

%Generally, it appears to converge in five or fewer iterations.
maxDepth=16;

%Tolerance for convergence.
tolVal=sqrt(eps()*0.01);

if(kc~=0)
    mcIsNeg=(kc<0);
    if(mcIsNeg)
       d=1-kc;
       kc=-kc/d;
       d=sqrt(d);
       xVals=d*xVals;
    end
    
    sizeX=size(xVals);
    snVals=zeros(sizeX);
    cnVals=zeros(sizeX);
    dnVals=zeros(sizeX);
    for curX=1:N
        mc=kc;
        x=xVals(curX);
        
        a=1;
        dn=1;
        m=zeros(maxDepth,1);
        n=zeros(maxDepth,1);
        for l=1:maxDepth
            m(l)=a;
            mc=sqrt(mc);
            n(l)=mc;
            c=(a+mc)/2;

            if(~(abs(a-mc)>tolVal*a))
                break;
            end
            mc=mc*a;
            a=c;
        end

        x=x*c;
        sn=sin(x);
        cn=cos(x);
        if(sn~=0)
            a=cn/sn;
            c=c*a;
            for i=l:-1:1
              b=m(i);
              a=a*c;
              c=c*dn;
              dn=(n(i)+a)/(b+a);
              a=c/b;
            end
            a=1/sqrt(c^2+1);
            if(sn<0)
                sn=-a;
            else
                sn=a;
            end
            cn=c*sn;

            if(mcIsNeg)
                temp=cn;
                cn=dn;
                dn=temp;
                sn=sn/d;
            end
        end
        
        snVals(curX)=sn;
        cnVals(curX)=cn;
        dnVals(curX)=dn;
    end
else%The mc=0 case.
    snVals=tanh(xVals);
    cnVals=1./cosh(xVals);
    dnVals=cnVals;
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
