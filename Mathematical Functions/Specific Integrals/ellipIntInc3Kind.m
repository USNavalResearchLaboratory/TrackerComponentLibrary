function intVal=ellipIntInc3Kind(phi,m,n)
%%ELLIPINT3KIND Evaluate Legendre's incomplete elliptical integral of
%               the third kind. This is the integral from 0 to phi of
%               (1-n*sin(theta)^2)^(-1)*(1-m*sin(theta)^2)^(-1/2) dtheta.
%               The integral is sometimes referred to as PI.
%
%INPUTS: phi   A scalar or matrix of the real upper bounds of integration.
%              These should be between -pi/2 and pi/2.
%           m  A value between 0 and 1. The integral diverges at m=1 and
%              phi=pi/2.
%           n  An arbitrary real value. The function is divergent at
%              phi=pi/2 and n=1. For n>2, the result is complex.
%
%OUTPUTS: intVal The value of the incomplete elliptic integral of the
%                third kind with the given parameters.
%
%The expression used to implement the function comes from [1], whereby
%negative values are handled by switching the sign of the result.
%
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numItems=length(phi(:));
intVal=zeros(size(phi));

%The definition of n above is the opposite that used in the paper. 
n=-n;

for curItem=1:numItems
    phiCur=phi(curItem);

    phiSign=sign(phiCur);
    phiCur=abs(phiCur);
    
    %The test deals with precision limitations that would prevent the other
    %integrals from converging.
    if(phiCur<eps)
        intVal(curItem)=0;
    else
        c=csc(phiCur)^2;
        intVal(curItem)=phiSign*symIntFirstKind(c-1,c-m,c)-(n/3)*symIntThirdKind(c-1,c-m,c,c+n);
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
