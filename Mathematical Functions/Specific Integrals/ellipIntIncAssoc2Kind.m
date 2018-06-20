function intVal=ellipIntIncAssoc2Kind(phi,m)
%%ELLIPINTINCASSOC2KIND Evaluate the associated incomplete elliptic
%                       integral of the second kind. This is the integral
%                       from 0 to phi of
%                       sin(theta)^2*(1-m*sin(theta)^2)^(-1/2) dtheta.
%                       This integral is sometimes referred to as D.
%
%INPUTS: phi  The real, positive upper bound of integration. This should
%              generally be between -pi/2 and pi/2.
%           m  A value between 0 and 1. The integral diverges when m=1 and
%              phi=pi/2.
%
%OUTPUT: intVal The value of the associated incomplete elliptic integral of
%               the second kind with the given parameters.
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

%The test deals with precision limitations that would prevent the other
%integral from converging.
if(abs(phi)<eps)
    intVal=0;
else
    c=csc(phi)^2;
    intVal=(1/3)*symIntThirdKindDegen(c-1,c-m,c);
    
    if(sign(phi)<0)
        intVal=-intVal;
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
