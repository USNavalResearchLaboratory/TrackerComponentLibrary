function intVal=ellipIntInc3Kind(phi,m,n)
%%ELLIPINT3KIND Evaluate Legendre's incomplete elliptical integral of
%               the third kind. This is the integral from 0 to phi of
%               (1-n*sin(theta)^2)^(-1)*(1-m*sin(theta)^2)^(-1/2) dtheta.
%               The integral is sometimes referred to as PI.
%
%INPUTS: phi A scalar or matrix of the real upper bounds of integration.
%          m A value between 0 and 1. The integral diverges at m=1 and
%            phi=pi/2.
%          n An arbitrary real value. The function is divergent at
%            phi=pi/2 and n=1. For n>2, the result is complex.
%
%OUTPUTS: intVal The value of the incomplete elliptic integral of the
%                third kind with the given parameters.
%
%The expression used to implement the function comes from [1], where
%Equation 61 is the basic identity. However, the first symmetric integral
%has been replaced using the identity in equation 59. A similar
%transformation was applied to the second symmetric integral so that all
%cosecant terms would be replaced by sines or cosines so as to avoid non-
%finite numbers.
%
%However, the resulting formula is only valid for values of phi from 0 to
%pi/2. However, the argument of the integral is just mirrored from pi/2 to
%pi. From then on, the function repeats again from 0. Thus, first, we
%determine how many intervals of pi/2 are present. We can easily evaluate
%the integral from 0 to pi/2, so we subtract out those values, leaving a
%residual to compute. The residual is  phi=phi-numMult*(pi/2), where
%numNult is the integer number of pi/2 values in the original phi. The
%residual phi is the fractional part of a pi/2-length section left. If
%numMult is even, then we can directly use the formula. However, if numMult
%is odd, it means that the residual part is on a mirror section. Thus, the
%residual integral is from pi/2-phi to pi/2, which can be evaluated as the
%integral from to to pi/2 minus the integral from 0 to pi/2-phi.
%
%However, the above method assumes that phi is positive. If phi is
%negative, then we evaluate the integral using abs(phi) and then flip the
%sign of it. This can be done due to the symmetry of the argument of the
%integral about the origin.
%
%REFERENCES:
%[1] B. C. Carlson, "Numerical computation of real or complex elliptic
%    integrals," Numerical Algorithms, vol. 10, no. 1, pp. 13-26, 1995.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numItems=length(phi(:));
intVal=zeros(size(phi));

for curItem=1:numItems
    phiCur=phi(curItem);

    phiSign=sign(phiCur);
    phiCur=abs(phiCur);
    
    if(n==1)
        %A limit value obtained from Mathematica...
        intVal(curItem)=ellipIntInc1Kind(phi,m)+(ellipIntInc2Kind(phi,m)-tan(phi)*sqrt(1-m*sin(phi)^2))/(m-1);
    else
        numMult=fix(phiCur/(pi/2));
        if(numMult>0)
            phiCur=phiCur-numMult*(pi/2);

            completeVal=(symIntFirstKind(0,1-m,1)+(n/3)*symIntThirdKind(0,1-m,1,1-n));

            intVal(curItem)=numMult*completeVal;
            if(mod(numMult,2))
               phiCur=pi/2-phiCur; 
            end
        end

        c2=cos(phiCur)^2;
        s=sin(phiCur);
        s2=s*s;

        %Equation 61, modified with Equation 59 and a similar modification for
        %the second symmetric integral.
        fracVal=(s*symIntFirstKind(c2,1-m*s2,1)+(n/3)*s*s2*symIntThirdKind(c2,1-m*s2,1,1-n*s2));

        if(mod(numMult,2))
            intVal(curItem)=phiSign*(intVal(curItem)+completeVal-fracVal);
        else
            intVal(curItem)=phiSign*(intVal(curItem)+fracVal);
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
