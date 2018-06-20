function intVal=ellipIntInc1Kind(phi,m)
%%ELLIPINTINC1KIND Evaluate Legendre's incomplete elliptical integral of
%                  the first kind. This is the integral from 0 to phi of
%                  (1-m*sin(theta)^2)^(-1/2) dtheta. Alternatively, this
%                  can be written as the integral from 0 to x of
%                  ((1-t^2)*(1-m*t^2))^(-1/2) dt where x=sin(phi). This
%                  integral is sometimes referred to as F.
%
%INPUTS: phi A scalar or matrix of the real upper bounds of integration.
%            These can be positive, negative or zero.
%          m A value between 0 and 1. The integral diverges when m=1 and
%            phi=pi/2.
%
%OUTPUTS: intVal The values of the incomplete elliptic integral of the
%                first kind evaluated at each entry in phi.
%
%The expression used to implement the function comes from [1], whereby
%negative values are handled by switching the sign of the result. The above
%algorithm also only considers values with magnitudes up to pi/2. To get
%the correct result, however, for larger values, the value up to pi/2
%(times the number of multiples of the value up to pi/2) is added to a
%fractional part with magnitude less than pi/2. Note however, that due to
%the odd symmetry of sin(x)^2, if the integer part is odd, then the
%computation subtracts the fractional part from pi/2 and the value added is
%the value of the complete integral MINUS the modified fractional part.
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
    
    %Separate out the multiples of pi/2 from the rest of the values.
    intPart=fix(phiCur/(pi/2));
    if(intPart>0)
        completeVal=symIntFirstKind(0,1-m,1);
        phiCur=phiCur-intPart*(pi/2);
    else
        completeVal=0;
    end
    
    %The test deals with precision limitations that would prevent the other
    %integral from converging.
    if(phiCur<eps)
        intVal(curItem)=phiSign*intPart*completeVal;
    else
        %Deal with the odd symmetry of sin(theta)^2
        if(mod(intPart,2)==0)
            c=csc(phiCur)^2;
            delta=symIntFirstKind(c-1,c-m,c);
            
            intVal(curItem)=phiSign*(intPart*completeVal+delta);
        else
            phiCur=pi/2-phiCur;
            
            c=csc(phiCur)^2;
            delta=symIntFirstKind(c-1,c-m,c);
            
            intVal(curItem)=phiSign*(intPart*completeVal+(completeVal-delta));
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
