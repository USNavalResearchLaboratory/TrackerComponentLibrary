function [c,s]=GivensCS(a,b)
%%GIVENSCS Given real or complex scalar values a and b, compute c and s
%          such that M'*[a;b]=[r;0], where M=[c, s;-s',c]' is a rotation
%          matrix, c and s are the cosine and sine of an angle, and r is a
%          scalar value. If a and b are real, then r=norm([a;b]);
%
%INPUTS: a,b Two scalar parameters of a vector [1;b] that are to be rotated
%            to [r;0]; 
%
%OUTPUTS: c,s The cosine and sine parameters of the rotation matrix.
%
%If a and b are real, then Algorithm 5.1.3 of Chapter 5.1.8 of [1] is used.
%If either is complex, then the algorithm of Chapter 5.1.13 is used, which
%calls Algorithm 5.1.13 recursively.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isreal(a) && isreal(b))
    %Algorithm 5.1.3
    if(b==0)
        c=1;
        s=0;
    else
        if(abs(b)>abs(a))
            tau=-a/b;
            s=1/sqrt(1+tau^2);
            c=s*tau;
        else
            tau=-b/a;
            c=1/sqrt(1+tau^2);
            s=c*tau;
        end
    end

else
    %The algorithm of Chapter 5.1.13. 
    [cAlpha,sAlpha]=GivensCS(real(a),imag(a));
    [cBeta,sBeta]=GivensCS(real(b),imag(b));
    [cTheta,sTheta]=GivensCS(norm(a),norm(b));
    
    ejTheta=(cAlpha*cBeta+sAlpha*sBeta)+1j*(cAlpha*sBeta-cBeta*sAlpha);
    c=cTheta;
    s=sTheta*ejTheta;
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
