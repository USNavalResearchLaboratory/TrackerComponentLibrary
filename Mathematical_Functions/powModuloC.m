function r=powModuloC(b,n,c)
%%POWMODULOC For integer values b, n and c, evaluate the value mod(b^n,c).
%          This function can work in many instances where b^n would greatly
%          exceed the finite precision accuracy of the processor.
%
%INPUTS: b The positive, scalar, integer valued base of the exponentiation. 
%        n The positive, scalar, integer valued exponent.
%        c The positive, scalar, integer valued modulo value.
%
%OUTPUTS: r The integer value mod(b^n,c).
%
%This function implements the algorithm described in Section 3 of [1],
%where it is noted that the maximum intermediate computed value does not
%exceed c^2.
%
%EXAMPLE 1:
% r=powModuloC(37,13,5)
%r will be 2, which is correct. It is, however, not the value that one gets
%directly with the value mod(37^13,5), because 37^13 has more digits than
%can be represented with a double precision floating point number.
%
%EXAMPLE 2:
% r=powModuloC(6,3,5)
%One will get r=1. In this case, it would match mod(6^3,5), because there
%is no loss of precision in the direct expression.
%
%REFERENCES:
%[1] D. Bailey, P. Borwein, and S. Plouffe, "On the rapid computation of
%    various polylogarithmic constants," Mathematics of Computation, vol.
%    66, no. 218, pp. 903-913, Apr. 1997.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%All inputs must be integers
t=prevPow2Val(n);

r=1;
while(1)
    if(n>=t)
        r=mod(b*r,c);
        n=n-t;
    end
    t=t/2;
    if(t>=1)
       r=mod(r*r,c);
    else
        break;
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
