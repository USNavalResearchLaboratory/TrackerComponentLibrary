function lambda=ChebyshevI2IIMoments(moments)
%%CHEBYSHEVI2IIMOMENTS Let moments hold values of the integrals
%           int_{-1}^1 T_r(x) f(x) dx for some value x and for r=0:n, where
%           T_n(x) is the nth order Chebyshev polynomial of the first kind.
%           This function converts the moments into those from the integral
%           int_{-1}^1 U_r(x) f(x) dx, where U_n(x) is the nth order
%           Chebyshev polynomial of the first kind.
%
%INPUTS: moments A 1X(n+1) or (n+1)X1 vector of moments of f(x) with
%                respect to Chebyshev polynomials of the first kind from
%                order 0 to n.
%
%OUTPUTS: lambda A(n+1)X1 vector of moments of f(x) with respect to
%                Chebyshev polynomials of the second kind.
%
%The conversion is given in remark 1 in Section 3 of [1].
%
%REFERENCES:
%[1] A. Sommariva, "Fast construction of Fejér and Clenshaw-Curtis rule for
%    general weight functions," Computers and Mathematics with
%    Applications, vol. 65, no. 4, pp. 682-693, Feb. 2013.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

n=length(moments);

lambda=zeros(n,1);
lambda(1)=moments(1);

%Odd entries
sumVal=0;
for k=1:2:(n-1)
    sumVal=sumVal+2*moments(k+1);
    lambda(k+1)=sumVal;
end

%Even entries
sumVal=2*lambda(1);
for k=2:2:(n-1)
    sumVal=sumVal+2*moments(k+1);
    lambda(k+1)=sumVal-lambda(1);
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
