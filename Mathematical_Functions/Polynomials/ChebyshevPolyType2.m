function U=ChebyshevPolyType2(tau,n,tauStart,tauEnd)
%%CHEBYSHEVPOLYTYPE2 Evaluate Chebyshev polynomials of the second kind from
%               order 0 to n order at a given point tau.
%
%INPUTS: tau An NX1 or 1XN vector of values from tauStart to tauEnd, or
%            from -1 to 1 if tauStart and tauEnd are omitted, where one
%            wishes to evaluate the Chebyshev polynomials of the second
%            kind.
%          n The non-negative integer maximum order of the Chebyshev
%            polynomials evaluated.
% tauStart,tauEnd The possible range of the inputs. Given this range, the
%            actual input to the polynomials is scaled to
%            tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
%            If omitted, a range of -1 to 1 is assumed --the normal range
%            for Chebyshev polynomials.
%
%OUTPUTS: T An (n+1)XN matrix of the Chebyshev polynomials of the second
%           kind from order 0 to n evaluated at each of the values of tau.
%           U(i,j) is the (i-1)th order Chebyshev polynomial evaluated at
%           tau(j).
%
%In remark 1 in Section 3 of [1], a formula to convert from Chebyshev
%polynomials of type I to type II is given. That is implemented here and
%the ChebyshevPoly function is called to get Chebyshev polynomials of type
%I.
%
%REFERENCES:
%[1] A. Sommariva, "Fast construction of Fejér and Clenshaw-Curtis rule for
%    general weight functions," Computers and Mathematics with
%    Applications, vol. 65, no. 4, pp. 682-693, Feb. 2013.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(tau);

if(nargin<3||isempty(tauStart))
    tauStart=-1;
end
   
if(nargin<4||isempty(tauEnd))
    tauEnd=1;
end
    
T=ChebyshevPoly(tau,n,tauStart,tauEnd);

U=zeros(n+1,N);

U(1,:)=T(1,:);

%Odd entries
sumVal=zeros(1,N);
for k=1:2:n
    sumVal=sumVal+2*T(k+1,:);
    U(k+1,:)=sumVal;
end

%Even entries
sumVal=2*U(1,:);
for k=2:2:n
    sumVal=sumVal+2*T(k+1,:);
    U(k+1,:)=sumVal-1;
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
