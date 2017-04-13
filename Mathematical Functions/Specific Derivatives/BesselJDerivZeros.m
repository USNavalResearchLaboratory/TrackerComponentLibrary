function [z,succeeded]=BesselJDerivZeros(m,N,epsVal,NMax)
%%BESSELJDERIVZEROS Determine the first N zeros (in magnitude) of the first
%              derivative of the Bessel function of the first kind. This
%              Bessel function is evaluated using the besselj function.
%              This function does not return 0 as a zero.
%
%INPUTS: m The order of the function. This can be any real number >-1.
%        N The number of zeros of the derivative of the Bessel function to
%          find.
%   epsVal An optional parameter for the approximate desired relative
%          accuracy. The default if this parameter is omitted or an emty
%          matrix is passed is 1e-14.
%     NMax The algorithm finds Bessel function zeros by finding the
%          eigenvalues of a square tridiagonal matrix. This is the maximum
%          matrix size to consider. If this parameter is omitted or an
%          empty matrix is passed, the NMax=1000 is used. Note that if NMax
%          is below N, then just N is used. 
%
%OUTPUTS: z A 1XN vector of the first N zeros to the first derivative of
%           the Bessel function of the first kind. These are always real
%           and zero is not included.
% succeeded This is true if the desired accuracy epsVal was achieved. If
%           false, one can try increasing NMax.
%
%This function implements the algorithm described in [1]. The solution to
%the problem is expressed as the eigenvalues of an infinite dimensional
%matrix. Here, we just need to determine what size matrix to truncate to to
%get the desired precision. For very large problems, finite precision
%issues can arise.
%
%EXAMPLE:
%The first 30 roots of the derivative of an order 7.5 Bessel function of
%the first kind:
% m=7.5;
% N=30;
% [z,succeeded]=BesselJDerivZeros(m,N);
%One will see that BesselJDeriv is less than 1e-13.
%
%REFERENCES:
%[1] Y. Ikebe, Y. Kikuchi, and I. Fujishiro, "Computing zeros and orders of
%   Bessel functions," Journal of Computational and Applied Mathematics,
%   vol. 38, no. 1-3, pp. 169-184, 23 Dec. 1991.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(epsVal))
    epsVal=1e-14;
end

if(nargin<4||isempty(NMax))
   NMax=1000; 
end

B=zeros(NMax,NMax);

g1=(4+3*m)/(m*(m+1)*(m+2));
B(1,1)=g1;

%Fill in the initial minGuess values. An initial guess is from the
%asymptotic approximation mentioned in [1] for the function, not its
%derivative. However, one would expect zero derivatives to occur two to
%three times as zeros, so this is probably a decent upper bound.
minGuess=max(min(ceil((pi/2)*N),NMax),N);

%Fill in the matrix.
for k=2:(minGuess-1)
    %The following equations are taken from Section 4 of [1].
    gk=2/((m+2*k-2)*(m+2*k));
    hk=1/((m+2*k-2)*sqrt((m+2*k-3)*(m+2*k-1)));
    B(k,k)=gk;
    B(k-1,k)=hk;
    B(k,k-1)=hk;
end

%Now, deternine the error. If it is small enough stop. Otherwise, keep
%expanding the matrix.
for k=minGuess:NMax
    %Increase the order by one.
    %The following equations are taken from Section 4 of [1].
    gk=2/((m+2*k-2)*(m+2*k));
    hk=1/((m+2*k-2)*sqrt((m+2*k-3)*(m+2*k-1)));
    B(k,k)=gk;
    B(k-1,k)=hk;
    B(k,k-1)=hk;
    
    %Obtain approximate roots.
    lambda=sort(eig(B(1:k,1:k)),'descend');
    %This is form Theorem 4.1
    z=2./sqrt(lambda(1:N));%Approximate roots
    
    %The upper bound is from Equation 4.5 in [1]. Here, we use the estimated
    %zeros z in place of the true ones for computing the bound.
    n=k;
    JTemp01=abs(besselj(m+2*n,z).*besselj(m+2*n+2,z));
    denom1=2*(1-(m./z).^2).*besselj(m,z).^2*(m+2*n);
    errBound=max(JTemp01./denom1+(JTemp01./(denom1*(m+2*n)*(m+2*n-1))).*(z.^3/(8*pi)));
    if(errBound<epsVal||k==NMax)
        if(errBound>=epsVal)
            succeeded=false;
        else
            succeeded=true;
        end
        return;
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
