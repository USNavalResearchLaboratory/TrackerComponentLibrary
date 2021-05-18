function [z,succeeded]=BesselJZeros(m,N,epsVal,NMax)
%%BESSELJZEROS Determine the first N zeros (in magnitude) of the Bessel
%              function of the first kind, which is given by the besselj
%              function.
%
%INPUTS: m The order of the function. This can be any real number >-1 or any
%          negative real number that is not an integer. Values >-1 will
%          produce all positive, real solutions. Values <-1 will produce a
%          number of complex solutions of low magnitude and then positive
%          real solutions.
%        N The number of zeros to find.
%   epsVal An optional parameter for the approximate desired relative
%          accuracy. The default if this parameter is omitted or an emty
%          matrix is passed is 1e-14. There is no need for an absolute
%          accuracy parameter as Bessel functions do not have zeros at
%          zero. 
%     NMax The algorithm finds Bessel function zeros by finding the
%          eigenvalues of a square tridiagonal matrix. This is the maximum
%          matrix size to consider. If this parameter is omitted or an
%          empty matrix is passed, the NMax=1000 is used. Note that if NMax
%          is below N, then just N is used. 
%
%OUTPUTS: z A 1XN vector of the first N zeros to a Bessel function of the
%           first kind.
% succeeded This is true if the desired accuracy epsVal was achieved. If
%           false, one can try increasing NMax.
%
%This function implements the algorithm described in [1], [2] and [3]. The
%solution to the problem is expressed as the eigenvalues of an infinite
%dimensional matrix. Here, we just need to determine what size matrix to
%truncate to to get the desired precision. For very large problems, finite
%precision issues can arise.
%
%EXAMPLE:
%The first 30 roots of an order 7.5 Bessel function of the first kind:
% m=7.5;
% N=30;
% [z,succeeded]=BesselJZeros(m,N)
% One will see that besselj(m,z) is less than 1e-13.
%
%REFERENCES:
%[1] Y. Ikebe, Y. Kikuchi, and I. Fujishiro, "Computing zeros and orders of
%   Bessel functions," Journal of Computational and Applied Mathematics,
%   vol. 38, no. 1-3, pp. 169-184, 23 Dec. 1991.
%[2] Y. Miyazaki, N. Asai, Y. Kikuchi, D. Cai, and Y. Ikebe, "Computation
%    of multiple eigenvalues of infinite tridiagonal matrices," Mathematics
%    of Computation, vol. 73, no. 246, pp. 719-730, Apr. 2004.
%[3] Y. Ikebe and Y. Kikuchi, "The eigenvalue problem for infinite compact
%    complex symmetric matrices with application to the numerical
%    computation of comp[lex zeros of J0(z)-i*J1(z) and of Bessel functions
%    Jm(z) of any real order m," Linear Algebra and its Applications, vol.
%    194, pp. 35-70, 15 Nov. 1993.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(epsVal))
    epsVal=1e-14;
end

if(nargin<4||isempty(NMax))
   NMax=1000; 
end

A=zeros(NMax,NMax);

alpha1=m+2*1;
d1=2./((alpha1-1).*(alpha1+1));
A(1,1)=d1;

%Fill in the initial minGuess values. This initial guess is from the
%asymptotic approximation mentioned in [1].
minGuess=max(min(ceil((pi/2)*N),NMax),N);
%Fill in the matrix.
for k=2:minGuess
    %The following equations are taken from Theorem 1.7 of [3].
    alphak=m+2*k;
    dk=2/((alphak-1)*(alphak+1));
    fk=1/((alphak-1)*sqrt(alphak-2)*sqrt(alphak));
    A(k,k)=dk;
    A(k-1,k)=fk;
    A(k,k-1)=fk;
end

%Now, deternine the error. If it is small enough stop. Otherwise, keep
%expanding the matrix.
for k=minGuess:NMax
    n=k;
    %Obtain approximate roots.
    lambda=sort(eig(A(1:k,1:k)),'descend');
    z=2./sqrt(lambda(1:N));%Approximate roots
    
    %The upper bound is from Equation 2.8 in [1]. Here, we use the
    %estimated zeros z in place of the true ones for computing the bound.
    JTemp1=besselj(m+2*n+2,z);
    JTemp2=besselj(m+1,z).^2;
    errBound=max(abs(besselj(m+2*n,z).*JTemp1)./(2*JTemp2*(m+2*n+1))+(JTemp1.^2./(2*JTemp2*(m+2*n+1)^2*(m+2*n))).*(z.^3/(8*pi)));
    if(errBound<epsVal||k==NMax)
        if(errBound>=epsVal)
            succeeded=false;
        else
            succeeded=true;
        end
        return;
    end
    
    %Increase the order by one.
    %The following equations are taken from Theorem 1.7 of [3].
    alphak=m+2*(k+1);
    dk=2/((alphak-1)*(alphak+1));
    fk=1/((alphak-1)*sqrt(alphak-2)*sqrt(alphak));
    A(k+1,k+1)=dk;
    A(k,k+1)=fk;
    A(k+1,k)=fk;
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
