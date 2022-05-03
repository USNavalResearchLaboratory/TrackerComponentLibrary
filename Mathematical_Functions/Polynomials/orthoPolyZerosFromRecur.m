function [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0)
%%ORTHOPOLYZEROSFROMRECUR Consider a family of polynomials p_n(x) such that
%              p_n(x) is nth order and the polynomials are orthonormal with
%              respect to some weighting function w(x). This means that
%              integral_lowL^upL w(x)*p_n(x)*p_m(x) dx =1 if m=n and is
%              zero otherwise. All such polynomials can be written in
%              recursive form as
%              p_i(x)=(a(i)*x+b(i))*p_{(i-1)}(x)-c(i)*p_{i-2}(x)
%              starting with p_{-1}(x)=0 and p_0(x)=1. This function finds
%              the zeros of the nth polynomial given expressions for the
%              recursion of the polynomials. Parameters for common
%              orthogonal polynomials are given below. The zeros of such
%              orthonormal polynomials can be used as quadrature points for
%              1D integration from lowL to upL of all polynomials up to
%              order 2*n-1 times the weighting function w(x). The
%              associated quadrature weights for such integrals can also be
%              returned by this function if two outputs are requested.
%
%INPUTS: n The order of the polynomial whose zeros are desired.
%    a,b,c There are two formulations for the input. In the first one,
%          a,b, and c are vectors or function handles such that a(i), b(i),
%          and c(i) give the appropriate values for the recursive
%          formulation of the function values in
%          p_i(x)=(a(i)*x+b(i))*p_{(i-1)}(x)-c(i)*p_{i-2}(x)
%          where the values of i used in a and b go from 1 to n, and
%          indices used in  c(i) go from 2 to n. c(1) is not used.
%          Alternatively, if c is omitted or an empty matrix is passed,
%          then a and b are assumed to be coefficients in the three-term
%          expansion of the form
%          b(i)*p_i(x)=(x-a(i))*p_{(i-1)}(x)-b(i-1)*p_{i-2}(x)
%          and can be vectors or function handles. Indices of a run from 1
%          to n and b run from 1 to (n-1).
%      mu0 This parameter is only needed if the second output of this
%          function is desired. The quadrature weights must be scaled by
%          mu0=integral_lowL^upL w(x) dx. This is that mu0.
%
%OUTPUTS: xi The 1Xn set of zeros of the given set of orthonormal
%            polynomials.
%          w The set of n cubature weights associated with the orthonormal
%            polynomials. This requires the input mu0 be specified.
%
%When xi and w are used for quadrature integration, one can write
%integral_lowL^upL w(x)*f(x) dx=sum_{i=1}^n w(i)*f(xi(i))
%where the equality holds for all polynomials up to degree 2*n-1. This is
%the fundamental theorem of Gaussian quadrature. Thus, finding the zeros of
%polynomials that are orthonormal over a region of interest with respect to
%a particular weighting function can be useful.
%
%The algorithm is that of [1]. Because an eigenvalue decomposition is used,
%the results can be numerically unstable for large n.
%
%Examples of values for common polynomials, are given in [2]. Some are:
%1) Hermite polynomials. These are orthonormal with respect to exp(-x^2)
%   integrated over (-Inf,Inf). The recursion from [3] is
% mu0=sqrt(pi);
% a=@(i)2;
% b=@(i)0;
% c=@(i)2*(i-1);
%
%2) Legendre polynomials. These are orthonormal with respect to w(x)=1 for
%   integration over [-1,+1]. The recursion from [3] is
% mu0=2;
% a=@(i)(2*(i-1)+1)/i;
% b=@(i)0;
% c=@(i)(i-1)/i;
%However, the function GaussLegendrePoints1D should provide
%higher-precision zeros than this function.
%
%3) Laguerre polynomials. These are orthonormal with respect to
%   w(x)=exp(-x) for integration over [0,Inf). The recursion from [3] is
% mu0=1;
% a=@(i)-1/i;
% b=@(i)(2*(i-1)+1)/i;
% c=@(i)(i-1)/i;
%
%4) Associated Laguerre polynomials. These are orthonormal with respect to
%   w(x)=x^c1*exp(-x) for integration over [0,Inf). When considering
%   c1>-1, the recursion from [3] is
% mu0=gamma(1+c1);
% a=@(i)(-1/i);
% b=@(i)(2*i-1+c1)/i;
% c=@(i)(i-1+c1)/i;
%
%5) Gegenbauer polynomials. These are orthonormal with respect to
%   w(x)=(1-x^2)^(c1-1/2) for integration over [-1,1]. When considering
%   c1>-(1/2), the recursion from [3] is
% mu0=(sqrt(pi)*gamma((1/2)+c1))/gamma(1+c1);
% a=@(i)2*(i+c1-1)/i;
% b=@(i)0;
% c=@(i)(i+2*c1-2)/i;
%
%6) Modified Hermite polynomials that are orthogonal with respect to the
%   standard normal N(0,1) distribution integrated over (-Inf,Inf). These
%   are for w(x)=1/sqrt(2*pi)*exp(-x^2/2). The recursion [3] is
% mu0=1;
% a=@(i)1;
% b=@(i)0;
% c=@(i)i-1;
%
%REFERENCES:
%[1] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%[2] Weisstein, Eric W. "Orthogonal Polynomials." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/OrthogonalPolynomials.html
%[3] Abramowitz, M. and Stegun, I. A. (Eds.). "Orthogonal Polynomials."
%    Table 22.7 in Ch. 22 in Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables, 9th printing. New York:
%    Dover, pp. 782, 1972.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

J=zeros(n,n);
if(nargin<4||isempty(c))
    %If alpha and beta are directly given as a and b.
    for i=1:(n-1)
        J(i,i)=a(i);
        J(i,i+1)=b(i);
        J(i+1,i)=b(i);
    end
    J(n,n)=a(n);
else
    %The first row
    %Rows 1 through n-1
    for i=1:(n-1)
        bi=b(i);
        ai=a(i);
        aiNext=a(i+1);
        ciNext=c(i+1);

        %From Equation 2.2 in [1].
        alpha=-bi/ai;
        beta=sqrt(ciNext/(ai*aiNext));

        J(i,i)=alpha;
        J(i,i+1)=beta;
        J(i+1,i)=beta;
    end
    %The final row
    i=n;
    bi=b(i);
    ai=a(i);
    alpha=-bi/ai;
    J(i,i)=alpha;
end

[V,D]=eig(J);

xi=diag(D)';
w=abs(V(1,:)').^2*mu0;

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
