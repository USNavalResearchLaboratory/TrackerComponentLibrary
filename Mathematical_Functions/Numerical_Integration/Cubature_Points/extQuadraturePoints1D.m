function [xiG,wG,xiK,wK,bounds]=extQuadraturePoints1D(n,algorithm,c1)
%%EXTQUADRATUREPOINTS1D Obtain quadrature points and weights to efficiently
%           numerically evaluate 1D integrals involving various weighting
%           functions. Unlike the function quadraturePoints1D, embedded
%           points and weights are obtained here (a Gauss-Kronrod
%           quadrature formula is created). That means that points and
%           weights of two different polynomial accuracies are returned,
%           whereby the higher accuracy points have the lower
%           accuracy points embedded in them (but different weights). This
%           is useful when performing adaptive integration to get an
%           estimate of the integration error while limiting the total
%           number of function evaluations necessary. For both sets of
%           points, The evaluation of a continuous integral using
%           quadrature points is
%           integral_lowL^upL w(x)*f(x) dx=sum_{i=1}^n w(i)*f(xi(i))
%           where the equality holds up for all polynomials up to a
%           particular degree. For high-order polynomials and other
%           functions, quadrature integration is just an approximation.
%
%%INPUTS: n A positive integer such that 2*n-1 is the highest degree to
%           which the lower-order quadrature points are accurate. The
%           lower-order quadrature points consist of n points, the higher-
%           order ones consists of 2*n+1 points with a polynomial accuracy
%           of 2*n+1.
% algorithm This specifies the type of weighting function and the range of
%           the integral. These points can generally be transformed for
%           integrals over other regions/ with other parameters. All of the
%           methods except 10 use parameters for the polynomial recursion
%           from Table 22.7 in Ch. 22 of [3], except where explicitly
%           stated. Possible values are:
%           0 The weighting function is w(x)=1/(sqrt(2*pi))*exp(-x^2/2),
%             the integration interval is (-Inf,Inf). Gauss-Kronrod points
%             exist for n=1 and n=2.
%           1 w(x)=exp(-x^2) on the interval (-Inf,Inf). Gauss-Kronrod
%             points exist for n=1 and n=2.
%           2 w(x)=1 on (-1,1). (The default if omittted or an empty matrix
%             is passed) 
%             Note that the transformation 
%             xiNew=xi*(b-a)/2+(b+a)/2;
%             wNew=w*(b-a)/2;
%             can be used to transform the points and weights to the
%             weighting function w(y)=1 on the range (a,b). Gauss-Kronrod
%             points should always exist.
%           3 w(x)=(1-x^2)^(c1-1/2) on (-1,1) with c1>-1/2. This becomes
%             numerically unstable for c1 close to zero.
%           4 w(x)=exp(-x) on (0,Inf). Gauss-Kronrod points exist for n=1.
%           5 w(x)=x^c1*exp(-x) on (0,Inf). Note that c1>-1. Gauss-Kronrod
%             points exist for n=1.
%           6 w(x)=(1+x)^c1 on (-1,1) for c1>-1. This is a special case of
%             the Jacobi polynomials. Gauss-Kronrod points exist for n=1.
%           7 w(x)=x^c1 on (0,1), c1>-1. The three-term recursion was
%             derived by explicitly evaluating the integrals with the
%             desired orthogonality constraints until a pattern could be
%             identified for an arbitrary order. See [2] for the
%             orthogonality constraint. Gauss-Kronrod points exist for n=1.
%           8 w(x)=|x|^c1 on (-1,1), c1>=0 and c1 is an integer. The three-
%             term recursion was derived by explicitly evaluating the
%             integrals with the desired orthogonality constraints until a
%             pattern could be identified for an arbitrary order. Separate
%             patterns for even and odd integers were found. Gauss-Kronrod
%             points exist for n=1,2,3.
%           9 w(x)=|x|^c1*exp(-x^2) on (-Inf,Inf), c1>=0 and c1 is an
%             integer. The three-term recursion was derived by explicitly
%             evaluating the integrals with the desired orthogonality
%             constraints until a pattern could be identified for an
%             arbitrary order. Separate patterns for even and odd integers
%             were found. Gauss-Kronrod points exist for n=1,2,3.
%          10 In this instance, the input c1 is a structure that provides
%             all of the necessary parameters for the recursion. See the
%             comments to c1 for more details.
%       c1 For values of algorithm not equal to 10, this input is a scalar
%          parameter in certain algorithms, as described above. For
%          algorithm 10, this parameter is a structure containing elements
%          that decibe the orthogonal polynomial associated with the
%          weighting function in the integral. The orthogonal polynomial
%          relation is discussed in more detail in the function
%          orthoPolyZerosFromRecur, which relies on the same recursion for
%          computing quadrature points and weights, albeit without a
%          Gauss-Kronrod extension. The elements of the structure are
%          a,b,c Function handles such that a(i), b(i) and c(i) give the
%              appropriate values for the recursive formulation of the
%              function values in
%              p_i(x)=(a(i)*x+b(i))*p_{(i-1)}(x)-c(i)*p_{i-2}(x)
%              where indexation starts at 1 and c(1) is not used.
%          mu0 The quadrature weights must be scaled by
%              mu0=integral_lowL^upL w(x) dx. This is that mu0.
%
%OUTPUTS: xiG The 1Xn set of quadrature points that with the weights wG
%             offer an order 2*n-1 polynomial accuracy for integration.
%          wG The nX1 set of weights associated with xiG. These are all
%             positive.
%         xiK The Gauss-Kronrod extension of xiG. The first n points are
%             those in xiG and the next n+1 points are different. This with
%             the weights wK have a polynomial accuracy of 2*n+1.
%          wG The weights of xiK. These are all different from the weights
%             in wG. These are all positive.
%      bounds The bounds over which the quadrature integrated. If
%             algorithm=10, this is just an empty matrix.
%
%Note that Gauss-Kronrod extensions of the points do not always exist for
%all n. For example, for algorithm 0, they only appear to exist for n=1 and
%n=2. If no extension exists, then an error is raised.
%
%This function implements the algorithm of [1]. To go from the general
%coefficients to the tri-diagonal coefficients needed for the algorithm,
%the diagonalization method of [2] is used.
%
%REFERENCES:
%[1] D. P. Laurie, "Calculation of Gauss-Kronrod quadrature rules,"
%    Mathematics of Computation, vol. 66, no. 219, pp. 1133-1145, Jul.
%    1997.
%[2] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%[3] Abramowitz, M. and Stegun, I. A. (Eds.). "Orthogonal Polynomials."
%    Table 22.7 in Ch. 22 in Handbook of Mathematical Functions with
%    Formulas, Graphs, and Mathematical Tables, 9th printing. New York:
%    Dover, pp. 782, 1972.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=2;
end

switch(algorithm)
    case 0%Weighting function w(x)=1/(sqrt(2*pi))*exp(-x^2/2)
          %on (-Inf,Inf)
        mu0=1;
        a=@(i)1;
        b=@(i)0;
        c=@(i)i-1;
        bounds=[-Inf;Inf];
    case 1%Weighting function w(x)=exp(-x^2) on (-Inf,Inf)
        mu0=sqrt(pi);
        a=@(i)2;
        b=@(i)0;
        c=@(i)2*(i-1);
        bounds=[-Inf;Inf];
    case 2%Weighting function w(x)=1 on (-1,1)
        mu0=2;
        a=@(i)(2*(i-1)+1)/i;
        b=@(i)0;
        c=@(i)(i-1)/i;
        bounds=[-1;1];
    case 3%Weighting function w(x)=(1-x^2)^(c1-1/2) on (-1,1) with c1>-1/2
          %and c1~=0. Given c1=0, Formula 3 is used.
        mu0=(sqrt(pi)*gamma((1/2)+c1))/gamma(1+c1);
        a=@(i)2*(i+c1-1)/i;
        b=@(i)0;
        c=@(i)(i+2*c1-2)/i;
        bounds=[-1;1];
    case 4%Weighting function w(x)=exp(-x) on (0,Inf)
        mu0=1;
        a=@(i)-1/i;
        b=@(i)(2*(i-1)+1)/i;
        c=@(i)(i-1)/i;
        bounds=[0;Inf];
    case 5%Weighting function w(x)=x^c1*exp(-x) on (0,Inf) for c1>-1.
        mu0=gamma(1+c1);
        a=@(i)(-1/i);
        b=@(i)(2*i-1+c1)/i;
        c=@(i)(i-1+c1)/i;
        bounds=[0;Inf];
    case 6%Weighting function w(x)=(1+x)^c1 on (-1,1) for c1>-1.
          %This uses a special case of the Jacobi polynomials.
        mu0=2^(c1+1)/(c1+1);
        a=@(i)(2*i+c1-1)*(2*i+c1)*(2*i+c1-2)/(2*i*(i+c1)*(2*i+c1-2));
        b=@(i)(2*i+c1-1)*(-c1^2)/(2*i*(i+c1)*(2*i+c1-2));
        c=@(i)2*(i-1)*(i+c1-1)*(2*i+c1)/(2*i*(i+c1)*(2*i+c1-2));
        bounds=[-1;1];
    case 7%Weighting function w(x)=x^c1 on (0,1) for c1>-1
        mu0=1/(1+c1);
        a=@(k)1;
        b=@(k)(-(c1^2+(2*(k-1)+1)*c1+2*k*(k-1))/((2*(k-1)+c1)*(2*k+c1)));
        c=@(k)((k-1)^2*(c1+k-1)^2/((2*k-3+c1)*(2*k-2+c1)^2*(2*k-1+c1)));
        bounds=[0;1];
    case 8%Weighting function w(x)=|x|^c1 on (-1,1), c1>=0 and c1 is an
           %integer.
        if(mod(c1,2)==0)%Even integer
            c1=c1/2;
            mu0=2/(1+2*c1);
            a=@(k)1;
            b=@(k)0;
            c=@(k)((k-1+2*c1*mod(k+1,2))^2/((2*k-3+2*c1)*(2*k-1+2*c1)));
        else%Odd integer
            c1=(c1-1)/2;
            mu0=1/(1+c1);
            a=@(k)1;
            b=@(k)0;
            c=@(k)(((k-mod(k,2))/2+mod(k+1,2)*c1)^2/((k-1+c1)*(k+c1)));
        end
        bounds=[-1;1];
    case 9%Weighting function w(x)=|x|^c1*exp(-x^2) on (-Inf,Inf), c1>=0
           %and c1 is an integer.
        if(mod(c1,2)==0)%Even integer
            c1=c1/2;
            mu0=gamma(c1+1/2);
            a=@(k)1;
            b=@(k)0;
            c=@(k)((k-1)/2+c1*mod(k+1,2));
        else%Odd integer
            c1=(c1-1)/2;
            mu0=factorial(c1);
            a=@(k)1;
            b=@(k)0;
            c=@(k)((k-mod(k,2))/2+c1*mod(k+1,2));
        end
        bounds=[-Inf;Inf];
    case 10%The user specified the functions.
        mu0=c1.mu0;
        a=c1.a;
        b=c1.b;
        c=c1.c;
        bounds=[];
    otherwise
        error('Unknown algorithm specified');
end

[alpha,beta]=getRecursionCoeffs(fix(3*n/2)+2,a,b,c,mu0);
[a,b]=getKonrodCoeffs(alpha,beta,n);
[xiK,wK,xiG,wG]=GaussKronrod4PolyRecur(n,a,b,mu0);

end

function [xiK,wK,xiG,wG]=GaussKronrod4PolyRecur(n,a,b,mu0)
%%GAUSSKRONROD4POLYRECUR Given the tri-diagonal weights of the Jacobi-
%               Kronrod matrix as in Equation 8 of [1], compute the
%               quadrature points and weights using the algorithm of [2].
%               In particular, both the full Gauss-Kronrod points as well
%               as the basic Gauss points and weights are computed. Since
%               the Gauss-Kronrod points are a subset of the basic points,
%               the values in the Gauss-Kronrod points are assigned to the
%               ones in the Gauss points (using a simple distance cost
%               function and assign2D) and are placed at the front of the
%               returned values.
%
%REFERENCES:
%[1] D. P. Laurie, "Calculation of Gauss-Kronrod quadrature rules,"
%    Mathematics of Computation, vol. 66, no. 219, pp. 1133-1145, Jul.
%    1997.
%[2] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTotal=2*n+1;
if(sum(b>0)~=numTotal)
    error('Gauss-Kronrod points for the given scenario do not exist.')
end

%Get the points and weights for the Gauss-Kronrod points.
J=zeros(numTotal,numTotal);
for k=1:(numTotal-1)
  J(k,k)=a(k);
  J(k,k+1)=sqrt(b(k+1));
  J(k+1,k)=J(k,k+1);
end
J(2*n+1,2*n+1)=a(2*n+1);
[V,D]=eig(J);
xiK=diag(D);
wK=abs(V(1,:)').^2*mu0;

%Get the points and weights for the Gauss points.
[V,D]=eig(J(1:n,1:n));
xiG=diag(D);
wG=abs(V(1,:)').^2*mu0;

%Sort the values so that in the assignment step, the corresponding values
%have the same ordering when they are placed to the front of xiK.
[xiK,idx]=sort(xiK);
wK=wK(idx);
[xiG,idx]=sort(xiG);
wG=wG(idx);

%We now want to determine which points in xiG correspond to those in xiK
%and separate them out. That way, one can evaluate the function just with
%the points in xiG and then only add the extension from xiK. We will just
%make a cost matrix and set up a 2D assignment problem.
C=zeros(n,numTotal);
for curG=1:n
   for curK=1:numTotal
       C(curG,curK)=abs(xiG(curG)-xiK(curK));
   end
end

k4G=assign2D(C,false);

sel=false(numTotal,1);
sel(k4G)=true;

xiK=[xiK(sel);xiK(~sel)]';
wK=[wK(sel);wK(~sel)];
xiG=xiG';

end

function [a,b]=getKonrodCoeffs(aInit,bInit,n)
%%GETKRONRODCOEFFS This implements the algorithm of [1] for a Gauss-Kronrod
%          extension based on the pseudocode given in the paper.

%REFERENCES:
%[1] D. P. Laurie, "Calculation of Gauss-Kronrod quadrature rules,"
%    Mathematics of Computation, vol. 66, no. 219, pp. 1133-1145, Jul.
%    1997.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTotal=2*n+1;

a=zeros(numTotal,1);
b=zeros(numTotal,1);

sel=0:floor(3*n/2);
a(sel+1)=aInit(sel+1);
sel=0:ceil(3*n/2);
b(sel+1)=bInit(sel+1);

halfFloorNum=floor(n/2)+2;
s=zeros(halfFloorNum,1);
t=zeros(halfFloorNum,1);
t(2)=b(n+2);
for m=0:(n-2)
  k=floor((m+1)/2):-1:0;
  l=m-k;
  s(k+2)=cumsum((a(k+n+2)-a(l+1)).*t(k+2)+b(k+n+2).*s(k+1)-b(l+1).*s(k+2));
  temp=s;
  s=t;
  t=temp;
end

j=floor(n/2):-1:0;
s(j+2)=s(j+1);
for m=(n-1):(2*n-3)
  k=m+1-n:floor((m-1)/2);
  l=m-k;
  j=n-1-l;
  s(j+2)=cumsum(-(a(k+n+2)-a(l+1)).*t(j+2)-b(k+n+2).*s(j+2)+b(l+1).*s(j+3));
  j=j(end);
  k=floor((m+1)/2);
  if(mod(m,2)==0)
      a(k+n+2)=a(k+1)+(s(j+2)-b(k+n+2)*s(j+3))/t(j+3);
  else
      b(k+n+2)=s(j+2)/s(j+3);
  end
  temp=s;
  s=t;
  t=temp;
end
a(numTotal)=a(n)-b(2*n+1)*s(2)/t(2);

end


function [alpha,beta]=getRecursionCoeffs(n,a,b,c,mu0)
%%GETRECURSIONCOEFFS The method of computing the Gauss-Kronrod quadrature
%           rule of [1] makes use of orthogonal polynomials with respect to
%           the weighting function w(x). This means that 
%           integral_lowL^upL w(x)*p_n(x)*p_m(x) dx =1 if m=n and is
%           zero otherwise. All such polynomials (p_n(x) for the nth
%           order) can be written in recursive form as
%              p_i(x)=(a(i)*x+b(i))*p_{(i-1)}(x)-c(i)*p_{i-2}(x)
%           This function takes the coefficients a, b, and c and transforms
%           them into two sets of coefficients as described for use in a
%           tri-diagonal matrix in [2], which corresponds to a recursion of
%           the form p_{i+1}(x)=(x-a(i))*p_i(x)-b_j*p_{j-1}(x)
%           where p_{-1}(x)=0, p_{0}(x)=1 and b_0 is defined to be
%           mu0=integral_lowL^upL w(x) dx
%
%REFERENCES:
%[1] D. P. Laurie, "Calculation of Gauss-Kronrod quadrature rules,"
%    Mathematics of Computation, vol. 66, no. 219, pp. 1133-1145, Jul.
%    1997.
%[2] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

alpha=zeros(n,1);
beta=zeros(n-1,1);

%Rows 1 through n-1
for i=1:(n-1)
    bi=b(i);
    ai=a(i);
    aiNext=a(i+1);
    ciNext=c(i+1);
    
    %From Equation 2.2 of [2].
    alpha(i)=-bi/ai;
    beta(i)=ciNext/(ai*aiNext);
end
%The final row
i=n;
bi=b(i);
ai=a(i);
alpha(n)=-bi/ai;

beta=[mu0;beta];
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
