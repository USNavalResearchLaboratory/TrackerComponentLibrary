function [xi,w]=quadraturePoints1D(n,algorithm,c1)
%%QUADRATUREPOINTS1D Obtain quadrature points and weights to efficiently
%           numerically evaluate 1D integrals involving various weighting
%           functions. The quadrature points and weights are based off
%           properties of orthogonal polynomials. The evaluation of a
%           continuous integral using cubature points is
%           integral_lowL^upL w(x)*f(x) dx=sum_{i=1}^n w_i*f(xi(i))
%           where the equality holds up for all polynomials up to a
%           particular degree. For high-order polynomials and other
%           functions, quadrature integration is just an approximation.
%
%INPUTS: n A positive integer such that 2n-1 is the highest degree to
%          which the quadrature points are accurate. n points are returned
%          by the function.
% algorithm This specifies the type of weighting function and the range of
%          the integral. These points can generally be transformed for
%          integrals over other regions/ with other parameters. Possible
%          values are
%          0 (The default if omitted or an empty matrix is passed) The
%            weighting function is w(x)=1/(sqrt(2*pi))*exp(-x^2/2), the
%            integration interval is (-Inf,Inf). The algorithm of [2] is
%            used with parameters from Table 22.7 in Ch. 22 of [1].
%          1 w(x)=exp(-x^2) on the interval (-Inf,Inf). The algorithm of
%            [2] is used with parameters from Table 22.7 in Ch. 22 of [1].
%          2 w(x)=1 on (-1,1). The function GaussLegendrePoints1D is used.
%            Note that the transformation 
%            xiNew=xi*(b-a)/2+(b+a)/2;
%            wNew=w*(b-a)/2;
%            can be used to transform the points and weights to the
%            weighting function w(y)=1 on the range (a,b).
%          3 w(x)=(1-x^2)^(c1-1/2) on (-1,1) with c1>-1/2. The algorithm of
%            [2] is used with the values from Table 22.7 of [1]. Formula 10
%            is used for the c1=0 case. This becomes numerically unstable
%            for c1 close to but not equal to zero.
%          4 w(x)=exp(-x) on (0,Inf). The algorithm of [2] is used with the
%            values from Table 22.7 of [1].
%          5 w(x)=x^c1*exp(-x) on (0,Inf). The algorithm of [2] is used
%            with the values from Table 22.7 of [1] for c1>-1.
%          6 w(x)=(1+x)^c1 on (-1,1) for c1>-1. This is a special case of
%            the Jacobi polynomials. The algorithm of [2] is used with
%            values from Table 22.7 of [1].
%          7 w(x)=x^c1 on (0,1), c1>-1. The algorithm of [2] is used. The
%            three-term recursion was derived by explicitly evaluating
%            the integrals with the desired orthogonality constraints
%            until a pattern could be identified for an arbitrary order.
%          8 w(x)=|x|^c1 on (-1,1), c1>=0 and c1 is an integer. The
%            algorithm of [2] is used. The three-term recursion was derived
%            by explicitly evaluating the integrals with the desired
%            orthogonality constraints until a pattern could be identified
%            for an arbitrary order. Separate patterns for even and odd
%            integers were found.
%          9 w(x)=|x|^c1*exp(-x^2) on (-Inf,Inf), c1>=0 and c1 is an
%            integer. The algorithm of [2] is used. The three-term
%            recursion was derived by explicitly evaluating the integrals
%            with the desired orthogonality constraints until a pattern
%            could be identified for an arbitrary order. Separate patterns
%            for even and odd integers were found.
%         10 w(x)=1/sqrt(1-x^2) on (-1,1). The formula in 25.4.38 in [1] is
%            used. Note that the transformation
%            xiNew=(b+a)/2+xi*(b-a)/2; can be used to change the
%            quadrature points to the weighting function
%            w(y)=1/sqrt((y-a)*(b-y)) on (a,b).
%         11 w(x)=sqrt(1-x^2) on (-1,1). The formula in 25.4.40 in [1] is
%            used. Note that the transformation xiNew=(b+a)/2+xi*(b-a)/2
%            can be used to change the quadrature points to the weighting
%            function w(y)=sqrt((y-a)*(b-y)) on (a,b).
%         12 w(x)=sqrt(x/(1-x)) on (0,1). The formula in 25.4.42 in [1] is
%            used. Note that the transformation xiNew=a+(b-a)*xi can be
%            used to change the quadrature points to the weighting
%            function w(x)=sqrt((x-a)/(b-x)) on (a,b).
%      c1 This is a parameter that is only used with certain algorithms
%         described above.
%
%OUTPUTS: xi A 1 X n vector containing the quadrature points.
%          w An nX1 vector of the weights associated with the quadrature
%            points. Depending on the algorithm, these may or may not sum
%            to one.
%
%In the algorithms implemented in this function, with the exception of 11
%and 12, this function obtains points by passing parameters for a three-
%point recursion to the orthoPolyZerosFromRecur function.
%
%Note that the function linCubPoints2MultiDim can be given a handle to this
%function to produce multi-dimensional cubature formula.
%
%REFERENCES:
%[1] Abramowitz, M. and Stegun, I. A. (Eds.). "Orthogonal Polynomials."
%    in Ch. 22, 25 in Handbook of Mathematical Functions with Formulas,
%    Graphs, and Mathematical Tables, 9th printing. New York: Dover, 1972.
%[2] G. H. Golub and J. H. Welsh, "Calculation of Gauss quadrature rules,"
%    Mathematics of Computation, vol. 23, pp. 221-230, 1969.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0%Weighting function w(x)=1/(sqrt(2*pi))*exp(-x^2/2)
          %on (-Inf,Inf)
        mu0=1;
        a=@(i)1;
        b=@(i)0;
        c=@(i)i-1;
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 1%Weighting function w(x)=exp(-x^2) on (-Inf,Inf)
        mu0=sqrt(pi);
        a=@(i)2;
        b=@(i)0;
        c=@(i)2*(i-1);
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 2%Weighting function w(x)=1 on (-1,1)
        [xi,w]=GaussLegendrePoints1D(n);
    case 3%Weighting function w(x)=(1-x^2)^(c1-1/2) on (-1,1) with c1>-1/2
          %and c1~=0. Given c1=0, Formula 10 is used.
        if(c1==0)
            [xi,w]=quadraturePoints1D(n,10);
            return;
        end
          
        mu0=(sqrt(pi)*gamma((1/2)+c1))/gamma(1+c1);
        a=@(i)2*(i+c1-1)/i;
        b=@(i)0;
        c=@(i)(i+2*c1-2)/i;
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 4%Weighting function w(x)=exp(-x) on (0,Inf)
        mu0=1;
        a=@(i)-1/i;
        b=@(i)(2*(i-1)+1)/i;
        c=@(i)(i-1)/i;
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 5%Weighting function w(x)=x^c1*exp(-x) on (0,Inf) for c1>-1.
        mu0=gamma(1+c1);
        a=@(i)(-1/i);
        b=@(i)(2*i-1+c1)/i;
        c=@(i)(i-1+c1)/i;
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 6%Weighting function w(x)=(1+x)^c1 on (-1,1) for c1>-1.
          %This uses a special case of the Jacobi polynomials.
        mu0=2^(c1+1)/(c1+1);
        a=@(i)(2*i+c1-1)*(2*i+c1)*(2*i+c1-2)/(2*i*(i+c1)*(2*i+c1-2));
        b=@(i)(2*i+c1-1)*(-c1^2)/(2*i*(i+c1)*(2*i+c1-2));
        c=@(i)2*(i-1)*(i+c1-1)*(2*i+c1)/(2*i*(i+c1)*(2*i+c1-2));
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 7%Weighting function w(x)=x^c1 on (0,1) for c1>-1
        mu0=1/(1+c1);
        a=@(k)1;
        b=@(k)(-(c1^2+(2*(k-1)+1)*c1+2*k*(k-1))/((2*(k-1)+c1)*(2*k+c1)));
        c=@(k)((k-1)^2*(c1+k-1)^2/((2*k-3+c1)*(2*k-2+c1)^2*(2*k-1+c1)));
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
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
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
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
        [xi,w]=orthoPolyZerosFromRecur(n,a,b,c,mu0);
    case 10%Weighting function w(x)=1/sqrt(1-x^2) on (-1,1) from 25.4.38 in
          %[1].
        xi=cos((2*(1:n)-1)*pi/(2*n));
        w=repmat(pi/n,[n,1]);
    case 11%Weighting function w(x)=sqrt(1-x^2) on (-1,1) from 25.4.40 in
          %[1].
        xi=cos((1:n)*pi/(n+1));
        w=(pi/(n+1))*sin((1:n)'*pi/(n+1)).^2;
    case 12%Weighting function w(x)=sqrt(x/(1-x)) on (0,1) from 25.4.42 in
          %[1].
        xi=cos((pi/2)*((2*(1:n)-1)/(2*n+1))).^2;
        w=(2*pi/(2*n+1))*xi';
    otherwise
        error('Unknown algorithm specified');
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
