function [xi,w]=FejerPoints1D(n,rule)
%%FEJERPOINTS1D Compute one of two types of Fejér quadrature points and
%          weights for integrating polynomials (or approximating integrals
%          of non-polynomials) over the range -1 to 1. Fejér's first rule
%          produces points n points allowing for an order n-1 quadrature
%          rule and Fejér's second rules produces n-1 points allowing for
%          an order n-2 quadrature rule. Contrast with with
%          GaussLegendrePoints1D, where  n+1 points can integration all
%          polynomials up to order 2*n+1 using n points. However, the Fejér
%          points perform well with non-polynomial functions, so they are
%          not necessarily worse than Gauss-Legendre points.
%
%INPUTS: n The parameter related to the number of points produced. For
%          rule=1, n points are produced and for rule=2, n-1 points are
%          produced.
%     rule A parameter specifying the Fejér interpolation rule to use.
%          The rules are described in [1]. Possible values are:
%          1 (The default if omitted or an empty matrix is passed) Use
%             Fejér's first rule.
%          2 Use Fejér's second rule.
%
%%OUTPUTS: xi A 1XnumPoints vector containing the quadrature points.
%          w A numPointsX1 vector of the weights associated with the
%            quadrature points. Note that these sum to 2, not 1, and are
%            all positive.
%
%The algorithms of [1] are implemented.
%
%EXAMPLE:
%Here, we integrate x.^2.*exp(-x.^2-x/2) from -1 to 1 using the Fejér
%quadrature points. We then compare the result to an analytic solution.
% [xi,w]=FejerPoints1D(30);
% f=@(x)(x.^2.*exp(-x.^2-x/2));
% numIntSol=sum(w'.*f(xi));
% analyticSol=(-12-20*exp(1)+9*exp(25/16)*sqrt(pi)*(erf(3/4)+erf(5/4)))/(32*exp(3/2));
% 
% relError=(numIntSol-analyticSol)/abs(analyticSol)
%One can see that the relative error is om the order of eps(1), so the
%analytic and quadrature solutions agree well.
%
%REFERENCES:
%[1] J. Waldvogel, "Fast construction of the Féjer and Clenshaw-Curtis
%    quadrature rules," BIT Numerical Mathematics, vol. 46, no. 1, pp. 195-
%    202, Mar. 2006.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(rule))
    rule=1;
end

switch(rule)
    case 1
        %Equation 4.4
        v=zeros(n,1);
        nm1d2Fix=fix((n-1)/2);
        k=0:nm1d2Fix;
        v(1:(nm1d2Fix+1))=(2./(1-4*k.^2)).*exp(1j*k*pi/n);
        %This starting index skips index n/2 when n is even, meaning that
        %the item in v(n/2) remains zero in that special case.
        startIdx=nm1d2Fix+2+(mod(n,2)==0);
        v(startIdx:n)=conj(v((nm1d2Fix+1):-1:2));
        w=ifft(v);%Fejer's first rile.

        xi=zeros(1,n);
        xi(1)=cos((1/2)*pi/n);
        
        %We use the angle addition formula to fill in the rest of the
        %values of cos((k-1/2)*pi/n), which are the points in x. We are
        %adding pi/n to the angle each time.
        sinThetaCur=sin((1/2)*pi/n);
        cosTheta=cos(pi/n);
        sinTheta=sin(pi/n);
        for k=2:n
            xi(k)=cosTheta*xi(k-1)-sinTheta*sinThetaCur;
            sinThetaCur=sinTheta*xi(k-1)+sinThetaCur*cosTheta;
        end
    case 2
        %We are creating v as in Equation 3.10
        n2Fix=fix(n/2);
        v=zeros(n,1);
        k=(0:1:(n2Fix-1));
        v(1:(n2Fix+1))=[2./(1-4*k.^2),(n-3)/(2*n2Fix-1)-1];
        startIdx=n2Fix+mod(n,2);
        v((n2Fix+2):n)=v(startIdx:-1:2);
        w=ifft(v);%Fejer's second rule.
        
        %The first point should always have zero weight, so we drop it.
        w=w(2:end);

        xi=zeros(1,n-1);
        cosTheta=cos(pi/n);
        %A double-angle formula.
        cos2Theta=2*cosTheta^2-1;
        xi(1)=cosTheta;
        xi(2)=cos2Theta;
        %We use the Chebyshev method to fill in the rest of the values of
        %cos(k*pi/n), which are the points in x.
        for k=3:(n-1)
            xi(k)=2*cosTheta*xi(k-1)-xi(k-2);
        end
    otherwise
        error('Unknown rule specified.')
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
