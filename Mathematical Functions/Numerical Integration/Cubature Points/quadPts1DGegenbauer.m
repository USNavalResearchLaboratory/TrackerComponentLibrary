function [xi,w]=quadPts1DGegenbauer(lambda,n,algorithm)
%%QUADPTS1DGEGENBAUER Obtain points for numeric quadrature so that one can
%               approximate the integral
%               int_{-1}^1 f(x)*(1-x^2)^(lambda-1/2) dx
%               with the sum(f(xi).*w.'). The term (1-x^2)^(lambda-1/2) is
%               the weight function and is orthogonal when integrated over
%               (-1,1) times Gegenbauer polynomials.
%
%INPUTS: lambda The real lambda term in the Gegenbauer polynomial. lambda
%               >-1/2.
%             n The order of the approximation. There will be n+1 points.
%     algorithm An optional parameter specifying the algorithm used.
%               Possible values are:
%               0 (The default if omitted or an empty matrix is passed).
%                 Use GegenbauerPolyChebMoments to get the moments and pass
%                 the result to FejerPtsWeighted1D.
%               1 USe the algorithm of [1].           
%
%OUTPUTS: xi A 1XnumPoints vector containing the quadrature points.
%           w A numPointsX1 vector of the weights associated with the
%             quadrature points.
%
%EXAMPLE:
%We compare a numeric value to an exact solution.
% n=32;
% lambda=pi;
% [xi,w]=quadPts1DGegenbauer(lambda,n,1);
% numIntVal=sum(xi.^(24).*w.');
% exactVal=(316234143225*sqrt(pi)*gamma(1/2+pi))/(4096*gamma(13+pi));
% relErr=(numIntVal-exactVal)/abs(exactVal)
%One will see that the relative error is near eps(1) indicating that there
%is a very good integral approximation.
%
%REFERENCES:
%[1] D.B. Hunter, H.V. Smith, A quadrature formula of Clenshaw-Curtis type
%    for the Gegenbauer weight-function, Journal of Computational and
%    Applied Mathematics 177 (2005) 389-400.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0
        moments=zeros(n+1,1);
        k=fix(n/2);
        moments((0:2:(2*k))+1)=GegenbauerPolyChebMoments(lambda,k);
        
        [xi,w]=FejerPtsWeighted1D(moments);
    case 1
        %Equation 1.12 in [1]
        xi=cos((0:n)*pi/n);

        s=floor(n/2);

        I=GegenbauerPolyChebMoments(lambda,s).';
        I(1)=I(1)/2;
        if(mod(n,2)==0)
            I(s)=I(s)/2; 
        end

        %We get w from Equation 1.13.
        w=zeros(n+1,1);
        for i=0:n
           w(i+1)=(2/n)*sum(I.*cos(2*i*(0:s)*pi/n));
        end
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
