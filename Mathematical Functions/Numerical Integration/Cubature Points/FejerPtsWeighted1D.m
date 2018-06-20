function [xi,w]=FejerPtsWeighted1D(moments,rule)
%%FEJERPTSWEIGHTED1D Compute one of two types of Fejér quadrature points
%          and weights for integrating polynomials (or approximating
%          integrals of non-polynomials) over the range -1 to 1 times a
%          weighting function w(x) whose moments are given with respect to
%          Chebyshev polynomials of the first kind from order 0 to n are
%          given by the input moments. This function is more general than
%          FejerPoints1D, but requires the user to provide the moments.
%
%INPUTS: moments A length n+1 vector containing values of the integral
%                int_{-1}^1 T_k(x) w(xx) dx for k=0 to n, where w(x) is the
%                weighting function over which integration using these
%                quadrature points should be performed.
%           rule A parameter specifying the Fejér interpolation rule to use.
%                The rules are described in [1]. Possible values are:
%                1 (The default if omitted or an empty matrix is passed)
%                  Use Fejér's first rule.
%                2 Use Fejér's second rule.
%
%OUTPUTS: xi A 1XnumPoints vector containing the quadrature points.
%          w A numPointsX1 vector of the weights associated with the
%            quadrature points.
%
%This function implements the algorithm described in [1].
%
%EXAMPLE:
%Here, we generate a formula for integration over Gegenbauer polynomials,
%using GegenbauerPolyChebMoments to get the necessary moments of the
%polynomials with respect to Chebyshev polynomials of the first kind.
% k=16;
% lambda=pi;
% %Note that odd moments are all zero.
% moments((0:2:(2*k))+1)=GegenbauerPolyChebMoments(lambda,k);
% [xi,w]=FejerPtsWeighted1D(moments);
% 
% %We compare the integral of x^24 times the Gegenbauer polynomial using the
% %cubature points to an analytic solution.
% numIntVal=sum(xi.^(24).*w.');
% exactVal=(316234143225*sqrt(pi)*gamma(1/2+pi))/(4096*gamma(13+pi));
% relErr=(numIntVal-exactVal)/abs(exactVal)
%One will see that the relative error is near eps(1) indicating that there
%is a very good integral approximation.
%
%REFERENCES:
%[1] A. Sommariva, "Fast construction of Fejér and Clenshaw-Curtis rule for
%    general weight functions," Computers and Mathematics with
%    Applications, vol. 65, no. 4, pp. 682-693, Feb. 2013.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(rule))
    rule=1;
end

moments=moments(:);

n=length(moments);
switch(rule)
    case 1
        %Equation 2 in [1].
        xi=cos((2*(1:n)-1)*pi/(2*n));
        %Equation 7 in [1], modified to have the proper normalization.
        w=(1/n)*discSinCosTrans(moments,'CIIIe');
    case 2
        %The first equation in Section 3.
        theta=(1:n)*pi/(n+1);
        xi=cos((1:n)*pi/(n+1));
        
        %From remark 1 in Section 3, we convert the moments as required.
        lambda=ChebyshevI2IIMoments(moments);
      
        %Equation 21, , modified to have the proper normalization.
        w=(sin(theta)/(n+1)).*discSinCosTrans(lambda,'SIe').';
        w=w';
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
