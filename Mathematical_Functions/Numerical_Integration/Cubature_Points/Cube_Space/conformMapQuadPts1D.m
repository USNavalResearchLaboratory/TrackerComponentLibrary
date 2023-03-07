function [xi,w]=conformMapQuadPts1D(n,mapping,pointType,param)
%%CONFORMMAPQUADPTS1D Compute quadrature points and weights that are a
%           transformation of either Gauss-Legendre points or
%           Clenshaw-Curtis points for integration of a function over the
%           range -1 to 1. n+1 points are returned. Though the conformal
%           mappings used can sometimes produce worse results, as discussed
%           in [1], there are times when they can produce improved results
%           too.
%
%INPUTS: n The order of the polynomial approximation. n>=2.
%  mapping This selects which conformal mapping is applied to the points.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed and
%            pointType~=0) Use the asin Taylor series expansion in [1]. The
%            order of the expansion is given by param, which must be an
%            integer value and has a default of 9. 
%          1 Use the Kosloff-Tal-Ezer map given in [1] with
%            alpha=2/(param+1/param).
%          2 (The default if omitted or an empty matrix is passed and
%            pointType==0) Use the strip mapping given in [1]. This cannot
%            be used if pointType==1.
%          3 Use the approximate mapping given in the Appendix of [1].
% pointType A parameter specifying which quadrature points to use in the
%          mapping. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            points produced by GaussLegendrePoints1D.
%          1 Use the points produced by ClenshawCurtisPoints1D.
%    param A parameter for each of the mappings. For mapping=1,2,3, param
%          defaults to 1.4 and affects the mapping. For mapping=2,3
%          param>1. For mapping=1, param>0. For mapping=0, param is the
%          integer order of the Taylor series expansion (and should be
%          odd).
%
%This function implements the conformal mapping transformations discussed
%in [1]. The mapping tries to avoid the clustering of too many points near
%(-1,1), which is common when using Gauss_Legendre or Clenshaw-Curtis
%points.
%
%EXAMPLE:
%This is similar to one of the examples in [1]. Here, we compare the
%conformally mapped quadrature points integrating exp(-40*x.^2) with the 
% n=29;
% f=@(x)exp(-40*x.^2);
% 
% [xi,w]=conformMapQuadPts1D(n,1,0);
% numIntSol=sum(w'.*f(xi));
% analyticSol=(1/2)*sqrt(pi/10)*erf(2*sqrt(10));
% relErrorConform=(numIntSol-analyticSol)/abs(analyticSol)
% 
% [xi,w]=GaussLegendrePoints1D(n);
% numIntSolOrig=sum(w'.*f(xi));
% relErrorOrig=(numIntSolOrig-analyticSol)/abs(analyticSol)
%One will see that relErrorConform<relErrorOrig, demonstrating that the
%conformal mapping transofmration improves performance in this scenario.
%
%REFERENCES:
%[1] N. Hale and L. N. Trefethen, "New quadrature formulas from conformal
%    maps," SIAM Journal on Numerical Analysis, vol. 46, no. 2, pp.
%    930-948, 2008.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(pointType))
    pointType=0; 
end

if(nargin<2||isempty(mapping))
    if(pointType==0)
        mapping=2;
    else
        mapping=0;
    end 
end

if(nargin<4||isempty(param))
    if(mapping==0)
        param=9;
    else
        param=1.4;
    end
end

rho=param;

switch(pointType)
    case 0
        [xi,w]=GaussLegendrePoints1D(n);
    case 1
        [xi,w]=ClenshawCurtisPoints1D(n);
    otherwise
        error('Unknown Points specified.')
end

switch(mapping)
    case 0%Degree d asin expansion.
        d=rho;
        coeffs=zeros(1,d+1);
        coeffs(d:-2:1)=(1./(1:2:d)).*[1, cumprod(1:2:d-2)./cumprod(2:2:d-1)];
        coeffs=coeffs/sum(coeffs);
        derivCoeffs=polyder(coeffs);
        dg=polyval(derivCoeffs,xi);
        xi=polyval(coeffs,xi);
        
        %Transplant the quadrature from Equation 2.6
        w=w.*dg.';
    case 1%The Kosloff-Tal-Ezer map
        alpha=2/(rho+1/rho);
        
        dg=alpha./(sqrt(1-xi.^2*alpha^2)*asin(alpha));
        
        xi=asin(alpha*xi)./asin(alpha);

        %Transplant the quadrature from Equation 2.6
        w=w.*dg.';
    case 2%The strip mapping.
        if(pointType==1)
            error('The strip mapping (mapping=3) cannot be used with Clenshaw-Curtis points'); 
        end
        
        %For the sums in the numerator and denominator of Equation 3.1
        num=0;
        den=0;
        %Find m from rho. The total number of terms is just what was used
        %in [1] without derivation.
        for j=1:round((1/2)+sqrt(10/log(rho)))
            num=num+rho^(-4*(j-(1/2))^2);
            den=den+rho^(-4*j^2);
        end
        %Equation 3.1
        m4=2*num/(1+2*den);%m^(1/4)
        m=m4^4;
        
        K=ellipIntInc1Kind(pi/2,m);
        u=asin(xi);
        omega=2*K*u/pi;
        [sn,cn,dn]=JacobiEllipticFuncVals(omega,m);

        %Equation 3.3
        dg=(2*K*m4./(pi*sqrt(1-xi.^2))).*(cn.*dn./(1-m4^2*sn.^2))/atanh(m4);
        
        %Equation 3.2
        xi=atanh(m4*sn)/atanh(m4);
        
        %Transplant the quadrature from Equation 2.6
        w=w.*dg.';
    case 3
        u=asin(xi);
        tau=pi/log(rho);

        %Terms in A.1 and A.2
        d=(1/2)+1/(exp(tau*pi)+1);
        pd2pu=pi/2+u;
        pd2mu=pi/2-u;

        C=1/(log(1+exp(-tau*pi))-log(2)+(pi/2)*tau*d);

        %Equation A.2
        dg=-C*(tau./sqrt(1-xi.^2)).*(1./(exp(tau*pd2pu)+1)+1./(exp(tau*pd2mu)+1)-d);
        %Equation A.3
        dg([1,end])=(C*tau^2/4)*tanh((pi/2)*tau)^2;

        %Transplant the quadrature from Equation 2.6
        w=w.*dg.';

        %Equation A.1 for g. We direct insert xi as in Equaion 2.6
        xi=C*(log(1+exp(-tau*pd2pu))-log(1+exp(-tau*pd2mu))+d*tau*u);
    otherwise
        error('Unknown Mapping specified.')
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
