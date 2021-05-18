function [PBarVals,dPBarValsdTheta,d2PBarValsdTheta2]=LegendreCos(theta,M,scalFactor)
%%LEGENDRECOS Evaluate the fully associated Legendre functions of
%           cos(theta) of degree n and order m for all n from 0 to M and
%           for each n, m goes from 0 to n. Such a function is typically
%           written as \bar{P}_{nm}(cos(theta)). Also evaluate
%           D{\bar{P}_{nm}(cos(theta))} and D2{\bar{P}_{nm}(cos(theta))},
%           where D{} is the first derivative with respect to theta and
%           D2{} is the second derivative with respect to theta. All of the
%           values can be scaled by a factor of scalFac, if desired, to
%           help prevent overflows with high degrees and orders.
%
%INPUTS: theta An angle in radians.
%       maxDeg The maximum degree and order of the output. This must be
%              >=0.
%   scalFactor A scale factor to help prevent overflow of the results. If
%              this parameter is omitted or an empty matrix is passed, the
%              default of 1 is used.
%
%OUTPUTS: PBarUVals An instance of the CountingClusterSet class such that
%                   PBarVals(n+1,m+1)=scalFac*\bar{P}_{nm}(cos(theta)).
%                   To extract all coefficients as a vector just call
%                   PBarUVals(:).
%   dPBarValsdTheta An instance of the CountingClusterSet class such that
%                   dPBarValsdTheta(n+1,m+1)=scalFac*D{\bar{P}_{nm}(cos(theta))}
% d2PBarValsdTheta2 An instance of the CountingClusterSet class such that
%                   d2PBarValsdTheta2(n+1,m+1)=scalFac*D2{\bar{P}_{nm}(cos(theta))}
%
%The definition of associated Legendre polynomial underlying this function
%is
% P_(nm)(x)=(1-x^2)^(m/2)*Dm{P(n,x)}
%where Dm{} represents the mth-order derivative and P(n,x) is the Legendre
%polynomial having the form
%P(n,x)=(1/2^n)*sum_{k=0}^n binomial(n,n)*(x-1)^(n-k)*(x+k)^k
%This definition lacks the Condon-Shortley phase. Adding the phase is the
%same as multiplying P_(nm)(x) by (-1)^m. The polynomials are fully
%normalized, meaning that each one is multiplied by
%N(n,m)=sqrt((2*n+1)*(2-KDelta(m))*factorial(n-m)/factorial(n+m))
%This is normalization 0 in the changeSpherHarmonicNorm function.
%
%Fully normalized associated Legendre polynomials often arise when
%computing spherical harmonic expansions. When simply evaluating an
%expansion, one will typically use fully normalized Legendre functions that
%have been divided by abs(sin(theta))^m as in [1]. However, when fitting
%coefficients for such an expansion, one will typically need the function
%values directly, which is what this function provides.
%
%This function just calls the function NALegendreCosRat and then multiplies
%out the denominators.
%
%EXAMPLE:
%One can see that this function returns the same values as the legendre
%function of cos(theta) when the proper normalization is selected. However,
%this function returns more values than legendre, as Matlab's legendre
%function just returns a vector for a single value of n.
% M=3;
% theta=-0.5;
% PBarVecMatLab=legendre(M,cos(theta));
% PBarVals=LegendreCos(theta,M);
% arePolyVals=true;
% origNormType=0;
% normType=8;
% PBarVals=changeSpherHarmonicNorm(PBarVals,origNormType,normType,arePolyVals);
% PBarVec=PBarVals(M+1,:);
% ratioCompare=PBarVecMatLab./PBarVec
%One will see that the ratio is essentially all ones.
%
%REFERENCES:
%[1] S. A. Holmes and W. E. Featherstone, "A unified approach to the
%    Clenshaw summation and the recursive computation of very high degree
%    and order normalised associated Legendre functions," Journal of
%    Geodesy, vol. 76, no. 5, pp. 279-299, May 2002.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(scalFactor))
    scalFactor=1; 
end

u=abs(sin(theta));

if(nargin<2)
    PBarVals=NALegendreCosRat(theta,M,scalFactor);
    for n=1:M
        um=u;
        for m=1:n
            PBarVals(n+1,m+1)=PBarVals(n+1,m+1)*um;
            um=um*u;
        end
    end
elseif(nargout==2)
    [PBarVals,dPBarValsdTheta]=NALegendreCosRat(theta,M,scalFactor);
    for n=1:M
        um=u;
        for m=1:n
            PBarVals(n+1,m+1)=PBarVals(n+1,m+1)*um;
            dPBarValsdTheta(n+1,m+1)=dPBarValsdTheta(n+1,m+1)*um;
            um=um*u;
        end
    end
else%nargout==3
    [PBarVals,dPBarValsdTheta,d2PBarValsdTheta2]=NALegendreCosRat(theta,M,scalFactor);
    for n=1:M
        um=u;
        for m=1:n
            PBarVals(n+1,m+1)=PBarVals(n+1,m+1)*um;
            dPBarValsdTheta(n+1,m+1)=dPBarValsdTheta(n+1,m+1)*um;
            d2PBarValsdTheta2(n+1,m+1)=d2PBarValsdTheta2(n+1,m+1)*um;
            um=um*u;
        end
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
