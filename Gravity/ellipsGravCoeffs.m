function [CBar,SBar,a,c]=ellipsGravCoeffs(maxOrder,isNormalized,omega,a,f,GM)
%%ELLIPSGRAVCOEFFS Obtain the unnormalized or fully normalized spherical
%                  harmonic coefficients for the gravitational potential
%                  implied by an ellipsoidal approximation to the geoid
%                  with a given rotation rate.
%
%INPUTS: maxOrder The order of the coefficients. The full ellipsoidal
%                 gravity potential as expressed with spherical harmonic
%                 coefficients is an infinite series. If this parameter is
%                 omitted, then a default of maxOrder=2 is used.
%    isNormalized If true, then fully normalized coefficients are returned
%                 (those for use with fully normalized associated Legendre
%                 polynomials). Otherwise, unnormalized coefficients are
%                 returned (those for use with associated Legendre
%                 polynomials). If this parameter is omitted, then
%                 normalized coefficients will be returned.
%           omega The average rotation rate of the planet. If this argument
%                 is omitted, the value in Constants.WGS84EarthRotationRate
%                 is used, having units of radians per second.
%               a The semi-major axis of the reference ellipsoid. If this
%                 argument is omitted, the value in
%                 Constants.WGS84SemiMajorAxis is used, having units of
%                 meters.
%               f The flattening factor of the reference ellipsoid. If this
%                 argument is omitted, the value in
%                 Constants.WGS84Flattening is used (and is unitless).
%              GM The universal gravitational constant times the mass of
%                 the planet. If this parameter is omitted, then the value
%                 in Constants.WGS84GMWithAtmosphere is used, having units
%                 of m^3/s^2.
%
%OUTPUTS:  CBar An array such that if passed to the CountingClusterSet
%               class, CBar(n+1,m+1) is the unitless coefficient of the
%               cosine term of order m in a spherical harmonic expression
%               for the gravitational field that uses the associated
%               Legendre functions that are fully normalized if
%               isNormalized=true.
%          SBar An array holding the unitless coefficients of the sine
%               terms in a spherical harmonic expression for the
%               gravitational field that uses the associated Legendre
%               functions that are fully normalized if isNormalized=true.
%               All elements of SBar are zero due to the symmetry of the
%               ellipsoidal approximation.
%             a The numerator in the (a/r)^n term in the spherical harmonic
%               sum, having units of meters.
%           c   The constant value by which the spherical harmonic series
%               is multiplied, having units of m^3/s^2.
%
%This function implements a number of equations from Chapters 2.5, 2.7,
%2.8, and 2.9 of [1]. The equation numbers are cited in the code below.
%
%On the surface of the reference ellipsoid, the gravity potential V plus
%the potential due to the centripetal acceleration of the Earth, that is
%0.5*omega^2*(x^2+y^2), where [x;y;z] is the point on the ellipsoid, is a
%constant U0. U0 is the normal potential at the ellipsoid and is the value
%returned by the function ellipParam2Potential. These coefficients and the
%function ellipParam2Potential can be used to debug the function
%spherHarmonicEval.
%
%In the literature, reference is often made to a "J2" model for the Earth's
%gravity field. The J2 model is equivalent to the coefficients
%obtained when maxOrder=2. The J4 model is equivalent to the coefficients
%obtained when maxOrder=4.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6)
    GM=Constants.WGS84GMWithAtmosphere;
end

if(nargin<5)
    f=Constants.WGS84Flattening;
end

if(nargin<4)
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<3)
    omega=Constants.WGS84EarthRotationRate;
end

if(nargin<2)
    isNormalized=true;
end

if(nargin<1)
    maxOrder=2;
end

%The degree n ranges from 0 to maxOrder and the order m ranges from 0 to n.
%Thus, the number of elements for each n ranges from 1 to maxOrder+1. This
%means that a total of (maxOrder+1)*(maxOrder+2)/2 elements are necessary
%for the coefficients.
numPBarU=(maxOrder+1)*(maxOrder+2)/2;
totalP=zeros(numPBarU,1);
CBar=CountingClusterSet(totalP);
if(nargout>1)
    SBar=totalP; 
end

%In the ellipsoidal Earth approximation, the S coefficients are all zero
%and the only nonzero normalized C coefficients are the \bar{C}_{2n,0}
%terms.
b=a*(1-f);%The semi-major axis of the reference ellipsoid.
E=sqrt(a^2-b^2);%The linear eccentricity (Equation 2-174)

e=E/a;%The first numerical eccentricity (Equation 2-174)
epsilon=E/b;%The second numerical eccentricity (Equation 2-174)

%From Equation 2-137
m=omega^2*a^2*b/GM;
%From Equation 2-113
q0=0.5*((1+3*b^2/E^2)*atan(E/b)-3*b/E);

innerTerm=(1/3)*(1-(2/15)*(m*epsilon/q0));
for n=0:floor(maxOrder/2)
    %From Equation 2-170 and 2-167
    CBar(2*n+1,0+1)=(-1)^n*3*e^(2*n)/((2*n+1)*(2*n+3))*(1-n+5*n*innerTerm);
    if(isNormalized)
        %The normalization constant is described in Equation 2-80. It is
        %consistent with the normalization used in the function
        %spherHarmonicEval. It is not consistent with the optional 
        %normalization that Matlab's built-in legendre function uses.
        CBar(2*n+1,0+1)=CBar(2*n+1,0+1)/sqrt(2*(2*n)+1);
    end
end

%The second constant to return for the sum. 
c=GM;

CBar=CBar.clusterEls;

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
