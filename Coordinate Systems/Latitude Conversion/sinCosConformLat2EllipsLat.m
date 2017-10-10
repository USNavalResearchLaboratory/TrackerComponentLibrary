function phi=sinCosConformLat2EllipsLat(sinChi,cosChi,f)
%%SINCOSCONFORMLAT2ELLIPSLAT Given the sine and cosine of a conformal
%           latitude, determine an ellipsoidal (geodetic) latitude. A
%           conformal latitude is a latitude in terms of a sphere that is
%           conformal with regard to the ellipsoid. That is, the
%           transformation from the ellipsoid to the sphere preserves
%           angles. This is such that the angle of intersection  between
%           any two lines on the ellipsoid is the same as the corresponding
%           angle on the sphere.
%
%INPUTS: sinChi The real sine of the conformal latitude. This can be a
%               scalar or a matrix of values to convert.
%        cosChi The real cosine of the conformal latitude. This can be a
%               scalar or a matrix of values to convert. This has to be the
%               same size as sinChi.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in Constants.WGS84Flattening
%               is used. This function is only accurate for values from 0
%               to about 0.999. Typically, this will be close to zero, as
%               in Constants.WGS84Flattening.
%
%OUTPUTS: phi The ellipsoidal latitudes in radians. This is the same size
%             as sinChi and cosChi.
%
%Conformal latitudes are discussed in Chapter 3 of [1]. This function
%implements the fixed-point iteration that is described in Section 2.9 of
%[2]. Smaller values of f lead to faster convergence. The maximum number of
%iterations in the function, which is typically never obtained, is set to
%be able to convert unusually large values of f up to about 0.999. A
%different fixed point iteration when given the conformal latitude itself
%is given in [3].
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] Office of Geomatics, "National geospatial-intelligence agency
%    standardization document: Implementation practice: The universal
%    grids and the transverse mercator and polar stereographic
%    map projections," 25 Mar. 2014. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf
%[3] Weisstein, Eric W. "Conformal Latitude." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ConformalLatitude.html
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

N=numel(sinChi);
phi=zeros(size(sinChi));

for curPoint=1:N
    s=sinChi(curPoint);
    %It should converge in 7 or fewer iterations for
    %f=Constants.WGS84Flattening. The upper bound of 5e7 should be enough for
    %f=0.999.
    for curIter=1:3e7
        P=exp(e*atanh(e*s));
        term1=(1+sinChi(curPoint))*P*P;
        term2=(1-sinChi(curPoint));

        sNew=(term1-term2)/(term1+term2);

        if(abs(s-sNew)<=eps(s))
            break
        end
        s=sNew;
    end
    P=exp(e*atanh(e*s));
    c=(1/2)*((1+s)/P+(1-s)*P)*cosChi(curPoint);
    phi(curPoint)=atan2(s,c);
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
