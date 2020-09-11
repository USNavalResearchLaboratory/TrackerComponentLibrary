function conformalLat=ellipsLat2ConformLat(ellipsLat,f,algorithm)
%%ELLIPSLAT2CONFORMLAT Convert an ellipsoidal (geodetic) latitude of a
%                    point on the surface of a reference ellipsoid into a
%                    latitude in terms of a sphere that is conformal with
%                    regard to the ellipsoid. That is, the transformation
%                    from the ellipsoid to the sphere preserves angles.
%                    This is such that the angle of intersection  between
%                    any two lines on the ellipsoid is the same as the
%                    corresponding angle on the sphere.
%
%INPUTS: ellipsLat A vector or matrix of ellipsoidal (geodetic) latitudes
%                  in radians. These range from -pi/2 to pi/2.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted, the value in
%                  Constants.WGS84Flattening is used.
%        algorithm An optional parameter specifying how the conformal
%                  latitude shall be determined. Both methods are
%                  algebraically equivalent and differ only in complexity
%                  and finite precision errors. Possible values are:
%                  0 (The default if omitted or an empty matrix is passed)
%                    Use the formula given in Section 2.4 of [2].
%                  1 Use the formula in Equation 3.1a of [1]. 
%
%OUTPUTS: conformalLat The latitudes in ellipsLat converted to conformal
%                  latitudes. This ranges from -pi/2 to pi/2.
%
%The conformal latitude is Defined in Chapter 3 of [1]. It can play a role
%in certain types of transverse mercator projections. It is also given in
%Section 2.4 of [2], though with less of an explanation.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] Office of Geomatics, "National geospatial-intelligence agency
%    standardization document: Implementation practice: The universal
%    grids and the transverse mercator and polar stereographic
%    map projections," 25 Mar. 2014. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

switch(algorithm)
    case 0
        isoLat=ellipsLat2IsoLat(ellipsLat,f);
        conformalLat=asin(tanh(isoLat));
    case 1
        sinLat=sin(ellipsLat);
        e=sqrt(2*f-f^2);

        conformalLat=2*atan(sqrt(((1+sinLat)./(1-sinLat)).*((1-e*sinLat)./(1+e*sinLat)).^e))-pi/2;
    otherwise
        error('Unknown algorithm specified.')
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
