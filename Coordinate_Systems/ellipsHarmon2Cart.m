function cartPoints=ellipsHarmon2Cart(pointsHarmon,E)
%%CART2ELLIPSHARMON Convert from ellipsoidal harmonic coordinates to
%                   Cartesian coordinates. Note that ellipsoidal harmonic
%                   coordinates are not the same as ellipsoidal
%                   coordinates.
%
%INPUTS: pointsHarmon One or more points given in terms of ellipsoidal
%                   harmonic reduced latitude and longitude in radians and
%                   a semi-major axis in meters that are to be converted to
%                   Cartesian coordinates. To convert N points,
%                   pointsHarmon is a 3XN matrix with each column
%                   having the format
%                   [reduced latitude;longitude; semi-major axis].
%                 E The linear eccentricity defining the ellipsoidal
%                   harmonic coordinate system. If this parameter is
%                   omitted, then the linear eccentricity of the WGS84
%                   reference ellipsoid is used.
%
%OUTPUTS: cartPoints For N points, cartPoints is a 3XN matrix of the
%                    converted points with each column having the format
%                    [x;y;z].
%
%Equation 1-150 in Chapter 1.15 of [1] is used for the conversion, where
%the complement of the reduced latitude has been replaced by the reduced
%latitude.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(E))
        a=Constants.WGS84SemiMajorAxis;
        f=Constants.WGS84Flattening;
        b=a*(1 - f);
        E=sqrt(a^2-b^2);
    end
    
    beta=pointsHarmon(1,:);
    lambda=pointsHarmon(2,:);
    u=pointsHarmon(3,:);

    x=sqrt(u.^2+E^2).*cos(beta).*cos(lambda);
    y=sqrt(u.^2+E^2).*cos(beta).*sin(lambda);
    z=u.*sin(beta);

    cartPoints=[x;y;z];
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
