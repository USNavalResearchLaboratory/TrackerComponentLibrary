function [ellipsPoints,coeffData]=MSL2EllipseHelmert(points,useNGAApprox,coeffData)
%%MSL2ELLIPSEHELMERT Given a point whose latitude and longitude are in 
%                    WGS-84 ellipsoidal coordinates, but whose height is a
%                    height above mean-sea level (MSL) (assuming a Helmert
%                    projection), convert the height components to WGS84
%                    ellipsoidal heights using Helmert's projection method
%                    with the EGM2008 gravitational (and terrain
%                    correction) models. This differs from the more precise
%                    but seldom-used Pizzetti's projection method in that
%                    the curvature of the plumb line is completely ignored.
%
%INPUTS: points     One or more points given in WGS-84 geodetic latitude
%                   and longitude, in radians, and a height above mean-seal
%                   level (in meters) for the specified latitude and
%                   longitude. To convert N points, points is a 3XN matrix
%                   with each column having the format
%                   [latitude;longitude; MSLheight].
%     useNGAApprox  If one wishes the intermediate tide-free geoid
%                   computation to match the results of the National
%                   Geospatial Intelligence Agency's code to three digits
%                   after the decimal point, then this should be true.
%                   Setting this to false might produce results that
%                   are marginally more accurate. The default is false if
%                   this parameter is omitted.
%       coeffData  A set of pre-loaded coefficients that can speed up the
%                  computation by eliminating the need to compute them on
%                  the fly. coeffData is a structures with elements C, S,
%                  CZeta and SZeta where C and S have been set by the
%                  output of getEGM2008WGS84TCoeffs and CZeta and SZeta
%                  have been set by the output of getEGM2008Zeta2NCoeffs.
%                  If this parameter is omitted, the coefficients are
%                  computed during the execution of the function.
%
%OUTPUTS ellipsPoints A 3XN array of the converted points, where each vector
%                   contains [latitude;longitude;height]; The latitude
%                   and longitude are unchanged from the input data, but
%                   the heights have been converted to from MSL heights to
%                   WGS-84 ellipsoidal heights.
%         coeffData The coeffData coefficients that can be passed to
%                   another call of ellips2MSLHelmert or MSL2EllipseHelmert
%                   to make it faster.
%
%As described in Chapter 5.5 of [1], the MSL height using Helmert's
%projection is just the ellipsoidal height minus the geoid height. This
%function just calls getEGM2008GeoidHeight to get the geoidal height and
%subtracts it off.
%
%Helmert's projection calculated the height using a line from the point to
%the reference ellipsoid, subtracting off the distance from the point to
%the geoid along the line. Pizzetti's projection follows the curved plumb
%line down to the geoid. Thus, the latitude and longitude of the projected
%point on the geoid are not the same as that of the point being projected.
%The MSL height obtained is also slightly different. However, the
%difference is small and Pizzetti's projection is difficult to use, so
%Helmert's projection is almost always used.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    useNGAApprox=false;
end

if(nargin>2)
    [geoidHeight,coeffData]=getEGM2008GeoidHeight(points(1:2,:),1,useNGAApprox,coeffData);
else
    [geoidHeight,coeffData]=getEGM2008GeoidHeight(points(1:2,:),1,useNGAApprox);
end

    ellipsHeight=points(3,:)+geoidHeight';
    numPoints=size(points,2);
    
    ellipsPoints=zeros(3,numPoints);
    ellipsPoints(1:2,:)=points(1:2,:);
    ellipsPoints(3,:)=ellipsHeight;

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
