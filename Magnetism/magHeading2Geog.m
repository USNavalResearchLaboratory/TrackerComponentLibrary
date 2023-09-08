function geoEastOfNorth=magHeading2Geog(points,PEastOfNorth,coeffType,year,a,f)
%MAGHEADING2GEOG Convert a heading in terms of radians clockwise
%                (East) of MAGNETIC North to a heading in terms of radians
%                clockwise from GEOGRAPHIC North as defined on a particular
%                reference ellipsoid. The model for the Earth's magnetic
%                field can be selected. In the rare spots (i.e. near the
%                magnetic poles), where the magnetic field vector points
%                directly towards the reference ellipsoid, the geographic
%                heading will be undefined and a NaN will be returned.
%
%INPUTS: points One or more points given in geodetic latitude and
%               longitude, in radians, and height, in meters where the
%               magnetic headings apply. To convert N headings, points is a
%               3XN matrix with each column having the format
%               [latitude;longitude; height]. All points should be
%               associated with the same time. At the geographic poles, the
%               longitude value determines the orientation of the local
%               coordinate system axes. Thus, geographic headings ARE
%               defined at the poles.
%  PEastOfNorth An NX1 array of N magnetic headings in radians clockwise
%               from North that should be turned into geographic headings.
%     coeffType This specifies the coefficient model for the coefficients.
%               If one wishes to explcitely pass a model, then this is a
%               structure with members C, S, aH, and cH, which are defined
%               in the same manner as the return values from the function
%               getWMMCoeffs. Otherwise, a string can be passed. Possible
%               values are:
%               'WMM' (The default if omitted or an empty matrix is
%                     passed). Use the World Magnetic Model via the
%                     function getWMMCoeffs.
%               'IGRF' Use the International Geomagnetic Reference Field
%                     model via the function getIGRFCoeffs.
%          year This is only used if coeffType is a string. A decimal
%               number indicating a year in the Gregorian calendar as
%               specified by UTC. For example, halfway through the year
%               2012 would be represented as 2012.5. The precision of the
%               model is not sufficiently fine that leap seconds matter. If
%               this parameter is omitted, then the last reference epoch of
%               the geomagnetic model is used.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84Flattening is used.
%
%OUTPUT: geoEastOfNorth The headings converted to radians clockwise of
%                       geographic North on the reference ellipse.
%
%First, spherical harmonic coefficients for the geomagnetic model are
%obtained at the desired time and points. Then, the magnetic flex vector B
%is determined at each of the points. Strictly speaking, the direction of
%magnetic North would be the component of the magnetic flux vector
%projected onto the local gravitationally horizontal plane at each point.
%However, since B is already rather imprecise, nothing is really lost by
%using a projection onto the local tangent plane of the reference
%ellipsoid of the Earth. Thus, rotating B about the negative of the local
%vertical to the reference ellipsoid by PEastOfNorth provides a vector
%whose projection into the local tangent plane points in the desired
%geographic heading.
%
%This function makes use of the functions getWMM2010Coeffs and 
%spherHarmonicEval to determine the magnetic flux vector at the specified
%points as well as the function getENUAxes to determine the local
%East-North-Up coordinate axes.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<3||isempty(coeffType))
    coeffType='WMM';
end

if(nargin<4)
    year=[];
end

%Get the fully normalized spherical harmonic coefficients for the
%geomagnetic model.
if(isstruct(coeffType))
    %If the user explicitly provided the coefficients.
    C=coeffType.C;
    S=coeffType.S;
    aH=coeffType.aH;
    cH=coeffType.cH;
else
    switch(coeffType)
        case 'WMM'
            [C,S,aH,cH]=getWMMCoeffs(year);
        case 'IGRF'
            [C,S,aH,cH]=getIGRFCoeffs(year);
        otherwise
            error('Unknown coefficients selected')
    end
end

[~,gradV]=spherHarmonicEval(C,S,ellips2Sphere(points,a,f),aH,cH);
B=-gradV;
%B is now a matrix of vectors of the magnetic flux of the Earth's field at
%the points. The direction of B in the local tangent plane to the reference
%ellipsoid approximately defined magnetic north.

numPoints=size(points,2);
geoEastOfNorth=zeros(numPoints,1);%Allocate space
for curPoint=1:numPoints
    %-B(:,curPoint) defines the direction of magnetic North at the current
    %point. We want the direction vector PEastOfNorth(curPoint). To go East
    %of North, one must rotate that number of degrees about the local
    %"Down" vector.
    u=getENUAxes(points(:,curPoint),false,a,f);
    
    %The negative of u(:,3) is a unit vector in the down direction.
    vRot=rotateVectorAxisAng(B(:,curPoint),-u(:,3),PEastOfNorth(curPoint));
    
    %Given the rotated field vector, only the components in the local
    %tangent plane matter.
    vEast=dot(vRot,u(:,1));
    vNorth=dot(vRot,u(:,2));
    %Find the angle East of North.
    geoEastOfNorth(curPoint)=atan2(vEast,vNorth);
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
