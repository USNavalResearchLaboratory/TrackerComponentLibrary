function PEastOfNorth=geogHeading2Mag(points,geoEastOfNorth,coeffType,year,a,f)
%GEOGHEADING2MAG Convert a heading in terms of radians clockwise (East)
%                from GEOGRAPHIC North to a heading in terms of radians
%                clockwise of MAGNETIC North as defined on a particular
%                reference ellipsoid. The model for the Earth's magnetic
%                field can be selected. In the rare spots (i.e. near the
%                magnetic poles), where the magnetic field vector points
%                directly towards the reference ellipsoid, the geographic
%                heading will be undefined and a NaN will be returned.
%
%INPUTS: points One or more points given in geodetic latitude and
%               longitude, in radians, and height, in meters where the
%               geographic headings apply. To convert N headings, points is
%               a 3XN matrix with each column having the format
%               [latitude;longitude; height]. All points should be
%               associated with the same time. At the geographic poles, the
%               longitude value determines the orientation of the local
%               coordinate system axes. Thus, geographic headings ARE
%               defined at the poles.
% geoEastOfNorth An NX1 array of N geographic headings in radians
%               clockwise from North that should be turned into magnetic
%               headings.
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
%          year A decimal number indicating a year in the Gregorian
%               calendar as specified by UTC. For example, halfway through
%               the year 2012 would be represented as 2012.5. The
%               precision of the model is not sufficiently fine that leap
%               seconds matter. If this parameter is omitted, then the
%               last reference epoch of the geomagnetic model is used.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84Flattening is used.
%
%OUTPUT: GeoEastOfNorth The headings converted to radians clockwise of
%                       magnetic North on the reference ellipse.
%
%The magnetic flux vector points to magnetic North. One finds the angle
%between geographic North and the projection of the negative magnetic flux
%vector in the local tangent plane. That angle subtracted from the
%geographic heading gives the magnetic heading.
%
%The consistency of the magnetic heading conversion routines can be
%validated using
% 
% pointEllips=[40*pi/180;170*pi/180;30];
% geoEastOfNorth=64*pi/180;
% year=2014.0;
% geoEastOfNorthConv=geogHeading2Mag(pointEllips,magHeading2Geog(pointEllips,geoEastOfNorth,'IGRF',year),'IGRF',year);
% abs(geoEastOfNorthConv-geoEastOfNorth)
%
%If everything is consistent, then abs(geoEastOfNorthConv-geoEastOfNorth)
%should be zero.
%
%This function makes use of the functions getWMMCoeffs and 
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
%ellipsoid approximately defined magnetic North.

numPoints=size(points,2);
PEastOfNorth=zeros(numPoints,1);%Allocate space
for curPoint=1:numPoints
    %-B(:,curPoint) defines the direction of magnetic North at the current
    %point. We want the direction vector PEastOfNorth(curPoint). To go East
    %of North, one must rotate that number of degrees about the local
    %"Down" vector.
    u=getENUAxes(points(:,curPoint),false,a,f);
        
    %Get the components of the magnetic flux vector (-B points North in the
    %local tangent plane).
    vEast=dot(B(:,curPoint),u(:,1));
    vNorth=dot(B(:,curPoint),u(:,2));
    
    %The angle in radians that the B vector makes from the North axis
    %(going clockwise) is 
    BAng=atan2(vEast,vNorth);
    
    %The angle offset (in radians) between the current point and the North
    %vector.
    PEastOfNorth(curPoint)=geoEastOfNorth(curPoint)-BAng;
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
