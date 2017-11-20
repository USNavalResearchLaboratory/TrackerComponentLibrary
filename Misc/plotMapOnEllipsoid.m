function h=plotMapOnEllipsoid(path2Map,a,f)
%PLOTMAPONELLIPSOID  Plot a map image, where the pixels are given in terms
%                    of an equirectangular projection, on an ellipsoidal
%                    Earth with a black background. This can be useful for
%                    plotting target trajectories on a 3D Earth.
%
%INPUTS: path2Map A string with the path to the mapfile to load. This
%                 should be an image file format that the imread command
%                 can read, such as .png or .jpg. It should also be an
%                 equirectangular projection. That is, each pixel has a
%                 constant width in latitude and longitude. If this
%                 parameter is omitted or an empty matrix is passed, then
%                 the file ./data/NASABlueMarble.jpg is used, where .
%                 refers to the directory of this file.
%               a The semi-major axis of the reference ellipsoid. If
%                 this argument is omitted, the value in
%                 Constants.WGS84SemiMajorAxis is used.
%               f The flattening factor of the reference ellipsoid. If
%                 this argument is omitted, the value in
%                 Constants.WGS84Flattening is used.
%
%OUTPUTS: h The handle to the surfaceplot object created.
%
%This function turns the axes off when plotting and sets the scaling to
%'equal'.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<1||isempty(path2Map))
    ScriptPath=mfilename('fullpath');
    ScriptFolder = fileparts(ScriptPath);
    path2Map=[ScriptFolder,'/data/NASABlueMarble.jpg'];
end

%Read in the map, converting the values to doubles.
theImage = imread(path2Map);
theImage=flipud(theImage);

%For plotting, the image data should be doubles whose values are between 0
%and 1. If theImage is not a bunch of floats or doubles but rather is a
%bunch of integer values, then they should be converted and normalized.
if(isa(theImage,'integer'))
    theImage = double(theImage)/double(intmax(class(theImage)));
else%This is in case it is singles, not doubles.
    theImage = double(theImage);
end

numLat=size(theImage,1);
numLon=size(theImage,2);

latPoints = linspace(-pi/2, pi/2, numLat);
lonPoints = linspace(-pi, pi, numLon);
[lonGrid, latGrid] = meshgrid(lonPoints, latPoints);
clear latPoints lonPoints
numPoints=length(lonGrid(:));

cartPoints=ellips2Cart([latGrid(:)';lonGrid(:)';zeros(1,numPoints)],a,f);
clear latGrid lonGrid
x=reshape(cartPoints(1,:),numLat,numLon);
y=reshape(cartPoints(2,:),numLat,numLon);
z=reshape(cartPoints(3,:),numLat,numLon);

%Plot the image
h=surf(x,y,z,'CData',theImage,'EdgeColor','none','CDataMapping','direct');
%Get rid of the axes (which would also block the black background)
axis('off');
%Make the axes equal to eliminate distortion.
axis('equal');
%Make the background black.
set(gcf, 'Color','k');

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
