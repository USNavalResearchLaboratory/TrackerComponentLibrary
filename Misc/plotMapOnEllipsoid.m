function h=plotMapOnEllipsoid(param1,param2,param3,param4)
%PLOTMAPONELLIPSOID Plot a map image, where the pixels are given in terms
%                   of an equirectangular projection, on an ellipsoidal
%                   Earth. This can be useful for plotting target
%                   trajectories on a 3D Earth.
%
%INPUTS: The function can be parameterized either as
%        plotMapOnEllipsoid(path2Map,a,f)
%        or
%        plotMapOnEllipsoid(hAxis,path2Map,a,f)
%        These inputs are
%        hAxis A matlab.graphics.axis.Axes or matlab.ui.control.UIAxes
%              handle for the axes on which the map should be plotted.
%     path2Map A string with the path to the mapfile to load. This should
%              be an image file format that the imread command can read,
%              such as .png or .jpg. It should also be an equirectangular
%              projection. That is, each pixel has a constant width in
%              latitude and longitude. If this parameter is omitted or an
%              empty matrix is passed, then the file
%              ./data/NASABlueMarble.jpg is used, where . refers to the
%              directory of this file. If the string "" is passed, then no
%              map will be loaded. Rather, a set of 50 latitude and
%              longitude lines will be drawn over a gray ellipsoid.
%            a The semi-major axis of the reference ellipsoid. If this
%              argument or an empty matrix is passed, the value in
%              Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%
%OUTPUTS: h The handle to the surfaceplot object created.
%
%This function turns the axes off when plotting and sets the scaling to
%'equal'. One might want to call set(gcf,'Color','k') after plotting to
%make the background black.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin==0)
    hAxes=[];
    path2Map=[];
    a=[];
    f=[];
else
    if(isa(param1,'matlab.graphics.axis.Axes')||isa(param1,'matlab.ui.control.UIAxes'))
        %If the axes on which to plot were passed.
        hAxes=param1;
        if(nargin>1)
            path2Map=param2;
        else
            path2Map=[];
        end
        
        if(nargin>2)
            a=param3;
        else
            a=[];
        end
        
        if(nargin>3)
            f=param4;
        else
            f=[];
        end
    else
        hAxes=[];
        %If the axes on which to plot were not passed.
        if(nargin==4)
           error('If 4 inputs are provided, the first must be of type matlab.graphics.axis.Axes.')
        end

        path2Map=param1;
        
        if(nargin>1)
            a=param2;
        else
            a=[];
        end
        
        if(nargin>2)
            f=param3;
        else
            f=[];
        end
    end
end

if(isempty(f))
    f=Constants.WGS84Flattening;
end

if(isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(isempty(path2Map))
    ScriptPath=mfilename('fullpath');
    ScriptFolder=fileparts(ScriptPath);
    path2Map=[ScriptFolder,'/data/NASABlueMarble.jpg'];
end

if(strcmp(path2Map,""))
    numLat=50;
    numLon=50;
    theImage=[];
else
    %Read in the map, converting the values to doubles.
    theImage=imread(path2Map);
    theImage=flipud(theImage);

    %For plotting, the image data should be doubles whose values are between 0
    %and 1. If theImage is not a bunch of floats or doubles but rather is a
    %bunch of integer values, then they should be converted and normalized.
    if(isa(theImage,'integer'))
        theImage=double(theImage)/double(intmax(class(theImage)));
    else%This is in case it is singles, not doubles.
        theImage=double(theImage);
    end

    numLat=size(theImage,1);
    numLon=size(theImage,2);
end

latPoints=linspace(-pi/2, pi/2, numLat);
lonPoints=linspace(-pi, pi, numLon);
[lonGrid,latGrid]=meshgrid(lonPoints,latPoints);
clear latPoints lonPoints
numPoints=length(lonGrid(:));

cartPoints=ellips2Cart([latGrid(:)';lonGrid(:)';zeros(1,numPoints)],a,f);
clear latGrid lonGrid
x=reshape(cartPoints(1,:),numLat,numLon);
y=reshape(cartPoints(2,:),numLat,numLon);
z=reshape(cartPoints(3,:),numLat,numLon);

%Plot the image
if(~isempty(theImage))
    if(~isempty(hAxes))
        h=surf(hAxes,x,y,z,'CData',theImage,'EdgeColor','none','CDataMapping','direct');
    else
        h=surf(x,y,z,'CData',theImage,'EdgeColor','none','CDataMapping','direct');
        hAxes=gca;
    end
else
    if(~isempty(hAxes))
        h=surf(hAxes,x,y,z,'EdgeColor','black');
    else
        h=surf(x,y,z,'EdgeColor','black');
        hAxes=gca;
    end
    colormap(hAxes,gray);
end
%Get rid of the axes (which would also block the black background)
axis(hAxes,'off');
%Make the axes equal to eliminate distortion.
axis(hAxes,'equal');

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
