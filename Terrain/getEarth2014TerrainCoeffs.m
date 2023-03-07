function [C,S]=getEarth2014TerrainCoeffs(M,coeffType)
%%GETEARTH2014TERRAINCOEFFS Get fully normalized spherical harmonic
%                   coefficients for the radial distance offset of the
%                   terrain out from the surface of the GRS80 reference
%                   ellipsoid in meters under various parts of the degree 
%                   2160 Earth2014 model.
%
%INPUTS: M The integer maximum order of the spherical harmonic
%          coefficients obtained. This is a value between 0 and 2160. If
%          this parameter is omitted, the default value is 2160.
% coeffType An optional parameter selecting the type of coefficients to
%          load. There are five variants of this model. Possible values
%          are
%          0 (The default if omitted or an empty matrix is passed) The
%            Earth's surface including water in lakes and major ice
%            sheets. oceans are placed on the GRS80 ellipsoid surface.
%          1 The Earth's bedrock.
%          2 The Earth's toporgraphy, bedrock plus major ice sheets.
%          3 The rock-equivalent topography of the Earth (ice and water
%            masses are condensed to layers of rock)
%          4 Major ice sheets. Everything else is set to the surface of the
%            GRS80 ellipsoid.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the spherical harmonic expansion. If given to a
%           CountingClusterSet class, C(n+1,m+1) is the coefficient of
%           degree n and order m. When a maximum degree of M is used, all C
%           have values for all n from 0 to M and for all m from 0 to n for
%           each n. The coefficients are unitless.
%         S An array holding the coefficient terms that are multiplied by
%           sines in the spherical harmonic expansion. The format of S is
%           the same as that of C.
%
%The model is described in [1]. This function is only for values up to
%degree and order 2160 as the spherical harmonic synthesis routine in the
%Tracker Component Library cannot currently handle the full degree 10,800
%model, because it does not use extended precision arithmetic.
%
%The Earth 2014 data can be downloaded from
%http://ddfe.curtin.edu.au/models/Earth2014/
%To use the data for all of the models, place the files
%Earth2014.SUR2014.degree2160.bshc, Earth2014.BED2014.degree2160.bshc
%Earth2014.TBI2014.degree2160.bshc, Earth2014.RET2014.degree2160.bshc
%and Earth2014.ICE2014.degree2160.bshc in the data folder that is in the
%same folder as this function.
%
%EXAMPLE:
%The fact that this only stores radial distance offsets can be a little
%confusing when one wants to find the location of a point on the Earth's
%surface in the WGS-84 coordinate system/ International Terrestrial
%Reference System (ITRS). Here, we give an example:
% %WGS-84 latitude and longitude in radians.
% latLon=[19.4721;-155.5922]*(pi/180);
% [C,S]=getEarth2014TerrainCoeffs();
% %To use the terrain coefficients, the latitude must be changed from
% %geodetic to geocentric. We will assume that the latitude is given in
% %WGS-84 coordinates, as is standard with GPS.
% latSpher=ellipsLat2SpherLat(latLon(1));
% azElSphere=[latLon(2);latSpher];
% terHeight=spherHarmonicEval(C,S,azElSphere);
% %Add the ellipsoidal radius of the GRS80 reference ellipsoid at the
% %point. Note that this actually changes the ellipsoidal latitude
% %slightly.
% a=Constants.GRS80SemiMajorAxis;
% f=Constants.GRS80Flattening;
% terRad=terHeight+ellipsoidalRadius(latSpher,0,a,f);
% %Convert from spherical to Cartesian coordinates to find the point.
% CartLoc=spher2Cart([terRad;azElSphere])
%
%REFERENCES:
%[1] C. Hirt and M. Rexer, "Earth2014: 1arc-min shape, topography,
%    bedrock and ice-sheet models - available as gridded data and degree-
%    10,800 spherical harmonics," International Journal of Applied Earth
%    Observation and Geoinformation, vol. 39, pp. 103-112, Jul. 2015.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(coeffType))
     coeffType=0;
end

%The Earth2014 terrain coefficient data file, should be located in a data 
%folder that is in the same folder as this file. This finds the path to
%this file.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

switch(coeffType)
    case 0
        fileID=fopen([ScriptFolder,'/data/Earth2014.SUR2014.degree2160.bshc'],'rb');
    case 1
        fileID=fopen([ScriptFolder,'/data/Earth2014.BED2014.degree2160.bshc'],'rb');
    case 2
        fileID=fopen([ScriptFolder,'/data/Earth2014.TBI2014.degree2160.bshc'],'rb');
    case 3
        fileID=fopen([ScriptFolder,'/data/Earth2014.RET2014.degree2160.bshc'],'rb');
    case 4
        fileID=fopen([ScriptFolder,'/data/Earth2014.ICE2014.degree2160.bshc'],'rb');
    otherwise
        error('Unknown coefficient type specified.')
end

data = fread(fileID,Inf,'double');
fclose(fileID); 

%Note that data(1) should be zero. data(1) is the degree of the lowest
%degree coefficient. data(2) is the degree of the maximum degree
%coefficient.
maxDeg = data(2);
maxNumCoeffs=(maxDeg+1)*(maxDeg+2)/2;

if(nargin<1||isempty(M))
    M=maxDeg;
end

totalNumCoeffs=(M+1)*(M+2)/2;

%Extract the data. The additive 2 skips the first two entries, which
%indicate the lowest and highest coefficient degrees.
C=data((2+1):(2+totalNumCoeffs));
S=data((2+maxNumCoeffs+1):(2+maxNumCoeffs+totalNumCoeffs));
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
