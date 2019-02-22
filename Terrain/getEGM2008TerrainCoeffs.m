function [C,S]=getEGM2008TerrainCoeffs(M)
%%GETEGM2008TERRAINCOEFFS Get fully normalized spherical harmonic
%                   coefficients for the terrain height (expressed in terms
%                   of meters above mean seal level) that is part of the 
%                   National Geospatial Intelligence Agency's (NGA's) Earth
%                   Gravitation Model 2008 (EGM2008). The MSL terrain
%                   height at a particular latitude and longitude is
%                   necessary for the conversion between Cartesian
%                   elevation and an MSL altitude. The coefficients are
%                   derived from the Digital Topographic Model 2006.0.
%
%INPUTS: M The integer maximum order of the spherical harmonic coefficients
%          obtained. This is a value between 2 and 2190. If this parameter
%          is omitted, the default value is 2190.
%
%OUTPUTS: C An array holding the coefficient terms that are multiplied by
%           cosines in the spherical harmonic expansion. If given to a
%           CountingClusterSet class, then C(n+1,m+1) is the coefficient of
%           degree n and order m. When a maximum degree of M is used, all C
%           have values for all n from 0 to M and for all m from 0 to n for
%           each n. The coefficients are unitless.
%         S An array holding the coefficient terms that are multiplied by
%           sines in the harmonic expansion. The format of S is the same as
%           that of C.
%
%The EGM2008 model is documented in [1] and [2].
%
%More on using the spherical harmonic coefficients is given in the comments
%for the function spherHarmonicEval and the format and use of
%the coefficients is also documented in 
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/README_FIRST.pdf
%The coefficients are for finding heights above mean sea level, which means
%that the geoid height must be added to get the true height above the
%WGS-84 reference ellipsoid. See the function getEGMGeoidHeight to obtain
%the height. The coefficients can be downloaded from
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html
%in a file called Coeff_Height_and_Depth_to2190_DTM2006.0 .
%
%As reading the data from the file can be slow, this function first checks
%for a .mat file with the coefficients in it in the folder named "data"
%that is in the same folder as this function. If
%such a file does not exist, then the coefficients are read in from the
%appropriate .dat file and if all 2190 coefficients are requested, then a
%.mat file is created so that the text in the .dat file need not be read
%again, because reading the text file is extremely slow.
%
%EXAMPLE:
%The fact that this stores elevations above the geoid can make it confusing
%to find a point on the terrain at a particular latitude and longitude.
%Here we give an example:
% %WGS-84 latitude and longitude in radians.
% latLon=[19.4721;-155.5922]*(pi/180);
% [C,S]=getEGM2008TerrainCoeffs();
% %To use the terrain coefficients, the latitude must be changed from
% %geodetic to geocentric. We will assume that the latitude is given in
% %WGS-84 coordinates, as is standard with GPS.
% latSpher=ellipsLat2SpherLat(latLon(1));
% azElSphere=[latLon(2);latSpher];
% orthoHeight=spherHarmonicEval(C,S,azElSphere);
% %The tide-free geoid should be the reference for the model. The function
% %getEGMGeoidHeight loads the coefficients for the gravitational model,
% %though we can also pass then in a structure to speed it up during
% %repeated calls.
% geoidHeight=getEGMGeoidHeight(latLon,0,true,0);
% CartLoc=ellips2Cart([latLon;geoidHeight+orthoHeight])
%
%REFERENCES:
%[1] N. K. Pavlis, S. A. Holmes, S. C. Kenyon, and F. J. K., "The
%    development and evaluation of the Earth gravitational model 2008
%    (EGM2008)," Journal of Geophysical Research, vol. 117, no. B4, Apr.
%    2012.
%[2] N. K. Pavlis, S. A. Holmes, S. C. Kenyon, and F. J. K., "Correction to
%    "The development and evaluation of the Earth gravitational model 2008
%    (EGM2008)"," Journal of Geophysical Research, vol. 118, 2013.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(M))
    M=2190;
end

if(M>2190)
    error('The EGM2008 terrain model only has coefficients up to degree 2190') 
end

totalNumCoeffs=(M+1)*(M+2)/2;

%The EGM20008 terrain coefficient data file, should be located in a data 
%folder that is in the same folder as this file. This find the path to this
%file.
ScriptPath=mfilename('fullpath');
ScriptFolder = fileparts(ScriptPath);

%First, see if a .mat file with all of the data exists. If so, then use
%that and ignore everything else.
if(exist([ScriptFolder,'/data/Coeff_Height_and_Depth_to2190_DTM2006.0.mat'],'file'))
    load([ScriptFolder,'/data/Coeff_Height_and_Depth_to2190_DTM2006.0.mat'],'C','S');

    C=C(1:totalNumCoeffs);
    S=S(1:totalNumCoeffs);
    return 
end

%If the .mat file does not exist, then assume that the coefficients must be
%read from the text file provided by the NGA.

%Allocate space for the coefficients.
emptyData=zeros(totalNumCoeffs,1);
C=CountingClusterSet(emptyData);
S=CountingClusterSet(emptyData);

%Read in the data up to the specified order.
fileID=fopen([ScriptFolder,'/data/Coeff_Height_and_Depth_to2190_DTM2006.0']);
data=textscan(fileID,'%d %d %f %f %f %f',totalNumCoeffs);
fclose(fileID);

%Put the data into the ClusterSet class instances.
numRows=totalNumCoeffs;
for curRow=1:numRows
    n=data{1}(curRow);%The degree of the coefficient.
    m=data{2}(curRow);%The order of the coefficient.
    C(n+1,m+1)=data{3}(curRow);
    S(n+1,m+1)=data{4}(curRow);
end

C=C.clusterEls;
S=S.clusterEls;

%Save the coefficients into a .mat file so that they can be read more
%quickly next time.
if(M==2190)
    save([ScriptFolder,'/data/Coeff_Height_and_Depth_to2190_DTM2006.0.mat'],'C','S')
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
