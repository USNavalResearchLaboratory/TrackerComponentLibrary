function xyPts=ellips2AzEquidistantProj(latLonPts,latLonRef,a,f)
%%ELLIPS2AZEQUIDISTANTPROJ Given points as [latitude;longitude] (or
%       [latitude;longitude;height])with respect to a reference ellipsoid
%       (or sphere), convert them to 2D (or 3D) points in an azimuthal
%       equidistant projection about a specified point.
%
%INPUTS: latLonPts A 2XN set of [latitude;longitude] points in radians to
%                  convert. Alternately, this can be a 3XN set of
%                  [latitude;longitude;height] points where the height
%                  remains unchanged after the conversion. 
%        latLonRef A 2X1 [latitude;longitude] reference point in raidans
%                  about which the projection is taken.
%                a The semi-major axis of the reference ellipsoid (in
%                  meters). If this argument is omitted or an empty matrix
%                  is passed, the value in Constants.WGS84SemiMajorAxis is
%                  used.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted or an empty matrix is passed,
%                  the value in Constants.WGS84Flattening is used.
%
%OUTPUTS: xyPts A 2XN set of the points converted to an azimuthal
%               equidistant projection about latLonRef or if a height was
%               given in latLonPts, then this is 3XN and the third row is
%               the same as the third row of latLonPts.
%
%The conversion is described in Chapter 25 of [1], where expressions for
%a spherical Earth are given. This is equivalent to finding the distance
%across the reference sphere or ellipsoid and launch angle, and then
%treating those values are polar coordinates and converting them to 2D
%Cartesian. Thus, we implement this function using the
%indirectGreatCircleProb if f=0 and indirectGeodeticProb if f is not zero
%and then we convert the values to 2D Cartesian, noting here that those
%functions return an azimuth East of North, but x is typically East, so the
%conversion given a distance d and an azimuth az is
%point=d*[sin(az);cos(az)];
%
%EXAMPLE:
%This function is demonstrated consistent with its inverse,
%azEquidistantProj2Ellipse, by converting a point from latitude and
%longitude into azimuthal equidistant coordinates and then converting back.
%The relative error of the back-converted point is on the order of finite
%precision limitation.
% latLonPt=[deg2rad([20.631756;-155.551818]),deg2rad([21.295516;-158.002386])];
% latLonRef=deg2rad([20.906029;-157.078918]);
% xyPt=ellips2AzEquidistantProj(latLonPt,latLonRef);
% RelErr=max(max(abs((azEquidistantProj2Ellipse(xyPt,latLonRef)-latLonPt)./latLonPt)))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

hasHeight=(size(latLonPts,1)==3);
numPts=size(latLonPts,2);
if(hasHeight)
    xyPts=zeros(3,numPts);
else
    xyPts=zeros(2,numPts);
end

if(f==0)
    %If using a spherical Earth approximation.   
    for k=1:numPts
        [az,distVal]=indirectGreatCircleProb(latLonRef,latLonPts(1:2,k),a,0);
        xyPts(1:2,k)=distVal*[sin(az);
                            cos(az)];
    end
else
    %If using an ellipsoidal Earth approximation.
    for k=1:numPts
        [az,distVal]=indirectGeodeticProb(latLonRef,latLonPts(1:2,k),a,f);
        xyPts(1:2,k)=distVal*[sin(az);
                            cos(az)];
    end
end

if(hasHeight)
    xyPts(3,:)=latLonPts(3,:);
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
