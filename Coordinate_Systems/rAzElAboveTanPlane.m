function rAzEl=rAzElAboveTanPlane(llhTar,llhRef,a,f,tarCoordSys)
%%RAZELABOVETANPLANE Given a reference location and a target location,
%   determine the azimuth in radians East of North in the local tangent
%   plane and the elevation of the target above the local tangent plane.
%
%INPUTS: llhTar A 3XN matrix of [latitude;longitude;ellipsoidal height] for
%               each target with the angles given in radians.
%               Alternatively, if tarCoordSys=1, this can be a set of
%               Cartesian target locations.
%        llhRef The 3X1 [latitude;longitude;ellipsoidal height] for the
%               reference point. Alternatively, a 2X1 vector can be passed,
%               in which case the height will be taken to be 0.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84Flattening is used.
%   tarCoordSys An optional parameter specifying the coordinate system used
%               for the first input of this function. Possible values are:
%               0 (The default if omitted or an empty matrix is passed)
%                 llhTar holds points as [latitude;longitude;height]
%               1 llhTar holds 3D Cartesian locations in Earth-centered
%                 Earth fixed coordinates.
%
%OUTPUTS: rAzEl A 3XN set of [one-way range;azimuth;elevation] with the
%              angles given in radians from the reference point to the
%              target, where the elevation is measured above the local
%              tangent plane to the reference ellipsoid and the azimuth is
%              in radians East of North in the local tangent plane.
%
%The azimuth is with respect to the tangent plane, it is not the azimuth
%that one would travel if one were to follow a geodesic across the Earth
%to get to the target location. Vectors in and orthogonal to the tangent
%plane are obtained using getENUAxes.
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(tarCoordSys))
    tarCoordSys=0;%Given [latitude;longitude;height].
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

zRef=ellips2Cart(llhRef,a,f);
uENU=getENUAxes(llhRef);
if(tarCoordSys==0)
    %Convert to Cartesian.
    zTar=ellips2Cart(llhTar,a,f);
else
    %Cartesian values were directly provided.
    zTar=llhTar;
end

numTar=size(llhTar,2);
rAzEl=zeros(3,numTar);
for k=1:numTar
    vec2Tar=zTar(:,k)-zRef;

    %Get the local offset of the target in local coordinates, where "up" is 
    %one coordinate axis and the other is in the local tangent plane. 
    yTar=uENU(:,3)'*vec2Tar;%How far "up" it goes.
    
    EastComp=uENU(:,1)'*vec2Tar;
    NorthComp=uENU(:,2)'*vec2Tar;

    %How far in the tangent plane it goes.
    xTar=sqrt((EastComp)^2+(NorthComp)^2);
    
    r=norm(vec2Tar);
    thetaAz=atan2(EastComp,NorthComp);%Azimuth, East of North.
    thetaEl=atan(yTar/xTar);%The elevation above the tangent plane.
    rAzEl(:,k)=[r;thetaAz;thetaEl];
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
