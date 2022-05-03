function [latLonEnd,azEnd]=directGeodeticProbGen(latLonStart,azStart,dist,height,useHeightApprox,a,f,numSteps4Circ)
%%DIRECTGEODETICPROBGEN Solve the direct geodetic problem. That is, given
%                     an initial point and an initial bearing an
%                     ellipsoidal Earth, find the end point and final
%                     bearing if one were to travel one or more given
%                     distances along a geodesic curve (the shortest curve
%                     between two points on a curved surface). This
%                     function is general in that it will also solve the
%                     direct geodetic problem at a fixed height above the
%                     reference ellipsoid, if desired.
%
%INPUTS: latLonStart A 2XN matrix of N starting points given in geodetic 
%                   latitude and longitude in radians. that is, each column
%                   is [latitude;longitude]; 
%           azStart An NX1 or 1XN vector of the forward azimuth (initial
%                   headings) at the starting point in radians East of true
%                   North on the reference ellipsoid.
%              dist An NX1 or 1XN vector of the distances in meters 
%                   that will be traveled on the geodesic curve starting at
%                   latLonStart with initial heading azStart where
%                   solutions are desired.
%            height An NX1 or 1XN vector of the heights above the reference
%                   ellipsoid at which the trajectories should be
%                   determined. If this parameter is omitted, then the
%                   default value of 0 is used for all of the trajectories.
%   useHeightApprox If true, and the height is not zero, then a distance
%                   scaling height approximation is used. An equatorial
%                   trajectory traveling a distance dist at a height h
%                   above the reference ellipsoid will only travel
%                   a/(a+h)*dist over the surface. If useHeightApprox is
%                   true, then this scaling applied to any trajectry
%                   (including non-equatorial). If useHeightApprox is
%                   false, a slow iterative optimization used. The default
%                   value if omitted or an empty matrix is passed is
%                   false. This parameter is ignored if height=0.
%                 a The semi-major axis of the reference ellipsoid (in
%                   meters). If this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%                 f The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%     numSteps4Circ If height!=0, and useHeightApprox=false, then an
%                   algorithm propagating a state in ECEF coordinates
%                   around the curved Earth is used. This parameter
%                   determines the number of steps that would be needed
%                   for a target that circumnavigates the globe around the
%                   equator. The default value if this parameter is not
%                   provided is 2000. A value of 6000 appears to be about
%                   the best number for overall precision. Reducing the
%                   number of steps will speed up the function. This
%                   parameter is not used if height=0.
%
%OUTPUTS: latLonEnd A 2XN matrix of geodetic latitude and longitudes of the
%                  final points of the geodesic trajectory given in radians
%                  as [latitude;longitude].
%            azEnd An NX1 vector of the forward azimuth (bearing) at the
%                  ending points in radians East of true North on the
%                  reference ellipsoid.
%
%For zero-altitude, the directGeodeticProb function is called. The solution
%for a non-zero height is based on [1]. The function is significantly
%faster if a zero height is used or if useHeightApprox=true, in which case
%the approximation described in [2] is used.
%
%The algorithm is generally the opposite of indirectgeodeticProbGen.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating Aerial Targets in 3D Accounting for the
%    Earth's Curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, pp. 31-57, Jun. 2015.
%[2] D. F. Crouse, "Worldwide Ground Target State Propagation," in
%    Proceedings of the 23rd International Conference on Information
%    Fusion, 6-9 Jul. 2020.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(numSteps4Circ))
%This is the number of steps that would be taken to travel a distance of
%2*pi*a (circumnavigate the globe the long way). This is a design parameter
%that should be large enough to ensure that the results are accurate.
%However, if it is too large, finite precision errors will accumulate.
    numSteps4Circ=2000;
end

if(nargin<7||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<6||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<5||isempty(useHeightApprox))
    useHeightApprox=false;
end

N=length(dist);

%If no points are provided, return an empty matrix.
if(N==0)
    latLonEnd=[];
    azEnd=[];
    return;
end

latLonEnd=zeros(2,N);
azEnd=zeros(N,1);

if(nargin<4||isempty(height))
	height=zeros(N,1); 
end

if(useHeightApprox)
    %Scale the distance based on the height and then solve the rumb problem
    %with the scaled distance.
    dist(:)=a./(a+height(:)).*dist(:);
    height(:)=0;
end

for curTraj=1:N
    %When on the surface of the reference ellipsoid, then just use the
    %function directGeodeticProb.
    if(height(curTraj)==0)
        [latLonEnd(:,curTraj),azEnd(curTraj)]=directGeodeticProb(latLonStart(:,curTraj),azStart(curTraj),dist(curTraj),a,f);
        continue;
    end

    %If the height is not zero, then use the algorithm for propagating a
    %level, non-maneuvering dynamic model forward in time.

    %This is the number of Runge-Kutta steps on a curved Earth that will be
    %performed.
    numSteps=ceil(numSteps4Circ*(dist(curTraj)/(2*pi*a)));
    stepSize=1/numSteps;

    %We have to calculate the initial location and heading in CARTESIAN
    %ECEF coordinates.
    xyzInit=ellips2Cart([latLonStart(:,curTraj);height(curTraj)],a,f);

    %The initial heading in the LOCAL East-North-Up tangent-plane
    %coordinate system for level flight.
    uh=[sin(azStart(curTraj));cos(azStart(curTraj));0];

    %The speed of the target is such that it travels a distance of dist in
    %1 unit of time.
    speed=dist(curTraj);

    %The initial target state with the position in the GLOBAL ECEF
    %coordinate system and the velocity in a LOCAL flat-Earth coordinate
    %system.
    xInit=[xyzInit;speed*uh];

    %The initial local set of East-North-Up coordinate axes.
    uInit=getENUAxes(latLonStart(:,curTraj));

    %The dynamic model for how the local coordinate system evolves as the
    %target modes.
    uDyn=@(u,x,t)uDotEllipsoid(u,x,a,f);
    aDyn=@(x,t)aPoly(x,3);%A constant velocity model.

    %Compute the Cartesian position and velocity of the target over the
    %entire trajectory. Use a fourth-order Runge-Kutta method.
    [xList,uList]=RungeKCurvedAtTimes(xInit,uInit,[0;1],aDyn,uDyn,stepSize,4);
    plhEnd=Cart2Ellipse(xList(1:3,end),[],a,f);
    %The Cartesian location at the end of the trajectory.
    latLonEnd(:,curTraj)=plhEnd(1:2);

    %The deviation of plhEnd(3) from the desired height can be used as a
    %measure of the error of the estimate.

    %Get the velocity at the end as a unit vector in ECEF Cartesian
    %coordinate.
    vEnd=getGlobalVectors(xList(4:6,end),uList(:,:,end));
    
    %Get a unit vector in the direction of the velocity at the end.
    headingEnd=vEnd/norm(vEnd);
    %Convert the unit vector into an azimuthal heading in radians East of
    %North.
    azEnd(curTraj)=uVec2GeogHeading(plhEnd,headingEnd,a,f);
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
