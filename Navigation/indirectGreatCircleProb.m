function [azStart,dist,azEnd,latLonWaypoints]=indirectGreatCircleProb(latLonStart,latLonEnd,r,N,algorithm,waypointType,noJumps)
%%INDIRECTGREATCIRCLEPROB Solve the indirect great circle problem. That is,
%                      given two points on a spherical Earth, find the
%                      initial bearing and distance one must travel to
%                      take the shortest (geodesic) path between the
%                      points. Additionaly, this function can compute
%                      waypoints on the path. This type of navigation is
%                      sometimes referred to as "great circle sailing."
%                      When used for navigation, it is less accurate than
%                      an ellipsoidal Earth approximation.
%
%INPUTS: latLonStart The 2XnumPts set of numPts initial points given in
%               latitude and longitude in radians in the format
%               [latitude;longitude] (on a reference sphere, latitude is
%               spherical elevation; longitude is spherical azimuth). If
%               all these points are the same but latLonEnd varies, then a
%               single 2X1 vector can be passed. Extra rows, if passed,
%               will be ignored.
%     latLonEnd The 2XnumPts final points given in latitude and longitude
%               in adians in the format [latitude;longitude]. If all these
%               points are the same but latLonStart varies, then a single
%               2X1 vector can be passed.
%             r The assumed radius of the spherical Earth model. If omitted
%               or an empty matrix is passed, the default of
%               r=osculatingSpher4LatLon(latLonStart) is used.
%             N The number of waypoints, besides the initial and final
%               points on the trajectory, to produce. If omitted or an
%               empty matrix is passed, then no waypoints are generated and
%               latLonWaypoints is an empty matrix.
%     algorithm An optional parameter selecting the algorithm to use.
%               Possible values are:
%               0 Use the algorithm of [2], which avoids issues with
%                 singularities at the poles.
%               1 (The default if omitted or an empty matrix is passed and
%                 waypointType=1) Use the COFI algorithm of [1]. Latitudes
%                 within 2^24*eps(pi/2) of +/-pi/2 (the North and South
%                 poles) will be clipped to that bound. This avoids a
%                 singularity at the poles.
%               2 (The default if omitted or an empty matrix is passed and
%                 waypointType=0) Obtain angles using Equation 5-4b of [3],
%                 which is stable at the poles. Obtaincthe distance using a
%                 cross product formula, which is also used in algorithm 1,
%                 as discussed below.
%    waypointType An optional parameter indicating how the waypoints will
%                 be spaced. Possible values are:
%                 0 (The default if omitted or an empty matrix is passed)
%                   Space the N waypoints between the start and end
%                   uniformly in distance.
%                 1 Space the N waypoints between the start and end
%                   uniformly in longitude. This option is only available
%                   for algorithm 0. If a trajectory is meridional, then
%                   only the starting and stopping points will be returned.
%
%OUTPUTS: azStart The 1XnumPts scalar forward azimuths at the starting
%                 point in radians East of true North on the reference
%                 sphere.
%            dist The 1XnumPts distances on the sphere-approximated Earth
%                 between the starting and stopping points.
%           azEnd The 1XnumPts forward azimuth at the ending point in
%                 radians East of true North on the reference sphere.
% latLonWaypoints A 2X(N+2) set of waypoints on the trajectory.
%                 latLonWaypoints(:,i) is the ith points as
%                 [latitude;longitude] in radians. The first point is
%                 latLonStart and the final points is latLonEnd, though the
%                 points are made to avoid jumps of 2*pi in longitude,
%                 which can mean that the last point in this can be off
%                 from the provided last point by a factor of 2*pi. The
%                 other points are equally spaced along the trajectory. If
%                 the trajectory is meridional and waypointType=1, then
%                 only the starting and ending points will be returned.
%
%Great circles are geodesics on a sphere. A geodesic path between two
%points on a sphere is much simpler to determine as compared to an
%ellipsoid (which is what one gets using the indirectGeodeticProb
%function). For the spherical model, the latitude is taken as the elevation
%and the longitude as the azimuth on the sphere (we assume coordinate
%system 0 in the spher2Cart function). The latitudes and longitudes
%returned as waypoints by this function can be used to piece together an
%_approximate_ geodesic trajectory on the reference ellipsoid. One can
%obtain headings and distances between the reference points using rhumb
%lines.
%
%The azimuthal values returned when one of the points is at the pole
%depends on the longitude value given for the points. The getENUAxes
%function can return East-North-Up axes anywhere on the globe (for a
%sphere, the ellipsoidal flattening is 0). When given a pole, the result is
%the equivalent one would get by taking a point with the same latitude and
%with the latitude magnitude a tiny epsilon less than pi/2.
%
%The implementation of the COFI algorithm of [1] has been modified so that
%the angular difference D is obtained using the angBetweenVecs given F and
%T rather than the technique derived from the dot product relation starting
%in Equation 13 in [1]. This is because the cross product relation used in
%angBetweenVecs is more accurate. This is also used for the distance in
%algorithm 2. Minor changes had to be made to algorithm 1 to deal with
%trajectories going West, since the derivation in [1] is for East-bound
%trajectories.
%
%EXAMPLE 1:
%This is example 2 from [1]. It is a trajectory that crosses the
%international date and and goes from the Northern hemisphere to the
%southern hemisphere. We also compute the reverse path and show that the
%start and end azimuth angles produced in each direction are consistent
%with each other (when the same Earth radius is used eahc time). We then
%plot the trajectory on an image of the spherical Earth. For better
%plotting, the radius of the Earth has been normalized to 1.
% N=500;
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% 
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% [azStartFwd,distFwd,azEndFwd,latLonWayPoints]=indirectGreatCircleProb(latLonStart,latLonEnd,1,N);
% [azStartRev,distRev,azEndRev,latLonWayPointsRev]=indirectGreatCircleProb(latLonEnd,latLonStart,1,N);
% 
% %All four of these should be approximately zero if the algorithm works
% %forwards and backwards:
% distFwd-distRev
% wrapRange(azStartFwd-wrapRange(azEndRev+pi,-pi,pi),-pi,pi)
% wrapRange(azStartRev-wrapRange(azEndFwd+pi,-pi,pi),-pi,pi)
% max(max(abs(wrapRange(latLonWayPoints-fliplr(latLonWayPointsRev),-pi,pi))))
% 
% xStartCart=ellips2Cart([latLonStart;0],1,0);
% xEndCart=ellips2Cart([latLonEnd;0],1,0);
% %Give the points a nonzero elevation so the lines are easier to see.
% pathPoints=ellips2Cart([latLonWayPoints;0.02*ones(1,N+2)],1,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],1,0);
% scatter3(xStartCart(1),xStartCart(2),xStartCart(3),100,'filled')
% scatter3(xEndCart(1),xEndCart(2),xEndCart(3),100,'filled')
% plot3(pathPoints(1,:),pathPoints(2,:),pathPoints(3,:),'-r','linewidth',4)
% view([-75.368365161957513,3.494544111878512])
%
%EXAMPLE 2:
%This is the same as example 1, except one point is given at the North
%pole.
% N=500;
% latStart=pi/2;%North pole.
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% 
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% [azStartFwd,distFwd,azEndFwd,latLonWayPoints]=indirectGreatCircleProb(latLonStart,latLonEnd,1,N);
% [azStartRev,distRev,azEndRev,latLonWayPointsRev]=indirectGreatCircleProb(latLonEnd,latLonStart,1,N);
% 
% %All four of these should be approximately zero if the algorithm works
% %forwards and backwards:
% distFwd-distRev
% wrapRange(azStartFwd-wrapRange(azEndRev+pi,-pi,pi),-pi,pi)
% wrapRange(azStartRev-wrapRange(azEndFwd+pi,-pi,pi),-pi,pi)
% max(max(abs(wrapRange(latLonWayPoints-fliplr(latLonWayPointsRev),-pi,pi))))
% 
% xStartCart=ellips2Cart([latLonStart;0],1,0);
% xEndCart=ellips2Cart([latLonEnd;0],1,0);
% pathPoints=ellips2Cart([latLonWayPoints;0.02*ones(1,N+2)],1,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],1,0);
% scatter3(xStartCart(1),xStartCart(2),xStartCart(3),100,'filled')
% scatter3(xEndCart(1),xEndCart(2),xEndCart(3),100,'filled')
% plot3(pathPoints(1,:),pathPoints(2,:),pathPoints(3,:),'-r','linewidth',4)
%
%REFERENCES:
%[1] C.-L. Chen, P.-F. Liu, and W.-T. Gong, "A simple approach to great
%    circle sailing: The COFI method," The Journal of Navigation, vol. 67,
%    no. 3, pp. 403-418, May 2014.
%[2] D. F. Crouse, "Singularity-free great-circle sailing," Naval Research
%    Laboratory, Washington, DC, Tech. Rep. NRL/5340/MR-2021/4, 26
%    Jul. 2021.
%[3] J. P. Snyder, "Map projections-a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(noJumps))
    noJumps=true;
end

if(nargin<6||isempty(waypointType))
    waypointType=0;
end

if(nargin<5||isempty(algorithm))
    if(waypointType==1)
        algorithm=1;
    else
        algorithm=2;
    end
end

if(nargin<4)
    N=[];%Generate no waypoints.
end

if(nargin<3||isempty(r))
    r=osculatingSpher4LatLon(latLonStart);
end

if(nargout<3)
    N=[];%Do not generate waypoints.
end

numPts1=size(latLonStart,2);
numPts2=size(latLonEnd,2);
numPts=max(numPts1,numPts2);

if(numPts1>numPts2)
    latLonEnd=repmat(latLonEnd,[1,numPts]);
elseif(numPts2>numPts1)
    latLonStart=repmat(latLonStart,[1,numPts]);
end

azStart=zeros(1,numPts);
dist=zeros(1,numPts);
azEnd=zeros(1,numPts);
if(~isempty(N))
    latLonWaypoints=zeros(2,N+2,numPts);
else
    latLonWaypoints=[];
end

switch(algorithm)
    case 0
        if(waypointType~=0)
            error('This algorithm only supports waypointType=0.');
        end
        
        for k=1:numPts
            [azStart(k),dist(k),azEnd(k),latLonWaypointsCur]=indirectGreatCircleProbCrouse(latLonStart(1:2,k),latLonEnd(1:2,k),r,N);
            if(~isempty(N))
                latLonWaypoints(:,:,k)=latLonWaypointsCur;
            end
        end
    case 1
        for k=1:numPts
            [azStart(k),dist(k),azEnd(k),latLonWaypointsCur]=indirectGreatCircleProbChen(latLonStart(1:2,k),latLonEnd(1:2,k),r,N,waypointType);
            if(~isempty(N))
                latLonWaypoints(:,:,k)=latLonWaypointsCur;
            end
        end
    case 2
        for k=1:numPts
            [azStart(k),dist(k),azEnd(k),latLonWaypointsCur]=indirectGreatCircleProbSnyder(latLonStart(1:2,k),latLonEnd(1:2,k),r,N);
            if(~isempty(N))
                latLonWaypoints(:,:,k)=latLonWaypointsCur;
            end
        end
    otherwise
        error('Unknown algorithm chosen.')
end

if(~isempty(N)&&noJumps&&N>0)
    cumOffset=0;
    for k=1:numPts
        for j=2:(N+2)
            diffVal=latLonWaypoints(2,j)+cumOffset-latLonWaypoints(2,j-1);
            if(diffVal>pi)
                cumOffset=cumOffset-2*pi;
            elseif(diffVal<-pi)
                cumOffset=cumOffset+2*pi;
            end
            latLonWaypoints(2,j)=latLonWaypoints(2,j)+cumOffset;
        end
    end
end
end

function [azStart,dist,azEnd,latLonWayPoints]=indirectGreatCircleProbChen(latLonStart,latLonEnd,r,N,waypointType)
%%INDIRECTGREATCIRCLEPROBCHEN Solve the indirect great circle problem using
%                   the algorithm of [1], except the angular difference D
%                   is obtained using the angBetweenVecs given F and T.
%
%Starting latitude values within 2^24*eps(pi/2); of the pole will be
%clipped to be 2^24*eps(pi/2); away from the pole. This reduces problems
%with singularities at the poles at the cost of a loss of accuracy compared
%to elsewhere on the globe.
%
%REFERENCES:
%[1] C.-L. Chen, P.-F. Liu, and W.-T. Gong, "A simple approach to great
%    circle sailing: The COFI method," The Journal of Navigation, vol. 67,
%    no. 3, pp. 403-418, May 2014.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Latitude at the start.
LF=latLonStart(1);
%Deal with locations very close to the North or South pole by offsetting
%them just enough to avoid finite precision issues.
epsVal=2^24*eps(pi/2);
if(abs(LF-pi/2)<epsVal)
    LF=pi/2-epsVal;
elseif(abs(LF+pi/2)<epsVal)
    LF=pi/2+epsVal;
end

%Latitude at the destination.
LT=latLonEnd(1);

%Difference in longitude. If DLo is positive, then the shortest path is
%going East. If it is negative, then the shortest path is going west.
DLo=wrapRange(latLonEnd(2)-latLonStart(2),-pi,pi);

goingWest=(DLo<0);

sinLF=sin(LF);
cosLF=cos(LF);
sinLT=sin(LT);
cosLT=cos(LT);

lonStart=latLonStart(2);
cosLS=cos(lonStart);
sinLS=sin(lonStart);

lonEnd=latLonEnd(2);
cosLE=cos(lonEnd);
sinLE=sin(lonEnd);

%The Cartesian starting point (unit sphere).
F=[cosLS*cosLF;
   sinLS*cosLF;
   sinLF];

%The Cartesian ending point (unit sphere).
T=[cosLE*cosLT;
   sinLE*cosLT;
   sinLT];

%This is more accurate than using the cosine formula that is implicit in
%Equations 13 and 15 in [1].
D=angBetweenVecs(F,T);
if(D==0)
   %Special case: The starting and ending points are the same. 
    azStart=NaN;%Undefined direction.
    dist=0;
    azEnd=NaN;%Undefined direction.
    latLonWayPoints=repmat(latLonStart,[1,N+2]);
    return
end

cosD=cos(D);
sinD=sin(D);

%From Equation 17 in [1].
cosC=(sinLT-sinLF*cosD)/(cosLF*sinD);

%Deal with finite precision issues.
cosC=min(1,max(-1,cosC));
C=acos(cosC);

%The sign adjustment for going west is not mentioned in [1].
%Alternatively, one could get sinC=cosLT*sin(DLo)/sinD; and then use the
%atan2 function. 
if(goingWest)
    azStart=-C;
else
    azStart=C;
end
dist=D*r;

%%%%%%%%%%%%%
%%%%The Azimuth at the End
%%%%%%%%%%%%%

%To get the azimuth at the end, use Equation 17 in [1] switching the
%beginning and the end; this flips the sign of D, so the sign of sinD
%flips. Then, flip the direction of the result 180 degrees if we are going
%West.
cosCEnd=-(sinLF-sinLT*cosD)/(cosLT*sinD);
%Deal with finite precision issues.
cosCEnd=max(-1,min(1,cosCEnd));
azEnd=acos(cosCEnd);
if(goingWest)
    azEnd=-azEnd;
end

%%%%%%%%%%%%%
%%%%Waypoints
%%%%%%%%%%%%%
if(~isempty(N))
    %If waypoints are desired.
    latLonWayPoints=zeros(2,N+2);
    latLonWayPoints(:,1)=latLonStart;
    latLonWayPoints(:,N+2)=latLonEnd;
    if(waypointType==0)
        %Waypoints uniformly spaced in distance.
        step=D/(N+1);
        DFX=step:step:(step*N);

        cosDFX=cos(DFX);
        sinDFX=sin(DFX);

        %Equation 19 in [1].
        sinLX=sinLF.*cosDFX+cosLF.*sinDFX.*cosC;
        cosLX=sqrt(1-sinLX.^2);%Equivalent to cosLX=cos(asin(sinLX));

        %Equation 21 in [1].
        cosDLoFX=(cosDFX-sinLF.*sinLX)./(cosLF.*cosLX);

        %Deal with finite precision issues.
        cosDLoFX=max(-1,min(1,cosDLoFX));
        if(goingWest)
            lonPoints=wrapRange(latLonStart(2)-acos(cosDLoFX),-pi,pi);
        else
            lonPoints=wrapRange(latLonStart(2)+acos(cosDLoFX),-pi,pi); 
        end
        latLonWayPoints(:,2:(N+1))=[asin(sinLX);lonPoints];
    else
        %Waypoints uniformly spaced in longitude.
        DLonFT=wrapRange(lonEnd-lonStart,-pi,pi);
       
        if(DLonFT==0)
            %Meridonal trajectory.
            latLonWayPoints=[latLonStart,latLonEnd];
            return 
        end
       
        deltaSign=sign(DLonFT);
        
        step=DLonFT/(N+1);
        DLoFX=abs(step:step:(step*N));
        cosDLoFX=cos(DLoFX);
        sinDLoFX=sin(DLoFX);
 
        sinC=sin(C);
        
        %Equation 38 in [1].
        tanLX=(cosC*sinDLoFX+sinLF*sinC*cosDLoFX)/(cosLF*sinC);
        
        latLonWayPoints(:,2:(N+1))=[atan(tanLX);wrapRange(lonStart+deltaSign*DLoFX,-pi,pi)];
    end
else
    latLonWayPoints=[];
end
end

function [azStart,dist,azEnd,latLonWayPoints]=indirectGreatCircleProbCrouse(latLonStart,latLonEnd,r,N)
%%INDIRECTGREATCIRCLEPROBCROUSE Solve the indirect great circle problem
%               using the algorithm of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Singularity-free great-circle sailing," Naval Research
%    Laboratory, Washington, DC, Tech. Rep. NRL/5340/MR-2021/4, 26
%    Jul. 2021.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Convert the latitude-longitude values to unit vectors.
u1=spher2Cart([latLonStart(2);latLonStart(1)]);
u2=spher2Cart([latLonEnd(2);latLonEnd(1)]);

%Rotate the vectors both onto the equator.
[vec1Rot,vec2Rot,R]=rot2Vecs2CommonPlane(u1,u2);

%The geodesic goes along the equator in the rotated coordinate system in
%the tangent plane. We need to find the shortest direction around the
%equator. To do this, we shall get the azimuth angles.
az1=atan2(vec1Rot(2),vec1Rot(1));
az2=atan2(vec2Rot(2),vec2Rot(1));
DLo=wrapRange(az2-az1,-pi,pi);
goingWest=(DLo<0);

%The distance traveled between the points.
dist=abs(DLo)*r;

%Unit vectors in the traveling directions are thus either the positive or
%the negative of the East unit direction vector at the starting and ending
%points.
uENU1=getENUAxes([0;az1],false,1,0);
uENU2=getENUAxes([0;az2],false,1,0);

if(goingWest)
    uDir1=-uENU1(:,1);
    uDir2=-uENU2(:,1);
else
    uDir1=uENU1(:,1);
    uDir2=uENU2(:,1);
end

%Now, we shall undo the rotations to get the unit direction vectors in the
%original Cartesian coordinate system.
uDir1=R'*uDir1;
uDir2=R'*uDir2;

%Given the unit vectors in the original Cartesian coordinate system, we
%need to get the azimuth angle. The azimuth angle is in terms of degrees
%East of North, so we need to project into the East-North plane.
%Get the new set of ENU basis vectors in the  original coordinate system
uENU1=getENUAxes(latLonStart,false,1,0);
uENU2=getENUAxes(latLonEnd,false,1,0);

EProj=dot(uDir1,uENU1(:,1));
NProj=dot(uDir1,uENU1(:,2));
azStart=atan2(EProj,NProj);

EProj=dot(uDir2,uENU2(:,1));
NProj=dot(uDir2,uENU2(:,2));
azEnd=atan2(EProj,NProj);

if(~isempty(N))%If waypoints are desired.
    distPts=linspace(0,dist,N+2);
    latLonWayPoints=directGreatCircleProb(latLonStart,azStart,distPts,r);
    %This ensures that the longitudes match even if the starting and/or
    %stopping points are at the poles.
    latLonWayPoints(:,1)=latLonStart;
    latLonWayPoints(:,N+2)=latLonEnd;
else
    latLonWayPoints=[];
end
end

function [azStart,dist,azEnd,latLonWayPoints]=indirectGreatCircleProbSnyder(latLonStart,latLonEnd,r,N)
%%INDIRECTGREATCIRCLEPRBSNYDER This function calls greatCircleDistance with
%   algorithm 0 selected, which should be the most stable solution. Then,
%   it calls greatCircleAzimuth with algorithm 2 selected, which is the
%   very simple formula in Equation 5-4b in [1].
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections-a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Latitude at the start.
phi1=latLonStart(1);
%Latitude at the destination.
phi2=latLonEnd(1);

lambda1=latLonStart(2);
lambda2=latLonEnd(2);
cosPhi1=cos(phi1);
sinPhi1=sin(phi1);
cosPhi2=cos(phi2);
sinPhi2=sin(phi2);

cosLambda1=cos(lambda1);
sinLambda1=sin(lambda1);
cosLambda2=cos(lambda2);
sinLambda2=sin(lambda2);
deltaLambda=lambda2-lambda1;
cosDeltaLambda=cos(deltaLambda);
sinDeltaLambda=sin(deltaLambda);

%The Cartesian starting point (unit sphere).
F=[cosLambda1*cosPhi1;
   sinLambda1*cosPhi1;
   sinPhi1];
%The Cartesian ending point (unit sphere).
T=[cosLambda2*cosPhi2;
   sinLambda2*cosPhi2;
   sinPhi2];
dist=r*angBetweenVecs(F,T);

azStart=atan2(cosPhi2*sinDeltaLambda,cosPhi1*sinPhi2-sinPhi1*cosPhi2*cosDeltaLambda);
azEnd=atan2(-cosPhi1*sinDeltaLambda,cosPhi2*sinPhi1-sinPhi2*cosPhi1*cosDeltaLambda);
azEnd=wrapRange(azEnd+pi,-pi,pi);

if(~isempty(N))%If waypoints are desired.
    distPts=linspace(0,dist,N+2);
    latLonWayPoints=directGreatCircleProb(latLonStart,azStart,distPts,r);
    %This ensures that the longitudes match even if the starting and/or
    %stopping points are at the poles.
    latLonWayPoints(:,1)=latLonStart;
    latLonWayPoints(:,N+2)=latLonEnd;
else
    latLonWayPoints=[];
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
