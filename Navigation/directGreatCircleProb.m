function latLonPoints=directGreatCircleProb(latLonStart,azimuth,dist,r,algorithm,distType)
%%DIRECTGREATCIRCLEPROB Solve the direct great circle problem. That is,
%                      given a starting point, an azimuth East of true
%                      North on a spherical Earth and a distance, determine
%                      the location of a ship traveleing the given distance
%                      along the great circle with the given azimuth. This
%                      type of navigation is sometimes referred to as
%                      "great circle sailing." When used for navigation, it
%                      is less accurate than an ellipsoidal Earth
%                      approximation.
%
%INPUTS: latLonStart The 2X1 initial point given in latitude and longitude
%                    in radians in the format [latitude;longitude]
%                    (on a reference sphere, latitude is spherical
%                    elevation; longitude is spherical azimuth).
%            azStart The scalar forward azimuth at the starting point in
%                    radians East of true North on the reference sphere.
%               dist The NX1 or 1XN set of distances on the sphere-
%                    approximated Earth starting from latLonStart where one
%                    wishes to find the stopping point. If distType=1, this
%                    specifies offsets in longitude from the starting point
%                    instead of distances.
%                  r The assumed radius of the spherical Earth model. If
%                    omitted or an empty matrix is passed, the default of
%                    osculatingSpher4LatLon(latLonStart)is used. If
%                    traveling a height h above the reference sphere, then
%                    one should use something like r=h+rEarth.
%          algorithm An optional parameter selecting the algorithm to use.
%                    Possible values are:
%                    0 Use the algorithm of [2], which avoids issues with
%                      singularities at the geographic poles and is
%                      consistent with the directions returned by
%                      indirectGreatCircleProb at the poles.
%                    1 (The default if omitted or an empty matrix is
%                      passed if distType=1) Use the COFI algorithm of [1].
%                      Latitudes within 2^24*eps(pi/2) of +/-pi/2 (the
%                      North and South poles) will be clipped to that
%                      bound. This avoids a singularity at the poles.
%                    2 (The default if omitted or an empty matrix is
%                      passed if distType=0) Use 20-14, 20-15, 20-18 and
%                      25-15 in Chapter 25 of [3]. 
%           distType This specifies the format of dist. Possible values
%                    are:
%                    0 (The default if omitted or an empty matrix is
%                      passed) dist holds distances along the sphere.
%                    1 dist holds offsets in longitude in the form of the
%                      starting points. If meridional trajectories are
%                      given, an error will be raised. This can only be
%                      used if algorithm=011.
%
%OUTPUTS: latLonPoints A 2XN matrix of geocentric latitude and longitudes
%                      of the final points of the spherical geodesic
%                      trajectory given in radians as [latitude;longitude].
%
%For the spherical model, the latitude is taken as the elevation and the
%longitude as the azimuth on the sphere (we assume coordinate system 0 in
%the spher2Cart function). 
%
%In algorithm 0, the starting point is rotated to [1;0;0]. Then a rotation
%is computed to make the direction point at [0;1;0]. This means that the
%trajectory goes East on the equator in the rotated coordinate system. The
%longitudinal offset is just the distance divided by the radius. Given
%Cartesian points in the transformed coordinate system, they are just
%rotated back to the global coordinate system.
%
%Algorithm 1 uses Equations 19 and 21 in [1] is used. 
%
%See the comments in indirectGreatCircleProb for additional comments on
%great circle navigation.
%
%EXAMPLE 1:
%This is the trajectory from Example 2 in [1]. We show that after solving
%the indirect great circle problem, the trajectory obtained leads to the
%desired endpoint when given to directGreatCircleProb.
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% [azStart,dist]=indirectGreatCircleProb(latLonStart,latLonEnd);
% latLonPoint=directGreatCircleProb(latLonStart,azStart,dist);
% max(abs(latLonPoint-latLonEnd))
%One will see that the two values are the same, as expected.
%
%EXAMPLE 2:
%This algorithm is an example of plotting by latitude and it also
%demonstrates that the trajectory can go further than halfway around the
%world. Here, it goes 350 degrees around and the points are plotted.
% N=500;
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latLonStart=[latStart;lonStart];
% 
% azStart=55*(pi/180);
% deltaLon=350*(pi/180);%Go almost all the way around the world.
% lonVals=linspace(0,deltaLon,N);
% 
% latLonPoints=directGreatCircleProb(latLonStart,azStart,lonVals,[],1,1);
% %Convert the points to Cartesian to plot. Plot then slightly above the
% %surface so that the line shows up better.
% xyzPoints=ellips2Cart([latLonPoints;0.02*ones(1,N)],1,0);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([],1,0);
% plot3(xyzPoints(1,:),xyzPoints(2,:),xyzPoints(3,:),'-r','linewidth',4)
%
%EXAMPLE 3:
%This example generates many random pairs of points and then uses
%indirectGreatCircleProb to get a distance and heading. Then, this function
%inverts the problem using directGreatCircleProb and the maximum l2 norm of
%the difference between the computed and original endpoings in latitude and
%longitude are shown. They should be on the order of 1-12, which is
%rasonable given finite precision limitiations.
% numMCRuns=1e4;
% maxAbsDiff=0;
% for k=1:numMCRuns
%     latLonRef=[UniformD.rand(1,[-pi/2;pi/2]);
%                 UniformD.rand(1,[-pi;pi])];
%     latLonPt=[UniformD.rand(1,[-pi/2;pi/2]);
%                 UniformD.rand(1,[-pi;pi])];
%     rE=osculatingSpher4LatLon(latLonRef);
%     [azStart,distVal]=indirectGreatCircleProb(latLonRef,latLonPt,rE);
%     latLonBack=directGreatCircleProb(latLonRef,azStart,distVal,rE);
%     diffVal=latLonBack-latLonPt;
%     %Deal with getting possible offsets near 2*pi.
%     diffVal(2)=wrapRange(diffVal(2),-pi,pi);
%     maxAbsDiff=max(maxAbsDiff,norm(diffVal));
% end
% maxAbsDiff
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

if(nargin<6||isempty(distType))
    distType=0;
end

if(nargin<5||isempty(algorithm))
    if(distType==1)
        algorithm=1;
    else
        algorithm=2;
    end
end

if(nargin<4||isempty(r))
    r=osculatingSpher4LatLon(latLonStart);
end

dist=dist(:).';

switch(algorithm)
    case 0
        if(distType==1)
            error('This algorithm only supports distances, not offsets in latitude.')
        end

        latLonPoints=directGreatCircleProbCrouse(latLonStart,azimuth,dist,r);
    case 1
        latLonPoints=directGreatCircleProbChen(latLonStart,azimuth,dist,r,distType);
    case 2
        latLonPoints=directGreatCircleProbSnyder(latLonStart,azimuth,dist,r);
    otherwise
        error('Unknown algorithm specified.')
end

%Deal with errors that arise with zero distances.
sel=(dist==0);
latLonPoints(:,sel)=latLonStart(:,sel);

end

function latLonPoints=directGreatCircleProbChen(latLonStart,azimuth,dist,r,distType)
%%DIRECTGREATCIRCLEPROBCHEN Solve the direct great circle problem using
%                           Equations 19 and 21 in [1] if given distances
%                           and using Equation 38 in 1 if given offsets in
%                           longitude.
%
%REFERENCES:
%[1] C.-L. Chen, P.-F. Liu, and W.-T. Gong, "A simple approach to great
%    circle sailing: The COFI method," The Journal of Navigation, vol. 67,
%    no. 3, pp. 403-418, May 2014.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

cosC=cos(azimuth);
if(wrapRange(azimuth,-pi,pi)<0)
    goingWest=true;
else
    goingWest=false;
end

LF=latLonStart(1);
%Deal with locations very close to the North or South pole by offsetting
%them just enough to avoid finite precision issues.
epsVal=2^24*eps(pi/2);
if(abs(LF-pi/2)<epsVal)
    LF=pi/2-epsVal;
elseif(abs(LF+pi/2)<epsVal)
    LF=pi/2+epsVal;
end

cosLF=cos(LF);
sinLF=sin(LF);

if(distType==0)
    %If the dist holds distances along the surface of the sphere.
    DFX=dist(:).'/r;

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
    latLonPoints=[asin(sinLX);lonPoints];
else
    %If dist holds offsets in longitude.
    
    sinC=sin(azimuth);
    
    if(sinC==0)
        error('Offsets in longitude cannot be used for meridional trajectories.') 
    end
    
    DLoFX=dist(:).';
    
    cosDLoFX=cos(DLoFX);
    sinDLoFX=sin(DLoFX);
        
    %Equation 38 in [1].
    tanLX=(cosC*sinDLoFX+sinLF*sinC*cosDLoFX)/(cosLF*sinC);
        
    latLonPoints=[atan(tanLX);wrapRange(latLonStart(2)+DLoFX,-pi,pi)];
end
end

function latLonPoints=directGreatCircleProbCrouse(latLonStart,azimuth,dist,r)
%%DIRECTGREATCIRCLEPROBCROUSE Solve the direct great circle problem using
%                           the algorithm of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Singularity-free great-circle sailing," Naval Research
%    Laboratory, Washington, DC, Tech. Rep. NRL/5340/MR-2021/4, 26
%    Jul. 2021.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

uENUOrig=getENUAxes(latLonStart,false,1,0);
%3D unit direction vector in the global coordinate system
uDirOrig=sin(azimuth)*uENUOrig(:,1)+cos(azimuth)*uENUOrig(:,2);

az=latLonStart(2);
el=latLonStart(1);
%Rotation matrix such that the point latLonStart is at [1;0;0].
Ryz=Euler2Ang2RotMat(el,-az,'yz');
%Rotate the heading direction vector.
uDirRot=Ryz*uDirOrig;

%Extract North and East in the rotated local coordinate system.
NProj=uDirRot(3);
EProj=uDirRot(2);
azLocal=atan2(NProj,EProj);

%Rotation matrix to rotate uDirRot to be [0;1;0]. In the transformed
%system, the path goes East long the equator.
Rx=Euler1Ang2RotMat(-azLocal,'x');
%Combined rotation matrix.
Rxyz=Rx*Ryz;

dist=dist(:).';
N=length(dist);

%Longitude points.
lon=dist/r;

%Get Cartesian unit vectors (spherical to Cartesian conversion)
uVecWaypoints=zeros(3,N);
uVecWaypoints(1,:)=cos(lon);
uVecWaypoints(2,:)=sin(lon);
uVecWaypoints(3,:)=0;

%Rotate the waypoints into the global coordinate system:
uVecWaypoints=Rxyz'*uVecWaypoints;

%Convert the waypoints into spherical coordinates.
az=atan2(uVecWaypoints(2,:),uVecWaypoints(1,:));
el=asin(uVecWaypoints(3,:));
latLonPoints=[el;az];

end

function latLonPoints=directGreatCircleProbSnyder(latLonStart,azimuth,dist,r)
%%DIRECTGREATCIRCLEPROBSNYDER Solve the direct great circle problem using
%          the inverse equations for azimuthal equidistant coordinates on a
%          sphere given in Chapter 25 of of [1].
%
%REFERENCES:
%[3] J. P. Snyder, "Map projections-a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

phi0=latLonStart(1);
lambda0=latLonStart(2);
numDist=length(dist);

xy=pol2Cart([dist;azimuth*ones(1,numDist)],1);
x=xy(1,:);
y=xy(2,:);
%The max should deal with finite precision limitations.
rho=sqrt(max(x.^2+y.^2,0));
c=rho./r;

cosC=cos(c);
sinC=sin(c);

latLonPoints=[asin(cosC.*sin(phi0)+(y.*sinC.*cos(phi0)./rho));
              lambda0+atan2(x.*sinC,rho.*cos(phi0).*cosC-y.*sin(phi0).*sinC)];

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
