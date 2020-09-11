function [azimuth, dist]=indirectRhumbProblem(latLonStart,latLonEnd,height,useHeightApprox,a,f,numSteps4Circ)
%INDIRECTRHUMBPROBLEM Given a starting and ending latitude and longitude on
%                     a reference ellipsoid, determine the heading (in
%                     radians East of North) and the distance one must
%                     travel on the shortest constant-heading course to go
%                     from the starting point to the stopping point. A
%                     constant heading course follows a rhumb line
%                     (loxodrome) and is usually not the shortest path
%                     between two points.
%
%INPUTS: latLonStart A 2X1 vector of the starting ellipsoidal latitude
%                    (North) and longitude (East) in radians. This cannot
%                    be a pole.
%          latLonEnd A 2X1 vector of the ending ellipsoidal latitude
%                    (North) and longitude (East) in radians.
%             height The height above the reference ellipsoid at which the
%                    trajectory should be determined. This changes the
%                    distance traveled, but not the azimuthal angle of
%                    departure. If this parameter is omitted, then the
%                    default value of 0 is used.
%    useHeightApprox If true, and the height is not zero, then an
%                    approximation is made for how dist scales with
%                    altitude. Specifically, an equiatorial trajectory will
%                    scale as (a+height)/a. Thus, this scaling factor is
%                    applied to any trajectory to scale dist with altitude.
%                    If height is false, a significantly slower iterative
%                    optimization technique is used. The default value is
%                    true. The difference made when useHeightApprox=false is
%                    can generally be assumed to be less than 80m. This
%                    parameter is ignored if height=0.
%                  a The semi-major axis of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84SemiMajorAxis is used.
%                  f The flattening factor of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84Flattening is used.
%      numSteps4Circ If height!=0 then an algorithm propagating a state
%                    in ECEF coordinates around the curved Earth is used to
%                    solve the direct geodetic problem as a step in solving
%                    the indirect geodetic problem. This parameter
%                    determines the number of steps that would be needed in
%                    the direct geodetic problem for a target that
%                    circumnavigates the globe around the equator. The
%                    default value if this parameter is not provided is
%                    2000. A value of 6000 appears to be about the best
%                    number for overall precision. Reducing the number of
%                    steps will speed up the function. This parameter is not
%                    used if height=0.
%
%OUTPUTS: azimuth The constant heading in radians East of North that one
%                 must travel to go on a constant-heading course from
%                 latLonStart to latLonEnd.
%            dist The distance that one must travel on a constant-heading
%                 course to go from latLonStart to latLonEnd.
%
%If height=0, the algorithm is mostly taken from [1]. However, a formula
%using isometric latitudes, which are described in Chapter 3 of [2] to get
%the azimuth angle was used, because it is simpler. The formula is also
%explicitly mentioned in Equation 3 of [3]. However, the expression for
%computing the distance from that paper is only for a sphere, not for an 
%ellipsoid, which is why the Carlton-Wippern distance computation using an
%incomplete elliptic integral of the second kind is preferred.
%
%When the azimuth found by the technique is very close to +/-pi/2 (when one
%is traveling at nearly a constant latitude), the distance computation
%switches to assume that the latitude is indeed constant even if
%latLonStart(1)!=latLonStart(2). This avoid precision problems that arise
%as a very small number is multiplied by a very large number. However, this
%reduces the accuracy of the method.
%
%Generally, calling directRhumbProblem or directRhumbProbGen with the
%azimuth and dist returned by this function will return latLonStart.
%However, if the stopping point is at a pole, then directRhumbProblem will
%correctly return a polar location, but the longitude will generally be
%wrong.
%
%When height!=0, the algorithm can be significantly slower if no
%approximation is used. Around the equator, the distance scales as
%(a+height)/a*dist as one changes the ellipsoidal height. This scaling
%applied to dist in  non-equatorial trajectories is the approximation if
%useHeightApprox=true. When useHeightApprox=false, the approximate value
%is used to determine the bounds around which the fminbnd function
%searches. The maximum error in the approximation is expected to be less
%than 80m. The search region for the value of dist used in the fminbnd
%function was set to 0.9*dist to 1.1*dist, where dist is the distance
%obtained after scaling the distance from the zero-altitude solution.
%
%EXAMPLE:
%A trajectory that crosses the international date line and and goes from
%the Northern hemisphere to the southern hemisphere. We also compute the
%reverse path and show that the azimuth angles in each direction are
%consistent with each other. We then plot the trajectory on an image of the
%spherical Earth.
% N=100;
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% 
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% [azimuth,dist]=indirectRhumbProblem(latLonStart,latLonEnd);
% [azEnd,distRev]=indirectRhumbProblem(latLonEnd,latLonStart);
% 
% %If the forward and reverse estimates agree, then these values will
% %ideally be zero.
% azimuth-(-azEnd)
% dist-distRev
% 
% distVals=linspace(0,dist,N);
% latLonWayPoints=directRhumbProblem(latLonStart,azimuth,distVals);
% 
% %Show that the approximate direct algorithm reaches nearly the same
% %endpoint as the indirect algorithm.
% xEndWay=ellips2Cart([latLonWayPoints(:,end);0]);
% xEnd=ellips2Cart([latLonEnd;0]);
% max(abs(xEndWay-xEnd))
% max(abs(wrapRange(latLonWayPoints(:,end)-latLonEnd,-pi,pi)))
% 
% xStartCart=ellips2Cart([latLonStart;0]);
% xEndCart=ellips2Cart([latLonEnd;0]);
% %The path is displayed slightly above the Earth's surface to make it
% %easier to see.
% pathPoints=ellips2Cart([latLonWayPoints;0.02*ones(1,N)]);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([]);
% scatter3(xStartCart(1),xStartCart(2),xStartCart(3),100,'filled')
% scatter3(xEndCart(1),xEndCart(2),xEndCart(3),100,'filled')
% plot3(pathPoints(1,:),pathPoints(2,:),pathPoints(3,:),'-r','linewidth',4)
%
%REFERENCES:
%[1] K. C. Carlton-Wippern, "On loxodromic navigation," Journal of
%    Navigation, vol. 45, no. 2, pp. 292-297, May 1992.
%[2] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[3] J. Alexander, "Loxodromes: A rhumb way to go," Mathematics Magazine,
%    vol. 77, no. 5, pp. 349-356, Dec. 2004.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(numSteps4Circ))
   numSteps4Circ=2000; 
end

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(useHeightApprox))
   useHeightApprox=true; 
end

if(nargin<3||isempty(height))
   height=0; 
end

%Extract the components
latStart=latLonStart(1);
lonStart=latLonStart(2);
latEnd=latLonEnd(1);
lonEnd=latLonEnd(2);

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

%Convert the ellipsoidal latitudes to reduced co-latitudes. A co-latitude
%is pi/2 minus the latitude. Also, 
nu1=pi/2-ellipsLat2ReducedLat(latStart,f);
nu2=pi/2-ellipsLat2ReducedLat(latEnd,f);

%Though not mentioned in the above papers, the difference in the longitudes
%must be wrapped to the range of -pi/pi or else one will get useless
%results when crossing the -pi/pi boundary.
num=wrapRange(lonEnd-lonStart,-pi,pi,false);
%Equation 11 in the paper provides an expression to get the azimuth.
%however, it is simpler if one just uses isometric latitudes.
val1=ellipsLat2IsoLat(latStart,f);
val2=ellipsLat2IsoLat(latEnd,f);

if(~isfinite(val1))
    warning('The starting point is located at a geographic pole. Azimuth values will be inaccurate.')
end

if(isfinite(val1)||isfinite(val2))
    %This will properly return 0 and pi for infinite values of val1 when
    %val2 is finite and vice versa, which corresponds to headings to or
    %from a pole.
    denom=val2-val1;
    azimuth=atan2(num,denom);
else%If neither is finite, then that is the case where one is going from
    %pole to pole. In such an instance, just set the heading to 0, if going
    %North, and pi, if going South.
    
    if(latStart>latEnd)
        azimuth=pi;
    else
        azimuth=0;
    end
end

%The distance 
if(abs(abs(azimuth)-pi/2)>2e-8)
    %Equation 12 in the paper.
    dist=a*abs(sec(azimuth))*abs(ellipIntInc2Kind(nu2,e^2)-ellipIntInc2Kind(nu1,e^2));
else
    %Equation 14b in the paper.
    dist=a*abs(sin(nu1))*abs(wrapRange(lonEnd-lonStart,-pi,pi));
end

%If a non-zero height is given, then iterate over the direct rhumb
%problem at height to determine the solution.
if(height~=0)
    endCart=ellips2Cart([latLonEnd;height],a,f);

    %The approximate scaling for the height.
    dist=((a+height)/a)*dist;
    
    %If a computationally-intensive but more precise algorithm to search
    %for the true distance at altitude should be used instead of a
    %simple approximation of scaling the distance.
    if(useHeightApprox==false)
        %Assume that the correct height-adjusted distance is within 10% of the
        %scaled distance value.
        distFun=@(distCur)distCostFunc(distCur,endCart,latLonStart,azimuth,height,a,f,numSteps4Circ);
        dist=fminbnd(distFun,0.9*dist,1.1*dist);
    end
end
end

function cost=distCostFunc(distCur,endCart,latLonStart,azStart,height,a,f,numSteps4Circ)
    latLonCalc=directRhumbProbGen(latLonStart,azStart,distCur,height,false,a,f,numSteps4Circ);
    cost=norm(ellips2Cart([latLonCalc;height],a,f)-endCart);
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
