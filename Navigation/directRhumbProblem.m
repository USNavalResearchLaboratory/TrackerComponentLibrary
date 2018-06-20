function latLonEnd=directRhumbProblem(latLonStart,azimuth,dist,height,a,f,numSteps4Circ)
%%DIRECTRHUMBLINEPROBLEM Given a starting location on the surface of the
%                        reference ellipsoid as well as a constant heading
%                        a distance to travel, determine the location of a
%                        vessel traveling at the constant heading after
%                        traveling the distance. A constant heading
%                        trajectory is a rhumb line (loxodrome) and is
%                        usually not the shortest path between two points.
%
%INPUTS: latLonStart A 2X1 vector of the starting ellipsoidal latitude
%                    (North) and longitude (East) in radians.
%            azimuth The constant heading that is to be traveled in
%                    radians East of North.
%               dist The distance that is to be traveled on a constant
%                    heading course.
%             height The heights above the reference ellipsoid at which the
%                    trajectoriey should be determined. If this parameter
%                    is omitted, then the default value of 0 is used.
%                  a The semi-major axis of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84SemiMajorAxis is used.
%                  f The flattening factor of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84Flattening is used.
%      numSteps4Circ If height!=0, then an algorithm propagating a state in
%                    ECEF coordinates around the curved Earth is used. This
%                    parameter determines the number of steps that would be
%                    needed for a target that circumnavigates the globe
%                    around the equator. The default value if this
%                    parameter is not provided is 2000. A value of 6000
%                    appears to be about the best number for overall
%                    precision. Reducing the number of steps will speed up
%                    the function. This parameter is not used if height=0.
%
%OUTPUTS: latLonEnd The ellipsoidal latitude and longitude that one will be
%                    after starting at latLonStart and traveling a distance
%                    of dist at a constant heading on azimuth.
%
%The algorithm for zero-altitude is taken from [1]. However, a formula
%using isometric latitudes, which are described in Chapter 3 of [2] to get
%the azimuth angle was used, because it is simpler. The formula is
%also explicitly mentioned in Equation 2 of [3]. However, the expression
%for computing the distance from [3] is only for a sphere, not for an
%ellipsoid, which is why the Carlton-Wippern distance computation using an
%incomplete elliptic integral of the second kind is preferred.
%
%The algorithm is generally the opposite of indirectRhumbProblem, except
%when a starting or ending point is at a geographic pole, as discussed in
%the comments to the indirectRhumbProblem function.
%
%The solution for a non-zero altitude is based on [4]. The function is
%significantly faster if a zero altitude is used.
%
%REFERENCES:
%[1] K. C. Carlton-Wippern, "On loxodromic navigation," Journal of 
%    Navigation, vol. 45, no. 2, pp. 292-297, May 1992.
%[2] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[3] J. Alexander, "Loxodromes: A rhumb way to go," Mathematics Magazine,
%    vol. 77, no. 5, pp. 349-356, Dec. 2004.
%[4] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(numSteps4Circ))
%This is the number of steps that would be taken to travel a distance of
%2*pi*a (circumnavigate the globe the long way). This is a design parameter
%that should be large enough to ensure that the results are accurate.
%However, if it is too large, finite precision errors will accumulate.
    numSteps4Circ=2000;
end

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(f))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(height))
   height=0; 
end

if(height==0)
    latStart=latLonStart(1);
    lonStart=latLonStart(2);

    %The first numerical eccentricity of the ellipsoid.
    e=sqrt(2*f-f^2);

    %Convert the ellipsoidal latitudes to a reduced co-latitude. A co-latitude
    %is pi/2 minus the latitude.
    nu1=pi/2-ellipsLat2ReducedLat(latStart,f);

    %Make sure that the azimuth is in the range of -pi to pi.
    azimuth=wrapRange(azimuth,-pi,pi,false);

    %If the azimuth is too close to due East/West, then use the solution for
    %due East/West to avoid numerical precision issues. Otherwise, use the
    %solution for other headings in general.
    if(abs(abs(azimuth)-pi/2)>2e-8)
        %The solution for general headings.

        %The following sets up and performs the iteration from Equation 13,
        %which inverts the incomplete elliptic function of the second kind
        %using Newton's method.
        nu2=nu1;
        if(-pi/2<azimuth&&azimuth<pi/2)
            signVal=1;
        else
            signVal=-1;
        end

        %It seems to converge in under 6 iterations
        numIter=6;
        for curIter=1:numIter
            num=abs(ellipIntInc2Kind(nu2,e^2)-ellipIntInc2Kind(nu1,e^2))-dist*abs(cos(azimuth))/a;
            denom=sqrt(1-e^2*sin(nu2)^2);

            nu2=nu2+signVal*num/denom;
        end

        %Convert the reduced co-latitude of the ending point to an ellipsoidal
        %latitude. 
        latEnd=reducedLat2EllipsLat(pi/2-nu2,f);

        %Next, we want to find the ending longitude. This could be done using
        %equation 11 in the Carlton-Wippern paper. However, we are going to do
        %it using isometric latitudes, since it is simpler.
        lonEnd=lonStart+tan(azimuth)*(ellipsLat2IsoLat(latEnd,f)-ellipsLat2IsoLat(latStart,f));

        %However, if the trajectory is too close to 0/ 180 degrees, then the
        %ellipsLat2IsoLat function might not be finite, in which case the
        %difference will be a NaN, rather than 0. Thus, we will check for that
        %case. In such an instance, the ending longitude is just set to the
        %starting latitude.
        if(~isfinite(lonEnd))
           lonEnd=lonStart; 
        end
    else
        %The solution for East-West headings. This inverts equation 14b to get
        %the ending longitude.

        %If it is going East (positive direction)
        if(abs(azimuth-pi/2)<abs(azimuth+pi/2))
            signVal=1;
        else
            signVal=-1;
        end
        lonEnd=lonStart+signVal*dist/(a*abs(sin(nu1)));

        %If one tries to have an East-West trajectory at the pole, then the
        %denominator above will be 0. This keeps it from returning a NaN.
        if(~isfinite(lonEnd))
           lonEnd=lonStart; 
        end

        latEnd=latStart;
    end
        %Make sure that the longitude is within the proper range.
        latLonEnd=[latEnd;wrapRange(lonEnd,-pi,pi,false)];
else
    %If the height is not zero, then use the algorithm for propagating a
    %level, non-maneuvering dynamic model along a rhumb line forward in
    %time.

    %This is the number of Runge-Kutta steps on a curved Earth that will be
    %performed.
    numSteps=ceil(numSteps4Circ*(dist/(2*pi*a)));
    stepSize=1/numSteps;
    
    %We have to calculate the initial location and heading in CARTESIAN
    %ECEF coordinates.
    xyzInit=ellips2Cart([latLonStart;height],a,f);
    
    %The initial heading in the LOCAL East-North-Up tangent-plane
    %coordinate system for level flight.
    uh=[sin(azimuth);cos(azimuth);0];
    
    %The speed of the target is such that it travels a distance of dist in
    %1 unit of time.
    speed=dist;
    
    %The initial target state with the position in the GLOBAL ECEF
    %coordinate system and the velocity in a LOCAL flat-Earth coordinate
    %system.
    xInit=[xyzInit;speed*uh];
    
    aDyn=@(x,t)aPoly(x,t,3);%A constant velocity model.

    %When traveling along a rhumb-line, the local coordinate system is
    %known at all times: it is the local ENU coordinate system. The local
    %basis vectors are deterministially known at all places on the Earth.
    uDet=@(x,t)getENUAxes(Cart2Ellipse(x(1:3),[],a,f));
    xList=RungeKCurvedAtTimes(xInit,[],[0;1],aDyn,uDet,stepSize,4);
    plhEnd=Cart2Ellipse(xList(1:3,end),[],a,f);
    latLonEnd=plhEnd(1:2);%The Cartesian location at the end of the
                          %trajectory.
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
