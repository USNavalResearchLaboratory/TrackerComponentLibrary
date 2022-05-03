function latLonEnd=directRhumbProbGen(latLonStart,azimuth,dist,height,useHeightApprox,a,f,numSteps4Circ)
%%DIRECTRHUMBLINEPROBGEN Given a starting location near the surface of the
%                    reference ellipsoid as well as a constant heading a
%                    distance to travel, determine the location of a vessel
%                    traveling at the constant heading after traveling the
%                    distance. A constant heading trajectory is a rhumb
%                    line (loxodrome) and is usually not the shortest path
%                    between two points. This function is general in that
%                    it will also solve the direct rhumb problem at a fixed
%                    height above the reference ellipsoid, if desired.
%
%INPUTS: latLonStart A 2X1 vector of the starting ellipsoidal latitude
%                    (North) and longitude (East) in radians.
%            azimuth The constant heading that is to be traveled in
%                    radians East of North.
%               dist The distance that is to be traveled on a constant
%                    heading course.
%             height The heights above the reference ellipsoid at which the
%                    trajectory should be determined. If this parameter
%                    is omitted, then the default value of 0 is used.
%    useHeightApprox If true, and the height is not zero, then a distance
%                    scaling height approximation is used. An equatorial
%                    trajectory traveling a distance dist at a height h
%                    above the reference ellipsoid will only travel
%                    a/(a+h)*dist over the surface. If useHeightApprox is
%                    true, then this scaling applied to any trajectry
%                    (including non-equatorial). If useHeightApprox is
%                    false, a slow iterative optimization used. The default
%                    value if omitted or an empty matrix is passed is
%                    false. This parameter is ignored if height=0.
%                  a The semi-major axis of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84SemiMajorAxis is used.
%                  f The flattening factor of the reference ellipsoid. If
%                    this argument is omitted, the value in
%                    Constants.WGS84Flattening is used.
%      numSteps4Circ If height!=0, and useHeightApprox=false, then an
%                    algorithm propagating a state in ECEF coordinates
%                    around the curved Earth is used. This parameter
%                    determines the number of steps that would be needed
%                    for a target that circumnavigates the globe around the
%                    equator. The default value if this parameter is not
%                    provided is 2000. A value of 6000 appears to be about
%                    the best number for overall precision. Reducing the
%                    number of steps will speed up the function. This
%                    parameter is not used if height=0.
%
%OUTPUTS: latLonEnd The ellipsoidal latitude and longitude that one will be
%                    after starting at latLonStart and traveling a distance
%                    of dist at a constant heading on azimuth.
%
%For zero-altitude, the directRhumbProb function is called. The solution
%for a non-zero height is based on [1]. The function is significantly
%faster if a zero hight is used or if useHeightApprox=true, in which case
%the approximation described in [2] is used.
%
%The algorithm is generally the opposite of indirectRhumbProblem, except
%when a starting or ending point is at a geographic pole, as discussed in
%the comments to the indirectRhumbProblem function.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%[2] D. F. Crouse, "Worldwide Ground Target State Propagation," in
%    Proceedings of the 23rd International Conference on Information
%    Fusion, 6-9 Jul. 2020.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
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

if(nargin<4||isempty(height))
    height=0; 
end

if(useHeightApprox)
    %Scale the distance based on the height and then solve the rumb problem
    %with the scaled distance.
    dist=a/(a+height)*dist;
    height=0; 
end

if(height==0)
    latLonEnd=directRhumbProblem(latLonStart,azimuth,dist,0,a,f);
    return
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
