function [azStart,dist,azEnd]=indirectGeodeticProbGen(latLonStart,latLonEnd,height,useHeightApprox,a,f,numSteps4Circ)
%%INDIRECTGEODETICPROBGEN  Solve the indirect geodetic problem. That is,
%           given two points on an ellipsoidal Earth, find the initial
%           bearing and distance one must travel to take the shortest
%           (geodesic) path between the points. This function is general in
%           that it will also solve the indirect geodetic problem at a
%           fixed height above the reference ellipsoid either using a
%           simple approximation, or using a slow exact technique.
%
%INPUTS: latLonStart The initial point given in geodetic latitude and
%           longitude in radians of the format [latitude;longitude].
% latLonEnd The final point for the geodesic path given in geodetic
%           latitude and longitude in radians. latLonEnd has the same
%           format at latLonStart.
%    height The height above the reference ellipsoid at which the
%           trajectory should be determined. This changes the distance
%           traveled, but not the azimuthal angle of departure. If this
%           parameter is omitted or an empty matrix is passed, then the
%           default value of 0 is used.
% useHeightApprox If true, and the height is not zero, then an
%           approximation is made for how dist scales with altitude.
%           Specifically, an equatorial trajectory will scale as
%           (a+height)/a. Thus, this scaling factor is applied to any
%           trajectory to scale dist with altitude. If useHeightApprox is
%           false, a significantly slower iterative optimization technique
%           is used. The default value if omitted or an empty matrix is
%           passed is true. The difference made when useHeightApprox=false
%           can often be assumed to be less than 80m on the Earth. This
%           parameter is ignored if height=0.
%         a The semi-major axis of the reference ellipsoid (in meters). If
%           this argument is omitted or an empty matrix is passed, the
%           value in Constants.WGS84SemiMajorAxis is used.
%         f The flattening factor of the reference ellipsoid. If this
%           argument is omitted or an empty matrix is passed, the value in
%           Constants.WGS84Flattening is used.
% numSteps4Circ If height!=0 then an algorithm propagating a state in ECEF
%           coordinates around the curved Earth is used to solve the direct
%           geodetic problem as a step in solving the indirect geodetic
%           problem. This parameter determines the number of steps that
%           would be needed in the direct geodetic problem for a target
%           that circumnavigates the globe around the equator. The default
%           value if this parameter is omitted or an empty matrix is passed
%           is 2000. A value of 6000 appears to be about the best number
%           for overall precision. Reducing the number of steps will speed
%           up the function. This parameter is not used if height=0.
%
%OUTPUTS: azStart The forward azimuth at the starting point in radians East
%                 of true North on the reference ellipsoid. This is the
%                 initial heading one would travel to go between the two
%                 points.
%            dist The geodetic distance between the starting and stopping
%                 points in meters.
%           azEnd The forward azimuth at the ending point in radians East
%                 of true North on the reference ellipsoid.
%
%The algorithm initially solves the indirect geodetic problem on the
%surface of the reference ellipsoid using the function
%indirectGeodeticProb. Then, if a non-zero height is used, the function
%directGeodeticProbGen is iterated to obtain the correct solution.
%
%Assuming that the effects of height do not change the initial heading,
%but rather just the distance that must be traveled to solve the problem,
%the solution at non-zero heights above the reference ellipsoid can be
%obtained by solving the indirect geodetic problem for the heading and then
%performing a search over a range of values to solve the direct geodetic
%problem at the desired height. The direct geodetic problem at a nonzero
%height is solved using the function directGeodeticProblem and the
%minimization is performed using the fminbnd function.  This is the
%algorithm used if useHeightApprox=false. In practice, a simple
%approximation generally has sufficient accuracy, which is why the default
%useHeightApprox=true.
%
%Around the equator, the distance scales as (a+height)/a*dist as one
%changes the ellipsoidal height. This scaling applied to dist in 
%non-equatorial trajectories is the approximation if useHeightApprox=true.
%When useHeightApprox=false, the approximate value is used to determine the
%bounds around which the fminbnd function searches. The greatest error in
%the approximation would be expected when going from one pole to another,
%as the flattening of the Earth is what invalidates the approximation. The
%maximum error observed going between poles on the Earth was about 72.7m.
%The search region for the value of dist used in the fminbnd function was
%set to 0.9*dist to 1.1*dist, where dist is the distance obtained after
%scaling the distance from the zero-altitude solution.
%
%The algorithm is fastest if height=0. It is moderately slower if height is
%not zero and useHeightApprox=true, and it is very slow if height!=0 and
%useHeightApprox=false.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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

[azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd,a,f);

%If a non-zero height is given, then iterate over the direct geodetic
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
        distFun=@(distCur)distCostFunc(distCur,endCart,latLonStart,azStart,height,a,f,numSteps4Circ);
        dist=fminbnd(distFun,0.9*dist,1.1*dist);
    end
    %Only determine the ending azimuth if necessary.
    if(nargout>2)
        [~,azEnd]=directGeodeticProbGen(latLonStart,azStart,dist,height,false,a,f,numSteps4Circ);
    end
end
end

function cost=distCostFunc(distCur,endCart,latLonStart,azStart,height,a,f,numSteps4Circ)
    latLonCalc=directGeodeticProbGen(latLonStart,azStart,distCur,height,false,a,f,numSteps4Circ);
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
