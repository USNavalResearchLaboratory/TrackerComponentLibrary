function xPosVel=makeConstSpeedLevelCurvedEarthTraj(llhStart,headingStart,speed,timeOffsets,useRhumb,trajParams)
%%MAKECONSTSPEEDLEVELCURVEDEARTHTRAJ Create a non-maneuvering, constant
%   (local tangent plane) speed, constant ellipsoidal height trajectory
%   over an ellipsoidal Earth, with position and velocity sampled at
%   specified times. By default, this is a geodesic path, but a rhumb
%   trajectory (constant heading) can also be used. Note that rhumb lines
%   will get stuck at the geographic poles.
%
%INPUTS: llhStart The [latitude;longitude;ellipsoidal height] of the
%                 starting point (time=0) of the trajectory.
%    headingStart The initial heading of the trajectory in radians East of
%                 North.
%           speed The positive speed of the target. This will typically be
%                 in meters.
%     timeOffsets A length numPts list of times where the position and
%                 velocity of the target should be obtained. This will
%                 typically be in seconds.
%        useRhumb If this is true, a rhumb (constant heading) trajectory is
%                 created as opposed to a geodesic (straightest)
%                 trajectory. The default if omitted or an empty matrix is
%                 passed is false.
%      trajParams An optional structure containing values that change how
%                 the trajectory is generated. Possible field and their
%                 defaults (if omitted or empty matrices are passed) are:
%                 a The semi-major axis of the reference ellipsoid. The
%                   default is Constants.WGS84SemiMajorAxis.
%                 f The flattening factor of the reference ellipsoid. The
%                   default is Constants.WGS84Flattening.
%                 useHeightApprox If true, an approximation for targets
%                   with nonzero height is used that is much faster than a
%                   numerical integration-based technique. The default is
%                   true.
%                 numSteps4Circ If useHeightApprox is false, then numerical
%                   integration is performed and this is the number of
%                   Runge-Kutta steps that would be taken when
%                   circumnavigating the equator. The default is 2000,
%                   though 6000 is often better.
%
%OUTPUTS: xPosVel The 6XnumPts set of [position;velocity] vectors of the
%                 target at the specified times.
%
%This function obtains positions using directRhumbProbGen or
%directGeodeticProbGen and then converts the local heading to ECEF, which
%when given magnitude speed is an approximation to the instantaneous
%velocity vector.
%
%October 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(useRhumb))
    useRhumb=false; 
end

useHeightApprox=true;
numSteps4Circ=2000;
a=Constants.WGS84SemiMajorAxis;
f=Constants.WGS84Flattening;

if(nargin>5&&~isempty(trajParams))
    if(isfield(trajParams,'a'))
        a=trajParams.a;
    end
    
    if(isfield(trajParams,'f'))
        f=trajParams.f;
    end
    
    if(isfield(trajParams,'useHeightApprox'))
        useHeightApprox=trajParams.useHeightApprox;
    end
    
    if(isfield(trajParams,'numSteps4Circ'))
        numSteps4Circ=trajParams.numSteps4Circ;
    end
end

numPts=length(timeOffsets);
height=llhStart(3);
latLonStart=llhStart(1:2);

xPosVel=zeros(6,numPts);
if(useRhumb)
    %The instantaneous local ENU velocity never changes, so find it outside
    %the loop.
    vENULocal=geogHeading2ENUVec(headingStart,0)*speed;
    
    for curOffset=1:numPts
        curDist=timeOffsets(curOffset)*speed;
        latLonEnd=directRhumbProbGen(latLonStart,headingStart,curDist,height,useHeightApprox,a,f,numSteps4Circ);
        vECEF=ENUVec2ECEFVec(latLonEnd,vENULocal,a,f);
        posCart=ellips2Cart([latLonEnd;height],a,f);
        xPosVel(:,curOffset)=[posCart;vECEF];
    end
else
    %A geodesic trajectory.
    for curOffset=1:numPts
        curDist=timeOffsets(curOffset)*speed;
        [latLonEnd,azEnd]=directGeodeticProbGen(latLonStart,headingStart,curDist,height,useHeightApprox,a,f,numSteps4Circ);

        vENULocal=geogHeading2ENUVec(azEnd,0)*speed;
        vECEF=ENUVec2ECEFVec(latLonEnd,vENULocal,a,f);
        posCart=ellips2Cart([latLonEnd;height],a,f);
        xPosVel(:,curOffset)=[posCart;vECEF];
    end
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
