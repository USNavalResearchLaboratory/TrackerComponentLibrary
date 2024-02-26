function azElDiff=distanceBetweenSpherDirs(azEl0,azEl1,systemType,offsetMinType)
%%DISTBETWEENSPHERDIRS Given two [azimuth;elelvation] locations on the unit
%       sphere, find the offset azElDiff such that the position represented
%       by azEl0+azElDiff equals the position represented by azEl1.
%       Multiple solutions are available including one minimizing the
%       change in the azimuth of azEl0, one minimizing the elevation offset
%       in azEl0 and one minimizing norm(azElDiff).
%
%INPUTS: azEl0, azEl1 These are two 2X1 [azimuth;elevation] points in
%           radians. It is assumed that azimuth ranges from -pi to pi and
%           that elevation is from -pi/2 to pi/2 if systemType~=2 and if
%           systemType==2, then it ranges from 0 to pi.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are:
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
% offsetMinType This specifies which type of solution to obtain. Possible
%             values are:
%             0 (The default if omitted or an empty matrix is passed)
%               Minimizes the magnitude of the azimuth offset.
%             1 Minimize the elevation offset
%             2 Minimize norm(azElDiff).
%
%OUTPUT: azElDiff The [2;1] offset that can be added to azEl0 to get a
%                  point in the same location as azEl1;
%
%EXAMPLE 1:
%This computes the offsets for two random points and then find the
%absolute difference in the position converted into Cartesian coordinates
%(r=1). This is done of 1e4 Monte Carlo runs and the maximum offset is what
%one might expect given finite precision limits.
% numMCRuns=1e4;
% systemType=0;
% offsetMinType=0;
% 
% maxAbsDiff=0;
% for k=1:numMCRuns
%     %Random points
%     if(systemType==2)
%         z0=[-pi+2*pi*rand(1);pi*rand(1)];
%         z1=[-pi+2*pi*rand(1);pi*rand(1)];
%     else
%         z0=[-pi+2*pi*rand(1);-pi/2+pi*rand(1)];
%         z1=[-pi+2*pi*rand(1);-pi/2+pi*rand(1)];
%     end
% 
%     deltaVal=distanceBetweenSpherDirs(z0,z1,systemType,offsetMinType);
%     maxAbsDiff=max(maxAbsDiff,norm(spher2Cart(z0+deltaVal,systemType)-spher2Cart(z1,systemType)));
% end
% maxAbsDiff
%
%EXAMPLE 2
%Here is a simple example using points initially specified in degrees and
%the offsets are converted into degrees. Here, we see that the solutions
%minimizing the azimuthal offset and minimizing the elevation offset
%magnitudes differ significantly.
% z0=deg2rad([45;89]);
% z1=deg2rad([-137;88]);
% systemType=0;
% offsetMinType=0;
% deltaVal=rad2deg(distanceBetweenSpherDirs(z0,z1,systemType,offsetMinType))
% offsetMinType=1;
% deltaVal=rad2deg(distanceBetweenSpherDirs(z0,z1,systemType,offsetMinType))
%
%March 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(offsetMinType))
    offsetMinType=0;%Minimize the azimuthal offset.
end

if(nargin<3||isempty(systemType))
    systemType=0;
end

if(systemType==2)
    el0=acos(cos(azEl0(2)));
    el1=acos(cos(azEl1(2)));

    DeltaEl=[-el0-el1;
             -el0+el1];
else
    el0=asin(sin(azEl0(2)));
    el1=asin(sin(azEl1(2)));
    DeltaEl=[pi-el0-el1;
            -el0+el1];
end

az0=azEl0(1);
az1=azEl1(1);

%Make the offsets as small as possible (check whether adding/subtracting
%2*pi leads to a smaller solution):
v=[DeltaEl+2*pi,DeltaEl,DeltaEl-2*pi];
[~,idx1]=min(abs(v(1,:)));
[~,idx2]=min(abs(v(2,:)));
DeltaEl(1)=v(1,idx1);
DeltaEl(2)=v(2,idx2);

%Having offset the elevation, we could have changed the azimuth of the
%original measurement.
azOffset=zeros(2,1);
for k=1:2
    elCur=el0+DeltaEl(k);
    azCur=getSpherAz(spher2Cart([az0;elCur],systemType),systemType);
    azOffset(k)=wrapRange(az1-azCur,-pi,pi);
end

switch(offsetMinType)
    case 0%Minimize the azimuthal offset.
        if(abs(azOffset(1))<=abs(azOffset(2)))
            azElDiff=[azOffset(1);DeltaEl(1)];
        else
            azElDiff=[azOffset(2);DeltaEl(2)];
        end
    case 1%Minimize the elevation offset.
        if(abs(DeltaEl(1))<=abs(DeltaEl(2)))
            azElDiff=[azOffset(1);DeltaEl(1)];
        else
            azElDiff=[azOffset(2);DeltaEl(2)];
        end
    case 2%Minimize the total offset
        v1=[azOffset(1);DeltaEl(1)];
        v2=[azOffset(2);DeltaEl(2)];

        if(norm(v1)<=norm(v2))
            azElDiff=v1;
        else
            azElDiff=v2;
        end
    otherwise
        error('Unknown offset type specified.')
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
