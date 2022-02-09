function typeOfHit=lineOfSightHitsEarth(llhStart,llhEnd,a,f)
%%LINEOFSIGHTHITSEARTH Using an ellipsoidal Earth model (no tides or
%       terrain), determine whether a particular line of sight is
%       obstructed by the Earth. That is, whether a line segment from the
%       sensor to the target intercepts the ellipsoidal Earth
%       approximation. This does not take atmospheric refraction into
%       account.
%
%INPUTS: llhStart A 3XnumLOS set of [latitude;longitude;height] values
%               with the angles in radians, for the starting points of the
%               lines of sight considered.
%        llhEnd A 3XnumLOS set of [latitude;longitude;height] values with
%               the angles in radians, for each of the ending points of the
%               lines of sight considered.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84Flattening is used.
%
%OUTPUTS: typeOfHit A numLOSX1 set of values indicating the type of
%               intersection of each line segment and the Earth. Possible
%               values are:
%               0 The line of sight never hits the ground.
%               1 The starting and ending points are not on the ground, but
%                 the line of sight hits the ground.
%               2 The line of sight starts on the ground and does not
%                 otherwise touch the ground.
%               3 The line of sight starts on the ground and does otherwise
%                 touch the ground.
%               4 The line of sight starts in the air and the target is on
%                 the ground, but it does not otherwise touch the ground.
%               5 The line of sight starts in the air and the target is on
%                 the ground, and it does otherwise touch the ground.
%               6 The line of sight starts and ends on the ground and thus
%                 is always underground.
%               7 The starting and/or ending point is underground.
%
%The implementation of the function is via a fairly straightforward use of
%the function findEllipsLineIntersect given that the reference ellipsoid is
%axis-aligned with the semi-major axis along the z-axis.
%
%EXAMPLE:
%In this example, seven line segments are passed satisfying all of the
%different intersection criteria. The returned typeOfHit vector should be
%[0;1;2;3;4;5;6;7].
% llhSensor=[[0;0;1e3],[0;0;1e3],[0;0;0],[0;0;0],[0;0;100e3],[0;0;1e3],[0;0;0]];
% llhTar=[[0;0.17;130e3],[0;0.17;13e3],[0;0.17;130e3],[0;0.17;13e3],[0;0.17;0],[0;0.17;0],[0;0.17;0]];
% typeOfHit=lineOfSightHitsEarth(llhSensor,llhTar)
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

numPts=size(llhStart,2);
typeOfHit=zeros(numPts,1);

b=a*(1-f);%Semi-minor axis.
A=diag([1/a^2,1/a^2,1/b^2]);
xc=[0;0;0];

for k=1:numPts
    typeOfHit(k)=lineOfSightHitsEarthOnePt(llhStart(:,k),llhEnd(:,k),A,xc,a,f);
end
end

function typeOfHit=lineOfSightHitsEarthOnePt(llhSensor,llhTar,A,xc,a,f)

if(llhSensor(3)<0||llhTar(3)<0)
    %The starting and/or ending point is underground.
    typeOfHit=7;
    return;
end


startIsZero=llhSensor(3)==0;
endIsZero=llhTar(3)==0;

if(startIsZero&&endIsZero)
    typeOfHit=6;
    return;
end

xRx=ellips2Cart(llhSensor,a,f);
xTar=ellips2Cart(llhTar,a,f);
[xInt,t]=findEllipsLineIntersect(xc,A,xRx,xTar,1);

if(startIsZero&&~endIsZero)
    %The start is on the ground.
    if(length(t)<2)
        %The only point touching the ground is the starting point. 
        typeOfHit=2;
    else
        %If there are two points, then get rid of the point closest to the
        %starting point and see whether the other point had a t value
        %between 0 and 1.
        
        if(norm(xInt(:,1)-xRx)<norm(xInt(:,2)-xRx))
            idx=2;
        else
            idx=1;
        end
        
        if(t(idx)>0&&t(idx)<1)
            %If the intersection point is between the starting and ending
            %points of the segment (i.e. it intersects between them).
            typeOfHit=3;
        else
            %The line of sight only hits at the starting point.
            typeOfHit=2;
        end
    end
elseif(~startIsZero&&endIsZero)
    %The end is on the ground.
    
    if(length(t)<2)
        %The only point touching the ground is the ending point. 
        typeOfHit=4;
    else
        %If there are two points, then get rid of the point closest to the
        %starting point and see whether the other point had a t value
        %between 0 and 1.
        
        if(norm(xInt(:,1)-xTar)<norm(xInt(:,2)-xTar))
            idx=2;
        else
            idx=1;
        end
        
        if(t(idx)>0&&t(idx)<1)
            %If the intersection point is between the starting and ending
            %points of the segment (i.e. it intersects between them).
            typeOfHit=5;
        else
            %The line of sight only hits at the ending point.
            typeOfHit=4;
        end
    end
else
    %Neither the starting nor the ending points are on the ground. 
    if(isempty(t)||(t(1)<0||t(1)>1)&&(t(2)<0||t(2)>1))
        %The line of sight never hits the ground.
        typeOfHit=0;
    else
        %The line of sight hits the ground somewhere.
        typeOfHit=1;
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
