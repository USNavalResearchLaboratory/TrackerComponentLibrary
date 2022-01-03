function [distVal,x,xL,t]=lineOfSightMinDist2Earth(llhStart,llhEnd,a,f,epsConst)
%%LINEOFSIGHTMINDIST2EARTH Using an ellipsoidal Earth model (no tides or
%       terrain), determine the closest distance between a line of sight
%       between two points and the surface of the Earth (ignoring
%       refraction). If one is looking over a long distance and the minimum
%       height is rather low, one might expect notable refractive effects
%       to arise.
%      
%INPUTS: llhStart A 3X1 [latitude;longitude;height] vector with the angles
%           in radians, for the starting point of the line of sight
%           considered.
%    llhEnd A 3X1 [latitude;longitude;height] vector with the angles in
%           radians, for the ending point of the line of sight considered.
%         a The semi-major axis of the reference ellipsoid. If this
%           argument is omitted or an empty matrix is passed, the value in
%           Constants.WGS84SemiMajorAxis is used.
%         f The flattening factor of the reference ellipsoid. If this
%           argument is omitted or an empty matrix is passed, the value in
%           Constants.WGS84Flattening is used.
%  epsConst In the rooting formulation of the solution, is is possible that
%           some candidate solutions are invalid (they do not satisfy the
%           constraint (x-xc)'*A*(x-xc)=1. Thus, the absolute value of a
%           transformed version of abs((x-xc)'*A*(x-xc)-1) must be
%           <=epsConst to be valid. The default if omitted or an empty
%           matrix is passed is 1e-9.
%
%OUTPUTS: distVal If the line segment doesn't intersect the ground, then
%             this is the scalar nearest point from the line to the surface
%             of the ellipsoidal Earth. If the line segment intersects the
%             Earth  once or twice, then this is a 1X1 or 1X2 vector of
%             zeros. If some type of error occurred, then this and the
%             other outputs will be empty matrices.
%           x The 3X1 nearest Cartesian point on the Earth or 3X2 points
%             if there are two intersections of the segment with the Earth.
%          xL The 3X1 nearest Cartesian point on the line or 3X2 points if
%             there are two solutions.
%           t The value along the line such that xL=xStart+(xEnd-xStart)*t
%             if xStart and xEnd are llhStart and llhEnd converted to
%             Cartesian coordinates. If there are two solutions, then this
%             is a 1X2 vector.
%
%This function makes use of the findClosestPointLineEllipsoid function and,
%if an endpoint might be closest, also the nearestPointOnEllipsoid
%function.
%
%EXAMPLE:
%This demonstrates some of the various return types of this function. For
%the first point, the starting point is the closest, and is 1km from the
%ground. For the second point, the  ending point is closest and is 1km from
%the ground. For the third one, the line cuts the ground at two points, so
%so two solutions with zero distance are returned. For the fourth one,
%there is no ground intersection - the line is always above the ground and
%neither one of the endpoints is the closest point.
% llhSensor=[[0;0;1e3],[0;0.17;130e3],[0;0;1e3],[0;0;40e3]];
% llhTar=[[0;0.17;130e3],[0;0;1e3],[0;0.17;13e3],[0;0.2;40e3]];
% numLOS=size(llhSensor,2);
% typeOfHit=cell(numLOS,1);
% for k=1:numLOS
%     typeOfHit{k}=lineOfSightMinDist2Earth(llhSensor(:,k),llhTar(:,k));
% end
% %Show the solutions.
% typeOfHit
%
%November 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(epsConst))
    epsConst=1e-9;
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

xStart=ellips2Cart(llhStart,a,f);
xEnd=ellips2Cart(llhEnd,a,f);

%Special cases of one or the other points being on the ground.
if(llhStart(3)==0&&llhEnd(3)==0)
    distVal=[0,0];
    x=[xStart,xEnd];
    xL=x;
    t=[0,1];
    return;
elseif(llhStart(3)==0)
    distVal=0;
    x=xStart;
    xL=xStart;
    t=0;
    return;
elseif(llhEnd(3)==0)
    distVal=0;
    x=xEnd;
    xL=xEnd;
    t=1;
    return;
end

b=a*(1-f);%Semi-minor axis.
A=diag([1/a^2,1/a^2,1/b^2]);
xc=[0;0;0];

[x,xL,t,intersects]=findClosestPointLineEllipsoid(xc,A,xStart,xEnd,1,epsConst);

endpointIsSol=false;
if(intersects)
    %If both intersection points along the line are outside of the
    %segment, then one of the endpoints is the closest point. Otherwise, 
    numSol=length(t);
    distInSegment=(t>=0&t<=1);
    if(all(distInSegment==1))
        %If all points hit the ground.
        distVal=zeros(1,numSol);
        return;
    elseif(distInSegment(1)==1)
        %Just the first point hits the ground.
        x=x(:,1);
        xL=xL(:,1);
        t=t(1);
        distVal=0;
        return;
    elseif(numSol>1&&distInSegment(2)==1)
        %Just the second point hits the ground.
        x=x(:,2);
        xL=xL(:,2);
        t=t(2);
        distVal=0;
        return;
    end
    
    endpointIsSol=true;
end

if(isempty(t))
    %If there was some type of an error.
    distVal=[];
    return;
end

if(endpointIsSol||t<0||t>1)
    %If the closest point is outside of the segment, then 
    %If here, then all of the intersection points are outside of the line
    %segment. We must return the solution for whichever of the endpoints is
    %closest.
    x1=nearestPointOnEllipsoid(xc,A,xStart,1,epsConst);
    x2=nearestPointOnEllipsoid(xc,A,xEnd,1,epsConst);
    
    dist1=norm(x1-xStart);
    dist2=norm(x2-xEnd);
    if(dist1<dist2)
        %The starting point is closest.
        x=x1;
        xL=xStart;
        t=0;
        distVal=dist1;
    elseif(dist2<dist1)
        %The ending point is closest.
        x=x2;
        xL=xEnd;
        t=1;
        distVal=dist2;
    end
    return;
end

distVal=norm(x-xL);
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
