function latLonPoint=rhumbIntersect(az1,az2,latLon1,latLon2,a,f)
%%RHUMBINTERSECT Given two points and headings, determine where the
%                constant-heading courses (rhumb lines) intersect. If the
%                course circles the Earth (in longitude) more than once
%                before an intersection, then an empty matrix is returned.
%
%INPUTS: az1 The direction of the target in radians East of North as
%            measured by the first sensor.
%        az2 The direction of the target in radians East of North as
%            measured by the second sensor.
%    latLon1 The 2X1 location of the first sensor given in latitude and
%            longitude in radians in the format [latitude;longitude].
%            This should not be the North or South pole.
%    latLon2 The 2X1 location of the second sensor given in latitude and
%            longitude in radians in the format [latitude;longitude].
%            This should not be the North or South pole.
%          a The semi-major axis of the reference ellipsoid (in meters). If
%            this argument is omitted or an empty matrix is passed, the
%            value in Constants.WGS84SemiMajorAxis is used.
%          f The flattening factor of the reference ellipsoid. If this
%            argument is omitted or an empty matrix is passed, the value in
%            Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonPoint The 2X1 [latitude;longitude] location in radians of
%                     the closest intersection that is not more than one
%                     revolution around the Earth. If no such intersection
%                     exists, then an empty matrix is returned.
%
%First, we convert the points into Mercator coordinates to get (x1,y1) and 
%(x2,y2). Rhumb lines are stright lines in mercator coordinates and az1 and
%az2 are preserved. However, we want the lines to intersect in the
%directions of the headings (not backwards). Also, the Mercator projection
%is a projection onto a cylinder. The x coordinates that differ by
%Deltax=2*pi*a alias back.
%Parametric equations for the first line are:
%x=r1*sin(az1)+x1+k1*Deltax
%y=r1*cos(az1)+y1;
%and parameteric equations for the second line are:
%x=r2*sin(az2)+x2+k2*Deltax
%y=r2*cos(az2)+y2;
%where r1 and r2 are the parameteric parameters and k1 and k2 are integers
%indicating the aliasing intervals. Intersections are where both sets of
%solutions are equal. That leads to the following expressions for r1 and
%r2:
%r1=(C2*(x2-x1)+S2*(y1-y2)+(k2-k1)*C2*Deltax)/(S1*C2-S2*C1)
%r2=(C1*(x2-x1)+S1*(y1-y2)+(k2-k1)*C1*Deltax)/(S1*C2-S2*C1)
%where C1=cos(az1), S1=sin(az1), C2=cos(az2), and S2=sin(az2). It can be
%seen that the r values depend only on the difference between k1 and k2,
%not on their actual values. We arbitrarily choose k1=0. This means that we
%have to scan k2 for solutions such that both r values are positive.
%
%This function only consider k2=-1, 0,and 1, because it is assumed that
%trajectories spiraling more than once around the world are not useful. Of
%the values of k2 that produce solutions where both r1 and r2 are positive,
%the one that minimizes r1+r2 is chosen, which is assumed to be the first
%intersection. If abs(x-x1)>Deltax or abs(x-x2)>Deltax for the solution
%found, then an empty matrix is returned, because that represents spiraling
%around the Earth.
%
%EXAMPLE 1:
%Here, we create rhumb lines with a known point of intersection, then run
%the algorithm and verify that the intersection point agrees with what we
%found. The lines and the point are plotted.
% N=100;
% latLonA=[34.685169;139.443632]*(pi/180);
% latLonC=[-33.8617;151.2117]*(pi/180);
% latLonX=[37.7917;-122.4633]*(pi/180);
% 
% %Get the azimuths.
% [azAX,distAX]=indirectRhumbProblem(latLonA,latLonX);
% [azCX,distCX]=indirectRhumbProblem(latLonC,latLonX);
% 
% dist=linspace(0,distAX,N);
% WPAX=directRhumbProblem(latLonA,azAX,dist);
% dist=linspace(0,distCX,N);
% WPCX=directRhumbProblem(latLonC,azCX,dist);
% 
% %Convert the waypoints to Cartesian to plot.
% WPAX=ellips2Cart([WPAX;zeros(1,N)]);
% WPCX=ellips2Cart([WPCX;zeros(1,N)]);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid([]);
% plot3(WPAX(1,:),WPAX(2,:),WPAX(3,:),'-r','linewidth',4)
% plot3(WPCX(1,:),WPCX(2,:),WPCX(3,:),'-b','linewidth',4)
% view(-90,0)
% 
% latLonPointX=rhumbIntersect(azAX,azCX,latLonA,latLonC);
% intersectPt=ellips2Cart([latLonPointX;0]);
% scatter3(intersectPt(1),intersectPt(2),intersectPt(3),200,'c','filled');
% max(abs(latLonPointX-latLonX))
%One will see that the point found agrees with the original point within
%finite precision limits.
%
%EXAMPLE 2:
%This is the same as example 1 without the plots, except one of the
%starting points is placed across the -pi/pi border in longitude, so the
%algorithm has to use a value of k2 that is not 0. Again, there is
%agreement within finite precision limits.
% latLonA=[34.685169;139.443632]*(pi/180);
% latLonC=[-33.8617;-151.2117]*(pi/180);
% latLonX=[37.7917;-122.4633]*(pi/180);
% 
% %Get the azimuths.
% [azAX]=indirectRhumbProblem(latLonA,latLonX);
% [azCX]=indirectRhumbProblem(latLonC,latLonX);
% 
% latLonPointX=rhumbIntersect(azAX,azCX,latLonA,latLonC);
% intersectPt=ellips2Cart([latLonPointX;0]);
% scatter3(intersectPt(1),intersectPt(2),intersectPt(3),200,'c','filled');
% max(abs(latLonPointX-latLonX))
%
%EXAMPLE 3:
%This is an example of trajectories that do not intersect within one
%revolution around the Earth. Thus, latLonPoint is an empty matrix.
% latLonA=[34.685169;139.443632]*(pi/180);
% latLonC=[-33.8617;-151.2117]*(pi/180);
% 
% azAX=0;
% azCX=pi/2;
% latLonPoint=rhumbIntersect(azAX,azCX,latLonA,latLonC)
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%Convert the points to Mercator coordinates.
xy1=ellips2Mercator(latLon1,a,f);
xy2=ellips2Mercator(latLon2,a,f);

x1=xy1(1);
y1=xy1(2);
x2=xy2(1);
y2=xy2(2);

S1=sin(az1);
C1=cos(az1);
S2=sin(az2);
C2=cos(az2);

term11=C2*(x2-x1)+S2*(y1-y2);
term12=C1*(x2-x1)+S1*(y1-y2);
denom=S1*C2-S2*C1;

%The k=0 case.
r1=(term11)/denom;
r2=(term12)/denom;

Deltax=a*2*pi;
k2=[-1;0;1];
    
r1=(term11+k2*C2*Deltax)/denom;
r2=(term12+k2*C1*Deltax)/denom;

sel=r1>0&r2>0;
r1=r1(sel);
r2=r2(sel);

if(isempty(r1))
    %There are no valid points, so we have to return an empty matrix.
    latLonPoint=[];
   return 
end

[~,minIdx]=min(r1+r2);
r1=r1(minIdx);
r2=r2(minIdx);

if(r1>Deltax||r2>Deltax)
    %More than one time around the world. 
    latLonPoint=[];
    return
end

%r1 and r2 are valid. Get the x and y points and convert back to Cartesian.
x=r1*S1+x1;
y=r1*C1+y1;

latLonPoint=Mercator2Ellipse([x;y],a,f);

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
