function latLonPoint=geodesicIntersect(az1,az2,latLon1,latLon2,convgCrit,a,f,AbsTol,maxIter)
%%GEODESICINTERSECT Earthbound sensors 1 and 2 located at latLon1 and
%           latLon2 and measure the azimuthal directions of a target given
%           in radians East of North, az1 and az2. This function determines
%           the latitude and longitude of the target on an ellipsoidal
%           Earth.
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
%  convgCrit This determines the convergence criterion to be used . This
%            also determines the distance measure used when choosing between 
%            solutions in the iteration. Possible values are:
%            0 (The default if omitted or an empty matrix is passed) The
%              norm of A difference in latitude, longitude (wrapped to -pi
%              to pi in each element) is used.
%            1 Cartesian distances between the points are used.
%          a The semi-major axis of the reference ellipsoid. If this
%            argument is omitted or an empty matrix is passed, the value in
%            Constants.WGS84SemiMajorAxis is used.
%          f The flattening factor of the reference ellipsoid. If this
%            argument is omitted or an empty matrix is passed, the value in
%            Constants.WGS84Flattening is used.
%     AbsTol The absolute tolerance in the minimum distance specified by
%            convgCrit to use for determining algorithmic convergence. the
%            default if omitted or an empty matrix is passed is 1e-10 if
%            convgCrit=0 and 7e-4 if convgCrit=1.
%    maxIter The maxium number of iterations to perform. The default if
%            omitted or an empty matrix is passed is 50. maxIter>=1.
%
%OUTPUTS: latLonPoint The 2X1 [latitude;longitude] location in radians of
%                     the target. If multiple solutions exist, this will
%                     typically be the one closest to latLon1.
%
%The algorithm of [1] is used. In the spherical approximation step, we use
%the mean radius of the Earth as implied by a and f. This is obtained using
%a formula from [2].
%
%EXAMPLE:
%This is an example of a bearing bearing taken near Japan and one being
%taken in the pacific being used to localize something in California. We
%shall plot the geodesic triangle on the ellipsoid. We shall estimate the
%target location using a spherical approximation with the functon
%greatCircleIntersect as well as using the ellipsoidal approximation with
%this function (on which the measurement model is based). Cartesian
%location estimation errors shall be displayed. It is seen that the
%spherical approximation produces an error of over 9 kilometers, whereas
%the geodesic approximation has a sub-meter approximation error.
% N=100;
% latLonA=[34.685169;139.443632]*(pi/180);
% latLonC=[22;-151.2117]*(pi/180);
% latLonX=[37.7917;-122.4633]*(pi/180);
% xCart=ellips2Cart([latLonX;0]);
% 
% %Get the azimuth
% [azAX,distAX]=indirectGeodeticProb(latLonA,latLonX);
% [azCX,distCX]=indirectGeodeticProb(latLonC,latLonX);
% [azAC,distAC]=indirectGeodeticProb(latLonA,latLonC);
% 
% dist=linspace(0,distAX,N);
% WPAX=directGeodeticProb(latLonA,azAX,dist);
% dist=linspace(0,distCX,N);
% WPCX=directGeodeticProb(latLonC,azCX,dist);
% dist=linspace(0,distAC,N);
% WPAC=directGeodeticProb(latLonA,azAC,dist);
% 
% %Convert the waypoints to Cartesian to plot. Plot slightly above the
% %surface so that the line shows up better.
% WPAX=ellips2Cart([WPAX;0.02*ones(1,N)]);
% WPAC=ellips2Cart([WPAC;0.02*ones(1,N)]);
% WPCX=ellips2Cart([WPCX;0.02*ones(1,N)]);
% 
% figure(1)
% clf
% hold on
% plotMapOnEllipsoid();
% plot3(WPAX(1,:),WPAX(2,:),WPAX(3,:),'-r','linewidth',4)
% plot3(WPAC(1,:),WPAC(2,:),WPAC(3,:),'-g','linewidth',4)
% plot3(WPCX(1,:),WPCX(2,:),WPCX(3,:),'-b','linewidth',4)
% 
% latLonPointX=greatCircleIntersect(azAX,azCX,latLonA,latLonC);
% xEstCart=ellips2Cart([latLonPointX;0]);
% norm(xEstCart-xCart)%Spherical approximation error.
% 
% latLonPointXG=geodesicIntersect(azAX,azCX,latLonA,latLonC);
% xEstCartG=ellips2Cart([latLonPointXG;0]);
% norm(xEstCartG-xCart)%Ellipsoidal approximation error.
% 
% scatter3(xEstCart(1),xEstCart(2),xEstCart(3),100,'r','filled')
% scatter3(xEstCartG(1),xEstCartG(2),xEstCartG(3),100,'c','filled')
% view(-80,20)
%
%REFERENCES:
%[1] S. Baselga, J. C. Martinez-Llario, "Intersection and point-to-line
%    solutions for geodesics on the ellipsoid," Studi Geophysica et
%    Geodaetica, vol. 62, no. 3, pp. 353-363, Jul. 2018.
%[2] H. Mortiz, "Geodetic Reference System 1980," Bulletin Geodesique,
%    vol. 54, no. 3, pp. 395-405, Sep. 1980. Given with corrections at
%    https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(convgCrit))
    convgCrit=0;
end

if(nargin<8||isempty(AbsTol))
    if(convgCrit==0)
        AbsTol=1e-10;
    else
        AbsTol=7e-4;
    end
end

if(nargin<9||isempty(maxIter))
    maxIter=50;
end

if(nargin<7||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<6||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%Semi-minor axis of the ellipsoid.
b=a*(1-f);

%For the radius to use with the spherical approximation, use the mean
%radius of the specified reference ellipsoid.
%The formula to derive the mean radius is given in [2].
rMean=(2*a+b)/3;

[~,dist1,dist2]=greatCircleIntersect(az1,az2,latLon1,latLon2,rMean);

for curIter=1:maxIter
    [latLon1,az1]=directGeodeticProb(latLon1,az1,dist1,a,f);
    [latLon2,az2]=directGeodeticProb(latLon2,az2,dist2,a,f);
    
    [latLonPoint,dist1,dist2]=greatCircleIntersect(az1,az2,latLon1,latLon2,rMean,true);
    if(convgCrit==0)
    	diff1=norm(wrapRange(latLonPoint(:,1)-latLon1,-pi,pi));
        diff2=norm(wrapRange(latLonPoint(:,2)-latLon2,-pi,pi));
    else
        sol1CartLoc=ellips2Cart([latLon1;0],a,f);
        diff1=norm(ellips2Cart([latLonPoint(:,1);0],a,f)-sol1CartLoc);
        diff2=norm(ellips2Cart([latLonPoint(:,2);0],a,f)-sol1CartLoc);
    end

    if(norm(diff1)<norm(diff2))
        dist1=dist1(1);
        dist2=dist2(1);
    else
        dist1=dist1(2);
        dist2=dist2(2);
    end
    
    if(min(diff1,diff2)<=AbsTol)
        break;
    end
end

latLonPoint=directGeodeticProb(latLon1,az1,dist1,a,f);
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
