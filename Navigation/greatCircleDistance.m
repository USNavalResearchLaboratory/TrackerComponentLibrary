function distVal=greatCircleDistance(latLon1,latLon2,rE,formula2Use)
%%GREATCIRCLEDISTANCE Given the spherical latitudes and longitudes of two
%       points on a reference sphere, determine the distance between them
%       across the surface of the sphere.
%
%INPUTS: latLon1 The 2XN set of N first points given in latitude and
%                longitude in radians in the format [latitude;longitude]
%                (on a reference sphere, latitude is spherical elevation;
%                longitude is spherical azimuth). If all these points are
%                the same but latLon2 varies, then a single 2X1 vector can
%                be passed.
%        latLon2 The 2XN second points given in latitude and longitude in
%                radians in the format [latitude;longitude]. If all these
%                points are the same but latLon1 varies, then a single 2X1
%                vector can be passed.
%             rE The assumed radius of the spherical Earth model. If
%                omitted or an empty matrix is passed, the default of
%                rE=osculatingSpher4LatLon(latLon1) is used.
%    formula2Use This optional value specified the formula to use for the
%                conversion. Possible values are:
%                0 (The default if omitted or an empty matrix is passed)
%                  Convert the points into unit vectors and use the Earth's
%                  radius times the angular distance between them (via
%                  angBetweenVecs) to get the distance across the Earth.
%                  This is more accurate than the other method.
%                1 Use the cosine formula that is implicit in Equations
%                  13-15 in [1].
%
%OUTPUTS: distVal The scalar distance between the points.
%   
%EXAMPLE:
%This example just shows that when the points are far apart, the two
%formula field essentially the same distance (relative error on the order
%of eps). However, when the points are very close, the relative difference
%between the distance formula increases significantly.
% latLon1=deg2rad([19.7241;-155.0868]);
% latLon2=deg2rad([48.8566;2.3522]);
% latLon3=deg2rad([19.7240;-155.0867]);
% d0=greatCircleDistance(latLon1,[latLon2,latLon3],[],0);
% d1=greatCircleDistance(latLon1,[latLon2,latLon3],[],1);
% RelDiff12=(d1-d0)./d0
%
%REFERENCES:
%[1] C.-L. Chen, P.-F. Liu, and W.-T. Gong, "A simple approach to great
%    circle sailing: The COFI method," The Journal of Navigation, vol. 67,
%    no. 3, pp. 403-418, May 2014.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(formula2Use))
    formula2Use=0;
end

if(nargin<3||isempty(rE))
    rE=osculatingSpher4LatLon(latLon1);
end

%Latitude at the start.
phi1=latLon1(1,:);
%Latitude at the destination.
phi2=latLon2(1,:);

lambda1=latLon1(2,:);
lambda2=latLon2(2,:);
cosPhi1=cos(phi1);
cosPhi2=cos(phi2);

if(formula2Use==0)
    %Use the distance between vectors method.
    
    %The Cartesian starting point (unit sphere).
    F=[cos(lambda1).*cosPhi1;
       sin(lambda1).*cosPhi1;
       sin(phi1)];

    %The Cartesian ending point (unit sphere).
    T=[cos(lambda2).*cosPhi2;
       sin(lambda2).*cosPhi2;
       sin(phi2)];

    distVal=rE*angBetweenVecs(F,T);
else
    %Use the cosine formula.
    angArg=sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(lambda2-lambda1);
    %Avoid invalid values due to finite precision limitiations.
    angArg=min(1,max(-1,angArg));
    distVal=rE.*acos(angArg);
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
