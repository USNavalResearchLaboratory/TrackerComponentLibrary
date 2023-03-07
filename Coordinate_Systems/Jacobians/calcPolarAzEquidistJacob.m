function J=calcPolarAzEquidistJacob(xCartLoc,latLonRef,rE)
%%CALCPOLARAZEQUIDISTJACOB Find the gradient of a point in polar azimuthal
%         equidistance coordinates with respect to the Cartesian [x;y;z]
%         components. This function assumes a spherical Earth model.
%
%INPUTS: xCartLoc The 3X1 Cartesian point at which the gradient should be
%                 evaluated.
%       latLonRef A 2X1 reference point about which the projection is
%                 taken. Note that latLonRef shouldn't coincide with
%                 latLonPt, because a singularity exists at that point.
%              rE The radius of the reference sphere. If this is omitted
%                 or an empty matrix is passed, then
%                 osculatingSpher4LatLon(latLonRef) is used.
%
%OUTPUTS: J The Jacobian matrix with derivatives with respect to Cartesian
%           position components. Each row is a component polar azimuthal
%           equidistant range, angle and height in that order with
%           derivatives of x, y, and z given across columns in that order.
%
%EXAMPLE:
%Here, the gradient obtained by this function is compared to that obtained
%from numeric diferentiation. The relative error between the two solutions
%implies about 9 digits of agreement
% latLonRef=deg2rad([20.906029;-157.078918]);
% rE=osculatingSpher4LatLon(latLonRef);
% llhPt=[deg2rad([21.295516;-158.002386]);10e3];
% xCartLoc=ellips2Cart(llhPt,rE,0);
% J=calcPolarAzEquidistJacob(xCartLoc,latLonRef,rE);
% f=@(x)Cart2PolarAzEquidistProj(x,latLonRef,rE,0);
% JNumDiff=numDiff(xCartLoc,f,3);
% RelErr=max(max(abs((J-JNumDiff)./JNumDiff)))
%
%January 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(rE))
    rE=osculatingSpher4LatLon(latLonRef);
end

x=xCartLoc(1);
y=xCartLoc(2);
z=xCartLoc(3);

phi0=latLonRef(1);
lambda0=latLonRef(2);

cosPhi0=cos(phi0);
sinPhi0=sin(phi0);
cosLambda0=cos(lambda0);
sinLambda0=sin(lambda0);

r2=dot(xCartLoc,xCartLoc);
r=sqrt(r2);

numVal=(x*cosLambda0*cosPhi0+y*cosPhi0*sinLambda0+z*sinPhi0)^2;
d1=r^3*sqrt(1-numVal/r2);
d2=(y*cosLambda0-x*sinLambda0)^2+(z*cosPhi0-sinPhi0*(x*cosLambda0+y*sinLambda0))^2;
d3=r;

J11=(rE/d1)*(x*y*cosPhi0*sinLambda0+x*z*sinPhi0-(y^2+z^2)*cosLambda0*cosPhi0);
J21=(rE/d1)*(x*y*cosLambda0*cosPhi0-(x^2+z^2)*cosPhi0*sinLambda0+y*z*sinPhi0);
J31=(rE/d1)*(z*cosPhi0*(x*cosLambda0+y*sinLambda0)-(x^2+y^2)*sinPhi0);
J12=(y*sinPhi0-z*cosPhi0*sinLambda0)/d2;
J22=(z*cosLambda0*cosPhi0-x*sinPhi0)/d2;
J32=cosPhi0*(x*sinLambda0-y*cosLambda0)/d2;
J13=x/d3;
J23=y/d3;
J33=z/d3;

J=[J11,J21,J31;
   J12,J22,J32;
   J13,J23,J33];

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
