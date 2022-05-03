function J=calcAzEquidistantCartJacob(latLonPt,latLonRef,rE)
%%CALCAZEQUIDISTANTCARTJACOB Find the gradient of an azimuthal equidistance
%           projection of points on a spherical Earth with respect to the 
%           Cartesian [x;y;z] components.
%
%INPUTS: latLonPt A 2X1 [latitude;longitude] point in radians at which
%                 the gradient should be evaluated. If a 3X1 vector is
%                 passed, then it is assumed that it is a height abov the
%                 reference sphere and a third row will be present in the
%                 output.
%       latLonRef A 2X1 reference point about which the projection is
%                 taken. Note that latLonRef shouldn't coincide with
%                 latLonPt, because a singularity exists at that point.
%              rE The radius of the reference sphere. If this is omitted
%                 or an empty matrix is passed, then
%                 osculatingSpher4LatLon(latLonRef) is used.
%
%OUTPUTS: J The 2X3 or 3X3 Jacobian matrix with derivatives with respect to
%           Cartesian position components. Each row is a component of
%           azimuthal equidistant x and y (and height iflatLonPt is 3X1) in
%           that order with derivaties of x, y, and z given across the
%           columns in that order.
%
%Starting with the conversion in Chapter 25 of [1], substitute in
%expressions for x, y and x in spherical coordinates and then evaluate
%derivatives with respect to them.
%
%EXAMPLE:
%The Jacobian as obtained by this function is compared to one via numeric
%differentiation. The relative error implies over 10 digits of precision,
%which is consistent with what one might expect with finite precision
%limitation.
% latLonPt=deg2rad([20.631756;-155.551818]);
% latLonRef=deg2rad([20.906029;-157.078918]);
% rE=osculatingSpher4LatLon(latLonRef);
% f=@(pt)Cart2AzEquidistantProj(pt,latLonRef,rE,0);
% xyzPt=ellips2Cart(latLonPt,rE,0);
% JNumDiff=numDiff(xyzPt,f,3,[],1e-4*norm(xyzPt));
% JNumDiff=JNumDiff(1:2,:);
% J=calcAzEquidistantCartJacob(latLonPt,latLonRef,rE);
% RelErr=max(abs((J(:)-JNumDiff(:))./JNumDiff(:)))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(rE))
    rE=osculatingSpher4LatLon(latLonRef);
end

xyzPt=ellips2Cart(latLonPt,rE,0);
x=xyzPt(1);
y=xyzPt(2);
z=xyzPt(3);

p=latLonPt(1);
sinP=sin(p);
cosP=cos(p);
l=latLonPt(2);

p0=latLonRef(1);
sinP0=sin(p0);
cosP0=cos(p0);
l0=latLonRef(2);
coslDiff=cos(l0-l);
sinlDiff=sin(l0-l);
cv=acos(sinP0*sinP+cosP0*cosP*coslDiff);
cosCV=cos(cv);
cscCV=csc(cv);
cotCV=cot(cv);

r2=x^2+y^2+z^2;
rxy2=x^2+y^2;
rxy=sqrt(rxy2);

dpdx=-((x*z)/(rxy*r2));
dpdy=-((y*z)/(rxy*r2));
dpdz=rxy/r2;
dldx=-(y/rxy2);
dldy=x/rxy2;
dldz=0;

dcdx=-(cosP0*cosP*sinlDiff*dldx+(cosP*sinP0-cosP0*sinP*coslDiff)*dpdx)/sqrt(1-cosCV^2);
dcdy=-(cosP0*cosP*sinlDiff*dldy+(cosP*sinP0-cosP0*sinP*coslDiff)*dpdy)/sqrt(1-cosCV^2);
dcdz=-(cosP0*cosP*sinlDiff*dldz+(cosP*sinP0-cosP0*sinP*coslDiff)*dpdz)/sqrt(1-cosCV^2);

dxAdx=rE*cscCV*(cosP*(-1+cv*cotCV)*sinlDiff*dcdx+cv*(coslDiff*cosP*dldx+sinlDiff*sinP*dpdx));
dxAdy=rE*cscCV*(cosP*(-1+cv*cotCV)*sinlDiff*dcdy+cv*(coslDiff*cosP*dldy+sinlDiff*sinP*dpdy));
dxAdz=rE*cscCV*(cosP*(-1+cv*cotCV)*sinlDiff*dcdz+cv*(coslDiff*cosP*dldz+sinlDiff*sinP*dpdz));

dyAdx=rE*cscCV*((-1+cv*cotCV)*(coslDiff*cosP*sinP0-cosP0*sinP)*dcdx+cv*(-cosP*sinP0*sinlDiff*dldx+(cosP0*cosP+coslDiff*sinP0*sinP)*dpdx));
dyAdy=rE*cscCV*((-1+cv*cotCV)*(coslDiff*cosP*sinP0-cosP0*sinP)*dcdy+cv*(-cosP*sinP0*sinlDiff*dldy+(cosP0*cosP+coslDiff*sinP0*sinP)*dpdy));
dyAdz=rE*cscCV*((-1+cv*cotCV)*(coslDiff*cosP*sinP0-cosP0*sinP)*dcdz+cv*(-cosP*sinP0*sinlDiff*dldz+(cosP0*cosP+coslDiff*sinP0*sinP)*dpdz));

if(length(latLonPt)==3)
    r=sqrt(r2);
    dhdx=x/r;
    dhdy=y/r;
    dhdz=z/r;

    J=[dxAdx,dxAdy,dxAdz;
       dyAdx,dyAdy,dyAdz
        dhdx, dhdy,dhdz];
else
    J=[dxAdx,dxAdy,dxAdz;
       dyAdx,dyAdy,dyAdz];
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
