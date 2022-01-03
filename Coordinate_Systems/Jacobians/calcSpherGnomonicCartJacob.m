function J=calcSpherGnomonicCartJacob(xyzPt,latLonRef)
%%CALCSPHERGNOMONICCARTJACOB Calculate the Jacobian of the transformation
%           from Cartesian coordinates on the surface of a unit sphere to
%           gnomonic coordinates. This is essentially the Jacobian of the
%           transformation given in the function Cart2RSpherGnomonicProj.
%
%INPUTS: xyzPt A 3X1 point on the reference sphere in [x;y;z] Cartesian
%              coordinates where the Jacobian should be evaluated.
%    latLonRef The 2X1 [latitude;longitude] reference point for the
%              gnomonic projection. If the reference point is the same as
%              xyzPt, then this input can be omitted.
%
%OUTPUTS: J The Jacobian matrix with derivatives with respect to Cartesian
%           position components. Each row is a component of r, Gnomonic x
%           and Gnomonic y in that order, with derivaties of x, y, and z
%           given across the columns in that order.
%
%Starting with the conversion in Chapter 22 of [1], substitute in
%expressions for x, y and x in spherical coordinates and then evaluate
%derivatives with respect to them. Since the radius is inferred from the
%points, the radius is taken as a component too.
%
%EXAMPLE 1:
%In this example, we show that the gradient returned by this function
%agrees with one obtained by numerically differentiating the
%Cart2RSpherGnomonicProj function. The reative error implies more than 11
%digits of agreement.
% latLonRef=deg2rad([20.756113;-156.010933]);
% latLonPt=deg2rad([19.456233;-154.823617]);
% r=osculatingSpher4LatLon(latLonRef);
% xyzCart=ellips2Cart(latLonPt,r,0);
% f=@(xyz)Cart2RSpherGnomonicProj(xyz,latLonRef);
% JacobNumDiff=numDiff(xyzCart,f,3,[],1e-4*norm(xyzCart));
% Jacob=calcSpherGnomonicCartJacob(xyzCart,latLonRef);
% RelErr=max(abs((Jacob(:)-JacobNumDiff(:))./JacobNumDiff(:)))
%
%EXAMPLE 2:
%This is the same as example 1 ,except the reference point is the same as
%the point at which the jacbian should be evaluated. Because one derivative
%is numerically zero, we display the absolute rather than relative error.
%Again, there is good agreement.
% latLonPt=deg2rad([20.756113;-156.010933]);
% r=osculatingSpher4LatLon(latLonPt);
% xyzCart=ellips2Cart(latLonPt,r,0);
% f=@(xyz)Cart2RSpherGnomonicProj(xyz,latLonPt);
% JacobNumDiff=numDiff(xyzCart,f,3,[],1e-4*norm(xyzCart));
% Jacob=calcSpherGnomonicCartJacob(xyzCart);
% AbsErr=max(abs(Jacob(:)-JacobNumDiff(:)))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

x=xyzPt(1,:);
y=xyzPt(2,:);
z=xyzPt(3,:);
r2=sum(xyzPt.*xyzPt,1);
r=sqrt(r2);

drdx=x/r;
drdy=y/r;
drdz=z/r;

if(nargin<2||isempty(latLonRef))
    %If the xyzPts are also taken as the reference.
    
    rxy=sqrt(x^2+y^2);
    
    dxGdx=-y/rxy;
    dxGdy=x/rxy;
    dxGdz=0;
    rrxy=r*rxy;

    dyGdx=-x*z/rrxy;
    dyGdy=-y*z/rrxy;
    dyGdz=rxy/r;

    J=[drdx,  drdy drdz;
       dxGdx,dxGdy,dxGdz;
       dyGdx,dyGdy,dyGdz];
    return
end

%The sine and cosine of the reference latitude.
sinPhi1=sin(latLonRef(1));
cosPhi1=cos(latLonRef(1));
%The sine and cosine of the reference longitude.
sinLambda0=sin(latLonRef(2));
cosLambda0=cos(latLonRef(2));
cos2Lambda0=cos(2*latLonRef(2));



denom=(r*(sinPhi1*z+cosPhi1*x*cosLambda0+cosPhi1*y*sinLambda0)^2);
dxGdx=(-cosPhi1*y*r2+cosPhi1*x^2*y*cos2Lambda0-sinPhi1*z*(2*x^2+y^2+z^2)*sinLambda0+x*cosLambda0*(sinPhi1*y*z+cosPhi1*(-x+y)*(x+y)*sinLambda0))/denom;
dxGdy=(cosPhi1*x*(x^2+y^2+z^2+y^2*cos2Lambda0)-sinPhi1*x*y*z*sinLambda0+cosLambda0*(sinPhi1*z*(x^2+2*y^2+z^2)+cosPhi1*y*(-x+y)*(x+y)*sinLambda0))/denom;
dxGdz=((y*cosLambda0-x*sinLambda0)*(-sinPhi1*(x^2+y^2)+cosPhi1*x*z*cosLambda0+cosPhi1*y*z*sinLambda0))/denom;
dyGdx=(-cosPhi1*sinPhi1*x^3*cosLambda0^2+x*(sinPhi1*z+cosPhi1*y*sinLambda0)*(cosPhi1*z-sinPhi1*y*sinLambda0)-cosLambda0*(cosPhi1^2*z*(y^2+z^2)+sinPhi1^2*z*(2*x^2+y^2+z^2)+2*cosPhi1*sinPhi1*x^2*y*sinLambda0))/denom;
dyGdy=(y*(sinPhi1*z+cosPhi1*x*cosLambda0)*(cosPhi1*z-sinPhi1*x*cosLambda0)-(cosPhi1^2*z*(x^2+z^2)+sinPhi1^2*z*(x^2+2*y^2+z^2)+2*cosPhi1*sinPhi1*x*y^2*cosLambda0)*sinLambda0-cosPhi1*sinPhi1*y^3*sinLambda0^2)/denom;
dyGdz=(cosPhi1*sinPhi1*z^3+(x*cosLambda0+y*sinLambda0)*((cosPhi1^2+sinPhi1^2)*(x^2+y^2)+2*cosPhi1^2*z^2-cosPhi1*sinPhi1*z*(x*cosLambda0+y*sinLambda0)))/denom;

J=[drdx,  drdy drdz;
   dxGdx,dxGdy,dxGdz;
   dyGdx,dyGdy,dyGdz];

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
