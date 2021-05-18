function [u,c]=getEllipsHarmAxes(pointHarmon,E)
%%GETELLIPSHARMAXES Find the Cartesian unit vectors corresponding to the
%                   normalized  gradients of the components of a point in
%                   ellipsoidal harmonic coordinates. Also, get the
%                   magnitudes of the gradients.
%
%INPUTS: pointHarmon A point given in terms of ellipsoidal harmonic reduced
%                    latitude and longitude in radians and a semi-major
%                    axis in meters.
%                  E The linear eccentricity defining the ellipsoidal
%                    harmonic coordinate system. If this parameter is
%                    omitted, then the linear eccentricity of the WGS84
%                    reference ellipsoid is used.
%
%OUTPUTS: u u(:,1), u(:,2) and u(:,3) are the unit vectors in the
%           respective directions of the respective gradients of reduced
%           latitude, longitude and the semi-major axis.
%         c c(1), c(2) and c(3) are the respective magnitudes of the
%           derivative of the Cartesian position with respect to reduced
%           latitude, longitude and the semi-major axis.
%
%The ellipsoidal harmonic coordinate system is described in Chapter 1.15
%of [1] and is an orthogonal coordinate system.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(E))
    a=Constants.WGS84SemiMajorAxis;
    f=Constants.WGS84Flattening;
    b=a*(1 - f);
    E=sqrt(a^2-b^2);
end

    beta=pointHarmon(1);
    lambda=pointHarmon(2);
    u=pointHarmon(3);
    
dBeta=[-sqrt(E^2+u^2)*sin(beta)*cos(lambda);
       -sqrt(E^2+u^2)*sin(beta)*sin(lambda);
        u*cos(beta)];
c1=norm(dBeta);
dBeta=dBeta/c1;
      
du=[u*cos(beta)*cos(lambda)/sqrt(E^2+u^2);
    u*cos(beta)*sin(lambda)/sqrt(E^2+u^2);
    sin(beta)];
c3=norm(du);
du=du/c3;

dLambda=[-sqrt(E^2+u^2)*cos(beta)*sin(lambda);
          sqrt(E^2+u^2)*cos(beta)*cos(lambda);
          0];
c2=norm(dLambda);

%If the point is too close to the poles, then it is possible that c2 is
%nearly equal to zero. However, a normalized dLambda can just be found by
%orthogonality: it is orthogonal to dBeta and du.
dLambda=cross(dBeta,du);

u=[dBeta,dLambda,du];
c=[c1;c2;c3];
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
