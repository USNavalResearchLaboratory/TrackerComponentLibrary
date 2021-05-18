function [u,c]=getENUAxes(plhPoint,justVertical,a,f)
%%GETENUAXES Find the East-North-Up unit vectors at the given point as well
%            as the magnitudes of the derivatives of a position vector with
%            respect to latitude, longitude and height.
%
%INPUTS: plhPoint The point at which the axes are to be found given in
%                 terms of [latitude;longitude] (with an assumed
%                 ellipsoidal height of 0) or [latitude;longitude;height]
%                 with the geodetic latitude and longitude in radians and
%                 the height in meters. The latitude should be between
%                 -pi/2 and pi/2. The height does not change the unit
%                 direction vectors u but it does change c.
%    justVertical An optional parameter. If this is given and is true, then
%                 u and c only for the Up direction will be returned. The
%                 default is false. For the up direction, c is always 1.
%               a The semi-major axis of the reference ellipsoid. This is
%                 only used if the c output is requested and
%                 justVertical=false. If this argument is omitted or an
%                 empty matrix is passed, the value in
%                 Constants.WGS84SemiMajorAxis is used.
%               f The flattening factor of the reference ellipsoid. This is
%                 only used if the c output is requested and
%                 justVertical=false. If this argument is omitted, the
%                 value in Constants.WGS84Flattening is used.
%
%OUTPUTS: u u(:,1), u(:,2), and u(:,3) are respectively the East, North and
%           Up unit vectors.
%         c c(1), c(2), and c(3) are the respective magnitudes of the
%           derivative of the Cartesian position with respect to latitude,
%           longitude, and height.
%
%The local East-North-Up coordinate system is an orthonormal coordinate
%system given by the normalized derivatives of a Cartesian position vector
%taken with respect to longitude, geodetic, latitude and height. It is
%described in [1], where expressions for the ENU vectors are given.
%However, when just considering the vector directions and not the
%magnitudes, the dependence on a and f goes away, as given in Appendix C of
%[2].
%
%The longitude coordinate determines the orientation of the axes at the
%poles.
%
%If only the up vector is desired, include the second parameter set to
%true.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating Aerial Targets in 3D Accouting for the
%    Earth's Curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, pp. 31-57, Jun 2015.
%[2] M. S. Grewal, R. Weill, Lawrence, and A. P. Andrews, Global
%    Positioning Systems, Inertial Navigation, and Integration, 2nd ed.
%    Hoboken: Wiley-Interscience, 2007.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The latitude
phi=plhPoint(1);
%The longitude
lambda=plhPoint(2);

if(nargin<2||isempty(justVertical))
    justVertical=false;
end

sinP=sin(phi);
cosP=cos(phi);
sinL=sin(lambda);
cosL=cos(lambda);

%u3 is dr/dh (Up)
u3=[cosP*cosL;
    cosP*sinL;
    sinP];
c3=1;

if(justVertical==false)
    %u1 is dr/dlambda, normalized (East).
    u1=[-sinL;
         cosL;
            0];
    %u2 is dr/dphi, normalized (North)
    u2=[-cosL*sinP;
        -sinL*sinP;
              cosP];
    u=[u1,u2,u3];
    if(nargout>1)
        %If the vector magnitudes are desired.
        if(nargin<4||isempty(f))
            f=Constants.WGS84Flattening;
        end

        if(nargin<3||isempty(a))
            a=Constants.WGS84SemiMajorAxis;
        end
        
        if(length(plhPoint)==2)
            h=0;
        else
            %The height
            h=plhPoint(3);
        end
        
        %The square of the first numerical eccentricity
        e2=2*f-f^2;
        %The normal radius of curvature.
        Ne=a/sqrt(1-e2*sinP^2);
        %The derivative of the normal radius of curvature with respect to phi.
        dNedPhi=a*e2*cosP*sinP/(1-e2*sinP^2)^(3/2);

        %u1 is dr/dlambda, normalized (East).
        c1=norm([-(Ne+h)*cosP*sinL;
                  (Ne+h)*cosP*cosL;
                                0]);
    
        %u2 is dr/dphi, normalized (North)
        c2=norm([(cosP*dNedPhi-(Ne+h)*sinP)*cosL;
                 (cosP*dNedPhi-(Ne+h)*sinP)*sinL;
                 (Ne*(1-e2)+h)*cosP+(1-e2)*dNedPhi*sinP]);

        c=[c1;c2;c3];
    end
else
    u=u3;
    c=c3;
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
