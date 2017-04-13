function J=spherAngGradient(xG,systemType,lRx,M)
%%SPHERANGGRADIENT Determine the gradient of a 3D spherical azimuth and
%          elevation measurement with respect to 3D position. Higher order
%          gradient terms are not provided and are zero. Relativity and
%          atmospheric effects are not taken into account.
%
%INPUTS: x The 3X1 target position vector in the global coordinate system
%          with [x;y;z] components.
% systemType An optional parameter specifying the axes from which the
%          angles are measured in radians. Possible vaues are
%          0 (The default if omitted) Azimuth is measured counterclockwise
%            from the x-axis in the x-y plane. Elevation is measured up
%            from the x-y plane (towards the z-axis). This is consistent
%            with common spherical coordinate systems for specifying
%            longitude (azimuth) and geocentric latitude (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%      lRx The 3X1  position vector of the receiver. If omitted, the
%          receiver is placed at the origin.
%        M A 3X3 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. If
%          omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: J A 2X3 Jacobian matrix where the rows are [azimuth;elevation] in
%           that order and the columns take the derivative of the row
%           components with respect to [x,y,z] in that order.
%
%The derivatives can be computed in a straightforward manner from
%the basic relation between spherical and Cartesian coordinates, which is
%given in Chapter 14.4.4.1 of [1], among other sources.
%
%Note that singularities exist at the poles; that is when the elevation is
%+/-(pi/2).
%
%REFERENCES:
%[1] R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%    to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%    Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(M))
   M=eye(3,3);
end

if(nargin<3||isempty(lRx))
   lRx=zeros(3,1); 
end

if(nargin<2||isempty(systemType))
   systemType=0; 
end

%Transform from global coordinates to local coordinates.
xLocal=M*(xG(1:3)-lRx(1:3));

x=xLocal(1);
y=xLocal(2);
z=xLocal(3);

r=norm(xLocal);

J=zeros(2,3);
switch(systemType)
    case 0
        %Derivatives with respect to x.
        J(1,1)=-y/(x^2+y^2);
        J(2,1)=-x*z/(r^2*sqrt(x^2+y^2));

        %Derivatives with respect to y.
        J(1,2)=x/(x^2+y^2);
        J(2,2)=-y*z/(r^2*sqrt(x^2+y^2));

        %Derivatives with respect to z.
        J(1,3)=0;
        J(2,3)=sqrt(x^2+y^2)/r^2;
    case 1
        %Derivatives with respect to x.
        J(1,1)=z/(z^2+x^2);
        J(2,1)=-x*y/(r^2*sqrt(z^2+x^2));
        
        %Derivatives with respect to y.
        J(1,2)=0;
        J(2,2)=sqrt(z^2+x^2)/r^2;
        
        %Derivatives with respect to z.
        J(1,3)=-x/(z^2+x^2);
        J(2,3)=-z*y/(r^2*sqrt(z^2+x^2));
    otherwise
        error('Invalid system type specified.')
end

%Rotate from local back to global coordinates.
J=J*M;

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
