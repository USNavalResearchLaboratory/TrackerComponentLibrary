function J=calcSpherJacob(point,systemType)
%%CALCSPHERJACOB  Compute the Jacobian matrix for a point in spherical
%                 [range;azimuth;elevation] coordinates. This function
%                 differs from the function calcSpherRRJacob in that only a
%                 1-way monostatic range in considered, range rate is
%                 not available, and the point here is given in spherical
%                 coordinates, whereas the function calcSpherRRJacob takes
%                 the point in global Cartesian coordinates.
%
%INPUTS: point   A point in the format [range;azimuth;elevation], where the
%                two angles are given in radians.
%     systemType An optional parameter specifying the axes from which
%                the angles are measured. Possible vaues are
%                   0 (The default if omitted) Azimuth is measured
%                     counterclockwise from the x-axis in the x-y plane.
%                     Elevation is measured up from the x-y plane (towards
%                     the z-axis). This is consistent with common spherical
%                     coordinate systems for specifying longitude (azimuth)
%                     and geocentric latitude (elevation).
%                   1 Azimuth is measured counterclockwise from the z-axis
%                     in the z-x plane. Elevation is measured up from the
%                     z-x plane (towards the y-axis). This is consistent
%                     with some spherical coordinate systems that use the z
%                     axis as the boresight direction of the radar.
%
%OUTPUTS: J     The 3X3 Jacobian matrix. Each row is a components of range,
%               azimuth and elevation (in that order by row) with
%               derivatives taken with respect to [x,y,z] by column.
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
%December 2013 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    systemType=0;
end

CartPoint=spher2Cart(point,systemType);
r=point(1);
J=zeros(3,3);

x=CartPoint(1);
y=CartPoint(2);
z=CartPoint(3);

switch(systemType)
    case 0
        %Derivatives with respect to x.
        J(1,1)=x/r;
        J(2,1)=-y/(x^2+y^2);
        J(3,1)=-x*z/(r^2*sqrt(x^2+y^2));

        %Derivatives with respect to y.
        J(1,2)=y/r;
        J(2,2)=x/(x^2+y^2);
        J(3,2)=-y*z/(r^2*sqrt(x^2+y^2));

        %Derivatives with respect to z.
        J(1,3)=z/r;
        J(2,3)=0;
        J(3,3)=sqrt(x^2+y^2)/r^2;
    case 1
        %Derivatives with respect to x.
        J(1,1)=x/r;
        J(2,1)=z/(z^2+x^2);
        J(3,1)=-x*y/(r^2*sqrt(z^2+x^2));
        
        %Derivatives with respect to y.
        J(1,2)=y/r;
        J(2,2)=0;
        J(3,2)=sqrt(z^2+x^2)/r^2;
        
        %Derivatives with respect to z.
        J(1,3)=z/r;
        J(2,3)=-x/(z^2+x^2);
        J(3,3)=-z*y/(r^2*sqrt(z^2+x^2));
    otherwise
        error('Invalid system type specified.')
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
