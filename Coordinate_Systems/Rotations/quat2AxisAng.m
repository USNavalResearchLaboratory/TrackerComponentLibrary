function [u,theta]=quat2AxisAng(q)
%%QUAT2AXISANG Convert a unit quaternion to the axis and angle of the
%              rotation represented by the quaternion. The handedness of
%              the angle matches the handedness of the quaternion.
%
%INPUTS: q A 4X1 unit quaternion corresponding to the rotation matrix. The
%          quaternion is ordered in terms of hypercomplex numbers as
%          q(1)+i*q(2)+j*q(3)+k*q(4). Since this is a unit-magnitude
%          quaternion, the ordering of the elements can also be expressed
%          as [cos(theta/2);sin(theta/2)u'] where u is a unit vector for
%          the axis of rotation and theta rotation angle about that unit
%          vector, according to whatever handedness is used with the
%          quaternion.
%
%OUTPUTS: u A unit vector representing a rotation axis.
%     theta The angle in radians by which a vector is to be rotated about
%           the axis u.
%
%The relation between an rotation about a given axis by a given angle and a
%quaternion is described in [1].
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. - Dec. 1993.
%    wherein the convention for the quaternions is left-handed.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Choose the sign to make the cosine term in the quaternion positive.
if(q(1)<0)
    q=-q;
end

%Go from the quaternion to the rotation angles.
uMag=norm(q(2:4));
u=q(2:4)/uMag;

%Deal with zero rotations.
if(any(~isfinite(u)))
    u=[1;0;0];
end

theta=2*acos(q(1));
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
