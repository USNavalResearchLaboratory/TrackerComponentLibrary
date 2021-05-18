function q=axisAng2Quat(u,theta)
%%AXISANG2QUAT Get a unit quaternion representative of a rotation of a 3D
%              vector an angle of theta counterlockwise (right-handed) or
%              clockwise (left-handed) about an axis given by u. The
%              handedness of theta matches the handedness of the
%              quaternion produced.
%
%INPUTS: u A unit vector representing the axis about which a vector should
%          be rotated using the rotation matrix.
%    theta The angle in radians by which a vector is to be rotated about
%          the axis u when multiplied by the rotation matrix. The
%          handedness of the rotation angle matches the handedness of the
%          quaternion produced.
%
%OUTPUTS: q A 4X1 unit quaternion corresponding to the supplied axis and
%           angle. The handedness of the quaternion (which is important to
%           know when using the quaternion with other functions, such as
%           quatMult) is the same as that of the angle. The quaternion is
%           ordered in terms of hypercomplex numbers as
%           q(1)+i*q(2)+j*q(3)+k*q(4).
%
%Quaternions and rotations are discussed in [1].
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. -Dec. 1993.
%    wherein the convention for the quaternions is left-handed.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

q=[cos(theta/2);sin(theta/2)*u];

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
