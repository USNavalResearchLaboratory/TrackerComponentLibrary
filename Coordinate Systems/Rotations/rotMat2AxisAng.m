function [u,theta,q]=rotMat2AxisAng(R,handed)
%%ROTMAT2AXISANG Any series of rotations in 3D can be expressed in terms of
%                a single rotation about a certain axis by a certain angle.
%                This function takes a rotation matrix and provides the
%                axis and angle of the rotation as well as a rotation
%                quaternion, if desired. The quaterion and rotation angle
%                can be expressed as either a right-handed or a left-handed
%                rotation. The default is right-handed.
%
%INPUTS:  R       A 3X3 orthonormal real rotation matrix.
%
%OUTPUTS: u,theta  The 3X1 unit vector u and rotation angle theta in
%                  radians such that to rotate a vector v, one could call
%                  rotateVector(v,u,theta,handed) instead of multiplying
%                  M*v. theta is the angle by which one would have to
%                  perform a rotation about the axis u (with the handedness
%                  of the rotation given by handed) to rotate a vector v
%                  the same as if one evaluated M*v.
%                  Also, note that R=axisAng2RotMat(u,theta,handed).
%         q        A quaternion corresponding to the axis-angle pair where
%                  q=[cos(theta/2);sin(theta/2)*u]; The first element is
%                  what many textbooks refer to as q0 or q4. The quaternion
%                  is ordered in terms of hypercomplex numbers as
%                  q(1)+i*q(2)+j*q(3)+k*q(4). The handedness of the
%                  quaternion matches handed (and the handedness of theta).
%           handed  The handedness of the rotation angle. If omitted, it
%                   is assumed that the rotation is right-handed
%                   (the standard). Possible values are
%            'right' The default if omitted. The rotation is right-handed.
%            'left'  The rotation is left-handed. The rotation angle is
%                    counterclockwise when one is looking in the direction
%                    that the rotation axis points.
%
%The function calls rotMat2Quat to get the quaternion for the rotation and
%then calls quat2AxisAng to get the axis and angle.
%
%March 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    handed='right';
end

q=rotMat2Quat(R,handed);
[u,theta]=quat2AxisAng(q);
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
