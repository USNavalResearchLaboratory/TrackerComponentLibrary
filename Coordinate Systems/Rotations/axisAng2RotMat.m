function M=axisAng2RotMat(u,theta,handed)
%%AXISANG2ROTMAT Get a rotation matrix to rotate a 3D vector an angle of
%                theta counterlockwise (right-handed) or clockwise
%                (left-handed) about an axis given by the unit vector u.
%
%INPUTS: u A unit vector representing the axis about which a vector should
%          be rotated using the rotation matrix.
%    theta The angle in radians by which a vector is to be rotated is to be
%          rotated about the axis u when multiplied by the rotation matrix.
%   handed The handedness of the rotation angle. If omitted, it is assumed
%          that the rotation is right-handed (the standard). Possible
%          values are
%          'right' The default if omitted. The rotation is right-handed.
%           'left' The rotation is left-handed. The rotation angle is
%                  counterclockwise when one is looking in the same
%                  direction that the rotation axis points.
%
%OUTPUTS: M The rotation matrix such that M*v rotates the vector v by the
%           desired rotation.
%
%To rotate a 3D vector x an angle of theta about the axis u, simply
%evaluate M*x.
%
%The rotation matrix is obtained by obtaining the unit quaternion for the
%rotation and turning it into the corresponding rotation matrix as
%described in  [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(handed))
    handed='right';
end

q=axisAng2Quat(u,theta);
M=quat2RotMat(q,handed);

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
