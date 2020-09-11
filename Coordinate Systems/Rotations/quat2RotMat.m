function R=quat2RotMat(q,handed)
%%QUAT2ROTMAT Turn a unit quaternion into an equivalent rotation matrix.
%             The multiplication rules for the quaternion algebra can be
%             chosen to support standard right-handed quaternion rotation,
%             or non-standard left-handed rotations that some authors use.
%
%INPUTS: q A 4X1 unit quaternion corresponding to the rotation matrix. The
%          quaternion is ordered [cos(theta/2);sin(theta/2)u'] where u is a
%          unit vector for the axis of rotation and theta is the rotation
%          angle about that unit vector according to the specified
%          handedness. The ordering of the elements corresponds to the
%          hypercomplex decomposition q(1)+i*q(2)+j*q(3)+k*q(4), where i,
%          j, and k are all roots of -1.
%   handed The handedness of the quaternion. If omitted, it is assumed that
%          the quaternion is right-handed (the standard). Possible values
%          are:
%          'right' The default if omitted. The quaternion multiplication is
%                  assumed right-handed (standard).
%          'left'  The quaternion multiplication is assumed left-handed.
%                  This is used in someplaces, including the reference from
%                  Shuster, below.
%
%OUTPUTS: R A 3X3 orthonormal rotation matrix.
%
%If q does not have unit magnitude, then R will not be a rotation matrix.
%
%The formula for turning a quaternion into a left-handed rotation matrix is
%given in [1]. However, right-handed quaternion algebra, as used in [2] is
%far more common.
%
%A quaternion form q(1)+i*q(2)+j*q(3)+k*q(4) that obeys right-handed
%multiplication rules supports the following rules for multiplication of i,
%j, and k, where an item in a row is multiplied by an item in the column to
%get the result:
%  i,  j, k
%i -1, k,-j
%j -k,-1, i
%k  j,-i,-1
%On the other hand, left-handed multiplication rules flip the signs of the
%off-diagonal terms:
%  i,  j, k
%i -1,-k, j
%j  k,-1,-i
%k -j, i,-1
%
%REFERENCES:
%[1] M. D. Shuster, "A survey of attitude representations," The Journal of
%    Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct. -Dec. 1993.
%[2] Weisstein, Eric W. "Quaternion." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Quaternion.html
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(handed))
    handed='right';
end

switch(handed)
    case 'left'
        %The rotation matrix for a left-handed quaternion algebra, as in
        %[1].
        R=[q(1)^2+q(2)^2-q(3)^2-q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)),    2*(q(2)*q(4)-q(1)*q(3));
           2*(q(2)*q(3)-q(1)*q(4)),     q(1)^2-q(2)^2+q(3)^2-q(4)^2,2*(q(3)*q(4)+q(1)*q(2));
           2*(q(2)*q(4)+q(1)*q(3)),     2*(q(3)*q(4)-q(1)*q(2)),    q(1)^2-q(2)^2-q(3)^2+q(4)^2];
    case 'right'
        %The rotation matrix for a right-handed quaternion algebra is the
        %transpose of the matrix from [1].
        R=[q(1)^2+q(2)^2-q(3)^2-q(4)^2, 2*(q(2)*q(3)-q(1)*q(4)),    2*(q(2)*q(4)+q(1)*q(3));
           2*(q(2)*q(3)+q(1)*q(4)),     q(1)^2-q(2)^2+q(3)^2-q(4)^2,2*(q(3)*q(4)-q(1)*q(2));
           2*(q(2)*q(4)-q(1)*q(3)),     2*(q(3)*q(4)+q(1)*q(2)),    q(1)^2-q(2)^2-q(3)^2+q(4)^2];
    otherwise
        error('Invalid handedness provided.')
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
