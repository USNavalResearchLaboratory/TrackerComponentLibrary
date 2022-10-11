function q=rotMat2Quat(R,handed)
%%ROTMAT2QUAT Get a unit quaternion corresponding to a particular rotation
%             matrix. The quaternion can be chosen to support standard
%             right-handed quaternion multiplication rules, or non-standard
%             left-handed rules that some authors choose to use.
%
%INPUTS:  M A 3X3 orthonormal real rotation matrix.
%    handed The handedness of the quaternion. If omitted, it is assumed
%           that the quaternion is right-handed (the standard). Possible
%           values are:
%           'right' The default if omitted. The quaternion multiplication
%                   is assumed right-handed (standard).
%           'left'  The quaternion multiplication is assumed left-handed.
%                   This is used in someplaces, including the reference
%                   from Shuster, below.
%
%OUTPUTS: q A 4X1 unit quaternion corresponding to the rotation matrix. The
%           quaternion is ordered [cos(theta/2);sin(theta/2)u'] where u is
%           a unit vector for the axis of rotation and theta is the
%           counterclockwise (right-handed) or clockwise (left-handed)
%           rotation angle about that unit vector.
%
%The formulae for converting a rotation matrix are from [1], where a minor
%change has been performed to support both right and left-handed
%quaternions. Both q and -q represent the same rotations. Here, the
%solution is chosen so that when considering a left-handed rotation, the
%sign of the largest magnitude element is positive. When considering a
%right handed rotation, the sign of the largest element may or may not be
%positive.
%
%A quaternion of form q(1)+i*q(2)+j*q(3)+k*q(4) that obeys right-handed
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
%[1] W. F. Phillips, C. E. Hailey, and G. A. Gebert, "Review of attitude
%    representations used for aircraft kinematics," Journal of Aircraft,
%    vol. 38, no. 4, pp. 718-737, Jul. - Aug. 2001.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(handed))
    handed='right';
end

R11=R(1,1);
R22=R(2,2);
R33=R(3,3);
R12=R(1,2);
R21=R(2,1);
R13=R(1,3);
R31=R(3,1);
R23=R(2,3);
R32=R(3,2);

%The square of the quaternion...
q2=(1/4)*[1+R11+R22+R33;
          1+R11-R22-R33;
          1-R11+R22-R33;
          1-R11-R22+R33];

q2Max=max(q2);

if(q2Max==q2(1))
    q0=0.5*sqrt(1+R11+R22+R33);
    q1=(1/(4*q0))*(R23-R32);
    q2=(1/(4*q0))*(R31-R13);
    q3=(1/(4*q0))*(R12-R21);
elseif(q2Max==q2(2))
    q1=0.5*sqrt(1+R11-R22-R33);
    q0=(1/(4*q1))*(R23-R32);
    q2=(1/(4*q1))*(R12+R21);
    q3=(1/(4*q1))*(R31+R13);
elseif(q2Max==q2(3))
    q2=0.5*sqrt(1-R11+R22-R33);
    q0=(1/(4*q2))*(R31-R13);
    q1=(1/(4*q2))*(R12+R21);
    q3=(1/(4*q2))*(R23+R32);
else
    q3=0.5*sqrt(1-R11-R22+R33);
    q0=(1/(4*q3))*(R12-R21);
    q1=(1/(4*q3))*(R31+R13);
    q2=(1/(4*q3))*(R23+R32);
end

%The above implementation provides a left-handed quaternion, so the sign of
%the vector part needs to be flipped if the standard right-handed system is
%used.
switch(handed)
    case 'right'
        q=[q0;-q1;-q2;-q3];
    case 'left'
        q=[q0;q1;q2;q3];
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
