function vRot=rotateVectorQuat(v,q,handed)
%%ROTATEVECTORQUAT Rotate three-dimensional vectors v by the rotation
%                  specified by the unit quaternion q. The quaternion can
%                  be specified in either a left or a right-handed system.
%                  A rotation quaternion can be obtained, for example,
%                  using the axisAng2Quat function.
%
%INPUTS: v The 3XN set of N vectors that are to be rotated.
%        q A 4X1 unit quaternion corresponding to a rotation. The ordering
%          of the elements corresponds to the hypercomplex decomposition
%          as q(1)+i*q(2)+j*q(3)+k*q(4), where i, j, and k are all roots of
%          -1.
%   handed The handedness of the quaternion. If omitted, it is assumed
%          that the quaternion is right-handed (the standard). Possible
%          values are
%          'right' The default if omitted. The quaternion multiplication is
%                  assumed right-handed (standard).
%          'left'  The quaternion multiplication is assumed left-handed.
%                  This is used in someplaces, including the reference
%                  from Shuster, below.
%
%OUTPUTS: vRot The 3XN set of vectors v rotated according to the rotation
%              implied by the unit quaternion with the given handedness.
%
%The use of unit-magnitude quaternions for rotation is discussed in [1],
%where left-handed quaternions are considered.
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
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(handed))
   handed='right'; 
end

N=size(v,2);
vRot=zeros(3,N);%Allocate space

%The conjugate quaternion.
qC=quatConj(q);

for curVec=1:N
    v4=quatMult(quatMult(q,[0;v(:,curVec)],handed),qC,handed);
    vRot(:,curVec)=v4(2:4);
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
