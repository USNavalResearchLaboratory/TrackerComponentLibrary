function Q=quatDivL(quat2,quat1,handed)
%%QUATDIVL Left-divide the quaternion quat1 by the quaternion quat2.
%          Quaternion algebra is not commutative, so a difference exists
%          between left and right division. Quaternions are an extension
%          of complex numbers. The quaternions can be either left or
%          right-handed.
%
%INPUTS: quat2 A 4XN set of N quaternions, where the first element in each
%             column is the scalar part of the quaternion (sometimes called
%             q0 or q4) and the next three elements are the (hypercomplex)
%             vector part. That is, the hypercomplex quaternion given by
%             quat2(:,1) can be written in hypercomplex, non-vector form
%             as quat2(1,1)+i*quat2(2,1)+j*quat2(3,1)+k*quat2(4,1), where
%             i, j, and k are all roots of -1.
%       quat1 A 4XN set of N quaternions by which each quaternion in quat1
%             is left-divided by the corresponding quaternion in quat2.
%     handed The handedness of the quaternions. If omitted or an empty
%            matrix is passed,, it is assumed that the quaternions are
%            right-handed (the standard). Possible values are
%            'right' The default if omitted. The quaternion multiplication
%                    is assumed right-handed (standard).
%            'left'  The quaternion multiplication is assumed left-handed.
%                    This is used in someplaces, including the reference
%                    from Shuster, below.
%
%OUTPUTS: Q The 4XN matrix of quaternions divided such that
%           Q(:,i)=(1/quat2(:,i))*quat1(:,i).
%
%Properties of quaternions including multiplication and division are
%described in  [1]. Additional aspects of quaternions, including their use
%in a left-handed system, are in [2].
%
%A quaternion of the form q(1)+i*q(2)+j*q(3)+k*q(4) that obeys right-handed
%multiplication and division rules supports the following rules for
%multiplication of i, j, and k, where an item in a row is multiplied by an
%item in the column to get the result:
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
%[1] Weisstein, Eric W. "Quaternion." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Quaternion.html
%[2] M. D. Shuster, "A survey of attitude representations," The Journal of
%    the Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct.-Dec.
%    1993.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(handed))
    handed='right';
end

Q=quatMult(quatInv(quat2),quat1,handed);
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
