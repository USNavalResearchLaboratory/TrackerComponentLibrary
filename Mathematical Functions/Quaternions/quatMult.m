function Q=quatMult(quat2,quat1,handed)
%%QUATMULT  Multiply two quaternions together. Quaternions are an extension
%           of complex numbers. The type of product performed here is
%           deemed the Grassman product. Unit-magnitude quaternions are
%           often used to represent rotations and the multiplication
%           quat2*quat1 is the orientation obtained by rotating by quat1
%           and then rotating by quat2. Quaternion multiplication is not
%           commutative. The quaternions can be either left or
%           right-handed.
%
%INPUTS:quat2 A 4XN set of N quaternions, where the first element in each
%             column is the scalar part of the quaternion (sometimes called
%             q0 or q4) and the next three elements are the (hypercomplex)
%             vector part. That is, the hypercomplex quaternion given by
%             quat2(:,1) can be written in hypercomplex, non-vector form
%             as quat2(1,1)+i*quat2(2,1)+j*quat2(3,1)+k*quat2(4,1), where
%             i, j, and k are all roots of -1. When considering unit
%             quaternions, this means that it has the form
%             [cos(theta/2);sin(theta/2)u'], where theta is an angle and u
%             is a 3X1 unit vector.
%       quat1 A 4XN set of N quaternions that are to be left-multiplied by
%             the corresponding quaternions in quat2.
%     handed The handedness of the quaternions. If omitted or an empty
%            matrix is passed,, it is assumed that the quaternions are
%            right-handed (the standard). Possible values are
%            'right' The default if omitted. The quaternion multiplication
%                    is assumed right-handed (standard).
%            'left'  The quaternion multiplication is assumed left-handed.
%                    This is used in someplaces, including the reference
%                    from Shuster, below.
%
%OUTPUTS: Q   The 4XN matrix of quaternions products such that
%             Q(:,i)=quat2(:,i)*quat1(:,1), where the multiplication
%             operation is quaternion multiplication.
%
%Properties of quaternions including multiplication are described in [1].
%Details of rotation using unit quaternions are given in [2].
%
%The difference in multiplication between right-handed and left-handed
%quaternions lies in the sign used for the cross product.
%
%A quaternion of the form q(1)+i*q(2)+j*q(3)+k*q(4) that obeys right-handed
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
%[1] Weisstein, Eric W. "Quaternion." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Quaternion.html
%[2] M. D. Shuster, "A survey of attitude representations," The Journal of
%    the Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct.-Dec.
%    1993.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of quaternions.
N=size(quat2,2);

if(nargin<3||isempty(handed))
    handed='right';
end

Q=zeros(4,N);
%The handedness switches the sign of the cross product
switch(handed)
    case 'right'
        for curQuat=1:N
            s1=quat1(1,curQuat);
            v1=quat1(2:4,curQuat);
            s2=quat2(1,curQuat);
            v2=quat2(2:4,curQuat);

            Q(:,curQuat)=[s1*s2-dot(v1,v2);
                          s1*v2+s2*v1+cross(v2,v1)];
        end
    case 'left'
        %Flip the algebra to left-handed for the formula below.
         for curQuat=1:N
            s1=quat1(1,curQuat);
            v1=quat1(2:4,curQuat);
            s2=quat2(1,curQuat);
            v2=quat2(2:4,curQuat);

            Q(:,curQuat)=[s1*s2-dot(v1,v2);
                          s1*v2+s2*v1-cross(v2,v1)];
        end
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
