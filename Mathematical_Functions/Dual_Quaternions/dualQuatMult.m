function dq=dualQuatMult(dq2,dq1,handed)
%%DUALQUATMULT Multiply two dual quaternions together in the order
%              dq=dq2*dq1. The multiplication of two dual quaternions is
%              not commutative. Dual quaternions are often used for
%              simultaneously representing orientation and position. A dual
%              quaternion consists of two parts dq=q1+eps*q2, where q1 and
%              q2 are quaternions and eps is a dual number. Dual numbers
%              are such that eps^2=0 and have various rules for
%              multiplication with complex numbers. The quaternions making
%              up a dual quaternion can either be left or right handed.
%
%INPUTS:dq2 A dual quaternion represented as a 4X2 matrix dq2(:,1) is the
%           non-dual quaternion and dq2(:,2) is the dual quaternion (the
%           ones times eps, the dual number). The elements of each
%           quaternion are ordered q(1,1)+i*q(2,1)+j*q(3,1)+k*q(4,1), where
%           i,j and k are the typical hypercomplex numbers.
%       dq1 The second (right side of multiplication) dual quaternion by
%           which dq2 will be multiplied.
%     handed The handedness of the quaternions within the dual quaternions.
%            If omitted or an empty matrix is passed,, it is assumed that
%            the quaternions in the dual quaternion are right-handed (the
%            standard). Possible values are
%            'right' The default if omitted. The quaternion multiplication
%                    is assumed right-handed (standard).
%            'left'  The quaternion multiplication is assumed left-handed.
%                    
%OUTPUT: dq The dual quaternion representing the product of dq2 and dq1.
%           dq=dq2*dq1.
%
%Dual quaternions and common operations including multiplication are
%discussed in [1].
%
%REFERENCES:
%[1] B. Kenwright, "A beginners guide to dual-quaternions: What they are,
%    how they work, and how to use them for 3D character hierarchies," in
%    Proceedings of the 20th International Conference on Computer Graphics,
%    Visualization and Computer Vision, Prague, Czech Republic, 24-27 Jun.
%    2012.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(handed))
    handed='right';
end

a1=dq2(:,1);
a2=dq2(:,2);
b1=dq1(:,1);
b2=dq1(:,2);

dq=[quatMult(a1,b1,handed),quatMult(a1,b2,handed)+quatMult(a2,b1,handed)];
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
