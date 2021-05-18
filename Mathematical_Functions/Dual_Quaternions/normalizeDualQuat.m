function dq1=normalizeDualQuat(dq)
%%NORMALIZEDUALQUAT Given a dual quaternion, normalize it. A normalized
%               dual quaternions is one such that dualQuatNorm produces a
%               dual quaternion result with a real component of 1 and all
%               other components 0. Dual quaternions are often used for
%               simultaneously representing orientation and position. A
%               dual quaternion consists of two parts dq=q1+eps*q2, where
%               q1 and q2 are quaternions and eps is a dual number. Dual
%               numbers are such that eps^2=0 and have various rules for
%               multiplication with complex numbers. The quaternions making
%               up a dual quaternion can either be left or right handed.
%
%INPUTS: dq A dual quaternion represented as a 4X2 matrix dq(:,1) is the
%           non-dual quaternion and dq(:,2) is the dual quaternion (the
%           ones times eps, the dual number). The elements of each
%           quaternion are ordered q(1,1)+i*q(2,1)+j*q(3,1)+k*q(4,1), where
%           i,j and k are the typical hypercomplex numbers.
%
%OUTPUTS: dq1 A normalized version of dq. This is a 4X2 matrix representing
%             a dual quaternion such that x=normalizeDualQuat(dq1) equals a
%             dual quaternion where x(1,1) is 1 and all other elements are
%             zero.
%
%Dual quaternions are discussed in [1]. The normalization routine in [1]
%does not get rid of the dual part of the norm. This function has a second
%step that elmiminates the dual part as well.
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

%This normalizes the nondual part of the dual quaternion, but the dual part
%is still there.
dq1=dq/norm(dq(:,1));
%This gets rid of the dual part.
dq1(:,2)=dq1(:,2)-dq1(:,1)*dot(dq1(:,1),dq1(:,2));

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
