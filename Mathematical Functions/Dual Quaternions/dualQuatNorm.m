function normVal=dualQuatNorm(dq)
%%DUALQUATNORM Compute the norm of a dual quaternion. This is defined as
%           dq*dq^*, where dq^* is the conjugate of the dual quaternion.
%           The result is a dual quaternion with no hypercomplex components
%           (but it might have both a real and a dual component). Dual
%           quaternions are often used for simultaneously representing
%           orientation and position. A dual quaternion consists of two
%           parts dq=q1+eps*q2, where q1 and q2 are quaternions and eps is
%           a dual number. Dual numbers are such that eps^2=0 and have
%           various rules for multiplication with complex numbers. Note
%           that the norm does not depend on the handedness of the
%           quaternion algebra, even though the dualQuatMult and quatMult
%           functions do.
%
%INPUTS: dq A dual quaternion represented as a 4X2 matrix dq(:,1) is the
%           non-dual quaternion and dq(:,2) is the dual quaternion (the
%           ones times eps, the dual number). The elements of each
%           quaternion are ordered q(1,1)+i*q(2,1)+j*q(3,1)+k*q(4,1), where
%           i,j and k are the typical hypercomplex numbers.
%
%OUTPUTS: normVal The value of the norm of the dual quaternion. This is a
%                 dual quaternion with no complex parts, but it might have
%                 a dual part (normVal(1,:) can be nonzero). A unit dual
%                 quaternion has only a real part and no dual part
%                 (normVal(1,1)=1,normVal(2,1)=0) nor complex parts.
%
%Dual quaternions and common operations including computation of the norm
%of the dual quaternion are discussed in [1]. The norm of a standard
%(non-dual) quaternion is generally just taken to be the same as the norm
%of a typical 4X1 vector.
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

normVal=dualQuatMult(dq,dualQuatConj(dq));

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
