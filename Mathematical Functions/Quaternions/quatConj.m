function qC=quatConj(q)
%%QUATCONJ  Find the conjugate of a quaternion. Quaternions are an
%           extension of complex numbers. The conjugate function does not
%           depend on the handedness of the quaternions used.
%
%INPUTS: q  A 4XN set of N quaternions, each whose conjugate is desired,
%           where the first element in each column is the scalar part of
%           the quaternion (sometimes called q0 or q4) and the next three
%           elements are the (hypercomplex) vector part. That is, the
%           hypercomplex quaternion given by q(:,1) can be written in
%           hypercomplex, non-vector form as
%           q(1,1)+i*q(2,1)+j*q(3,1)+k*q(4,1), where i, j, and k are all
%           roots of -1.
%
%OUTPUTS: qC The 4XN set of complex conjugates of the N quaternions
%
%Properties of quaternions including conjugation are described in [1]. When
%a quaternion has unit magnitude, its conjugate it also its inverse.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Quaternion." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/Quaternion.html
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The conjugate quaternion just flips the sign of the hypercomplex part.
qC=q;
qC(2:4,:)=-qC(2:4,:);
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
