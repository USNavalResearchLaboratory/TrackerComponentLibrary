function [kVec,cNew]=quadEq2CenteredForm(A,bVec,c)
%%QUADEQ2CENTEREDFORM Given a real quadratic equation of the form
%        z=xVec'*A*xVec+bVec'*xVec+c, where xVec and bVec are vectors and
%        A is a symmetric matrix, this function returns kVec and cNew such
%        that z=(xVec-kVec)'*A*(xVec-kVec)-cNew. A cannot be singular. This
%        function can be useful for transforming various expressions for
%        ellipsoids into a centered quadratic form. However, the function
%        is not limited to ellipsoids.
%
%INPUTS: A A real symmetric NXN matrix.
%     bVec A real NX1 vector.
%        c A real scalar value.
%
%OUTPUTS: cNew, kNew These are values such that
%                (xVec-kVec)'*A*(xVec-kVec)-cNew==xVec'*A*xVec+bVec'*xVec+c
%
%Shortening xVec, bVex and kVec to x, b, and k, we say
%kVec=-(1/2)*(A\bVec) and cNew=bVec.'*inv(A)*bVec/4-c. Substituting into
%the expression z=(xVec-kVec)'*A*(xVec-kVec)-cNew, one ends up with the
%expression z=xVec'*A*xVec+bVec'*xVec+c after simplification.
%
%EXAMPLE:
%Here, we convert one quadratic form into another and show that when
%evaluated at the same point, they have the same result.
% A=[16,  4, -10;
%     4, 10, 16;
%    -10, 16, 4];
% bVec=[10;-4;6];
% c=2;
% [kVec,cNew]=quadEq2CenteredForm(A,bVec,c)
% xVec=2*randn(3,1)-1;
% z0=xVec'*A*xVec+bVec'*xVec+c
% z1=(xVec-kVec)'*A*(xVec-kVec)-cNew
%One will see that z0 and z1 are equal within finite precision limits.
%
%November 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

kVec=-(1/2)*(A\bVec);
cNew=bVec.'*inv(A)*bVec/4-c;
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
