function c=cross7(a,b)
%%CROSS7 Find the cross product of two 7-dimensional vectors. The cross
%        product has the properties:
%        1) cross7(a,b) is a bilinear function of a and b.
%        2) cross7(a,b) is perpendicular to both a and b.
%           Thus, dot(cross7(a,b),a)=0 and dot(cross7(a,b),b)=0
%        3) norm(cross7(a,b))^2=norm(a)^2*norm(b)^2-dot(a,b)^2
%        Just as two forms of the 3D cross product exist (left-handed and
%        right-handed) 480 forms of the 7-dimensional cross product exist.
%        This function implements a cross product defined in terms of an
%        orthonormal basis using asymmetry as in [1].
%
%INPUTS: a A 7 by numVecs matrix of numVecs vectors.
%        b A 7 by numVecs matrix of numVecs vectors.
%
%OUTPUTS: c The 7 by numVecs cross products of a and b, aXb, for each of
%           the vectors in a and b.
%
%As noted in [2], it is only possible for a cross product between two
%vectors to have all of the properties listed above in three and seven
%dimensions. 
%
%As noted in 1, the Jacobi identity does not hold. This means that
%cross7(cross7(a,b),c)+cross7(cross7(b,c),a)+cross7(cross7(c,a),b) does
%not necessarily equal zero. However, the following identities can be
%proven:
%cross7(a,b)=-cross7(b,a)
%dot(a,cross7(b,c))=dot(b,cross7(c,a))=dot(c,cross7(a,b))
%cross7(a,cross7(a,b))=dot(a,b)*a-norm(a)^2*b
%cross7(cross7(a,b),cross7(a,c))=cross7(cross7(cross7(a,b),c),a)+cross7(cross7(cross7(b,c),a),a)+cross7(cross7(cross7(c,a),a),b)
%
%REFERENCES:
%[1] P. Lounesto, "Octonians and triality," Advances in Clifford Algebras,
%    vol. 11, no. 2, pp. 191-213, Dec. 2001.
%[2] W. S. Massey, "Cross products of vectors in higher dimensional
%    Euclidean spaces," The American Mathematical Monthly, vol. 90, no. 10,
%    pp. 697-701, Dec. 1983.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVecs=size(a,2);
c=zeros(7,numVecs);

c(1,:)=+0              +a(2,:).*b(4,:)  +a(3,:).*b(7,:)  -a(4,:).*b(2,:)  +a(5,:).*b(6,:)  -a(6,:).*b(5,:)  -a(7,:).*b(3,:);
c(2,:)=-a(1,:).*b(4,:) +0               +a(3,:).*b(5,:)  +a(4,:).*b(1,:)  -a(5,:).*b(3,:)  +a(6,:).*b(7,:)  -a(7,:).*b(6,:);
c(3,:)=-a(1,:).*b(7,:) -a(2,:).*b(5,:)  +0               +a(4,:).*b(6,:)  +a(5,:).*b(2,:)  -a(6,:).*b(4,:)  +a(7,:).*b(1,:);
c(4,:)=+a(1,:).*b(2,:) -a(2,:).*b(1,:)  -a(3,:).*b(6,:)  +0               +a(5,:).*b(7,:)  +a(6,:).*b(3,:)  -a(7,:).*b(5,:);
c(5,:)=-a(1,:).*b(6,:) +a(2,:).*b(3,:)  -a(3,:).*b(2,:)  -a(4,:).*b(7,:)  +0               +a(6,:).*b(1,:)  +a(7,:).*b(4,:);
c(6,:)=+a(1,:).*b(5,:) -a(2,:).*b(7,:)  +a(3,:).*b(4,:)  -a(4,:).*b(3,:)  -a(5,:).*b(1,:)  +0               +a(7,:).*b(2,:);
c(7,:)=+a(1,:).*b(3,:) +a(2,:).*b(6,:)  -a(3,:).*b(1,:)  +a(4,:).*b(5,:)  -a(5,:).*b(4,:)  -a(6,:).*b(2,:)  +0;

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
