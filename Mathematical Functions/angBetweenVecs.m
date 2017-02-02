function angDiff=angBetweenVecs(v1,v2)
%%ANGBETWEENVECS  Given vector pairs in threespace, find the angular
%                 distances between them (non-negative values). A
%                 cross-product formula is used instead of the more logical
%                 dot-product formula so as to maximize numerical precision
%                 when vectors are nearly parallel. The vectors do not need
%                 to have the same magnitudes.
%
%INPUTS: v1       A 3XN matrix of N real vectors.
%        v2       A 3XN matrix of N real vectors.
%
%OUTPUTS: angDiff An NX1 matrix of the angular differences (in radians)
%                 between the vectors in v1 and the correspinding vectors
%                 in v2. The distances can range from 0 to pi.
%
%If one or both of the vectors has zero magnitude, then the angular
%difference will be zero. Note that angBetweenVecs(v1,v2) is the same as
%angBetweenVecs(v2,v1).
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The number of vectors.
N=size(v1,2);

angDiff=zeros(N,1);
for curVec=1:N
    angDiff(curVec)=atan2(norm(cross(v1(:,curVec),v2(:,curVec))),dot(v1(:,curVec),v2(:,curVec)));
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
