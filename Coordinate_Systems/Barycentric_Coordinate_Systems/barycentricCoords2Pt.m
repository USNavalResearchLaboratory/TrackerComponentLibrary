function x=barycentricCoords2Pt(phi,v)
%%BARYCENTRICCOORDS2PT Given a set of weights in a barycentric coordinate
%       system, obtain the point to which those weights correspond. That
%       is, just compute the weighted sum of the vertices in v. This
%       function is a convenient way of transforming many points at once.
%
%INPUTS: phi A numPtsXN matrix of N sets of weights for the vertices.
%          v A numDimXnumPts set of the vertices. Each vertices is a
%            numDimX1 vector.
%
%OUTPUTS: x The numDimXN set of synthesized points.
%
%A number of different barymetric coordinate systems exist. For example,
%consider those in [1].
%
%REFERENCES:
%[1] M. S. Floater, "Generalized barycentric coordinates and applications,"
%    Acta Numerica, vol. 24, pp. 161-214, 1 May 2015.
%
%September 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(v,1);
numPts=size(phi,2);

x=zeros(numDim,numPts);
for k=1:numPts
    x(:,k)=sum(bsxfun(@times,v,phi(:,k).'),2);
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
