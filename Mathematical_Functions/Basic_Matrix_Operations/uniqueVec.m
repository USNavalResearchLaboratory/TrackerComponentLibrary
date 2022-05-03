function vUnique=uniqueVec(v)
%%UNIQUEVEC This function takes an mXn matrix of n m-dimensional vectors
%           and returns only the vectors that are unique. This is currently
%           suited for a small to moderate  numbers of vectors, because a
%           brute-force comparison algorithm is used (rather than a
%           sorting-based approach that might lessen the number of
%           comparisons).
%
%INPUTS: v An mXn matrix of n vectors.
%
%OUTPUTS: vUnique An mXnumUnique matrix of the unique vectors in v.
%
%EXAMPLE:
% v=[1, 1, 1, 1, 1, 12;
%    1, 2, 1, 2, 2, 24;
%    1, 3, 1, 4, 3, 36]
% uniqueVec(v)
%
%February 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(v))
    vUnique=v;
    return;
end

numRows=size(v,1);
numCols=size(v,2);

vUnique=zeros(numRows,numCols);

vUnique(:,1)=v(:,1);
numUnique=1;
for k=2:numCols
    isUnique=true;
    for kPrev=1:numUnique
        if(all(v(:,k)==vUnique(:,kPrev)))
            isUnique=false;
            break;
        end
    end

    if(isUnique)
        numUnique=numUnique+1;
        vUnique(:,numUnique)=v(:,k);
    end
end

%Size to fit
vUnique=vUnique(:,1:numUnique);

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
