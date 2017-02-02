function A=randBinMat(dims,numOnes)
%%RANDBINMAT Create a random binary matrix. The matrix can either be
%            created with each entry having a 50% probability of being 0 or
%            1, or it can be created with a fixed number of ones in random
%            positions.
%
%INPUTS: dims A vector specifying the dimensions of the desired matrix,
%             the same as what the size function would return. For example,
%             for a 6X8 matrix, dims=[6,8];
%     numOnes An optional scalar specifying the number of ones that should
%             be in the random matrix. If this parameter is omitted, then
%             each entry in the matrix has a 50% probability of being one.
%
%The algorithms are straightforward. If numOnes is provided, then a random
%subset of numOnes of the elements is set to 1. Otherwise, random variables
%are generated for all of the elements as either 1 or 0.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A=zeros(dims(:)');
numElA=numel(A);

if(nargin<2)
%If the number of ones is not specified, then each entry has a 50%
%probability of being 1.
    A(:)=(rand(numElA,1)>0.5);
else
%If the number of ones is specified, then set them randomly.
    A(randperm(numElA,numOnes))=1;
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
