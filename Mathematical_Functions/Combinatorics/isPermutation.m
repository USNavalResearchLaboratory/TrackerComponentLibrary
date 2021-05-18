function boolVal=isPermutation(possiblePerm,startIdx)
%%ISPERMUTATION Given an array, determine whether it is a valid
%        permutation of the values startIdx:(startIdx+permLength-1).
%
%INPUTS: possiblePerm A length permLength vector to test for being a
%                     permutation.
%            startIdx The starting index of the permutation values to
%                     consider. If this input is omitted or an empty matrix
%                     is passed, then the default of 1 is used.
%
%OUTPUTS: boolVal This is true if possiblePerm is a permutation vector and
%                 false otherwise. Empty matrices are considered
%                 permutation vectors.
%
%After initial checks, the function just allocates a vector having the
%length of possiblePerm and adds 1 to each index given in possiblePerm. If
%any index is repeated. Afterwards, it checks whether all of the indices
%contain a "1".
%
%EXAMPLE:
% boolVal1=isPermutation([1;4;3;5])
% boolVal2=isPermutation([1;4;3;2;5])
% boolVal3=isPermutation([13;16;15;17;14],13)
%One will find boolVal1=false, boolVal2=true, boolVal3=true.
%
%September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(startIdx))
    startIdx=1;
end

possiblePerm=possiblePerm(:);

if(~isnumeric(possiblePerm)||~isreal(possiblePerm)||any(~isfinite(possiblePerm))||any(possiblePerm~=fix(possiblePerm)))
    boolVal=false;
    return;
end

numEls=length(possiblePerm);

possiblePerm=possiblePerm-startIdx+1;
if(any(possiblePerm>numEls)||any(possiblePerm<1))
    boolVal=false;
    return;
end

numValsPresent=zeros(numEls,1);
for k=1:numEls
    idx=possiblePerm(k);
    numValsPresent(idx)=numValsPresent(idx)+1;
end

boolVal=all(numValsPresent==1);

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
