function numUnique=numUniqueElsInArray(vecList,skipEl)
%%NUMUNIQUEELSINARRAY Count the number of unique elements in a vector.
%
%INPUTS: vecList A numElsXnumVecs list of vectors for which the number of
%                unique elements in each vector (not considering across
%                vectors) is desired. NaN values are not counted.
%         skipEl If there is an element that should not count toward being
%                unique (i.e. these are skipped), then it can be given. If
%                no elements should be skipped, then this input can be
%                omitted or an empty matrix can be passed.
%
%OUTPUTS: numUnique A numVecsX1 array, where the number of unique elements
%               of each (column) vector in the matrix vecList is given.
%
%This function sorts the array, which puts makes duplicates consecutive.
%then, it scans through the array, skipping duplicates. The sorting places
%NaN values, if any, at the highest end of the array, past Inf. NaN values
%are just skipped.
%
%EXAMPLE:
% vecList=[17, 24,  1,  8, 15;
%          17,  2,  7, 14, 16;
%           4,  6,  2, 20, 22;
%          10, 12, 19,  2,  3;
%          11, 18, 25,  2,  0];
% numUnique=numUniqueElsInArray(vecList,0)
%One should get
% numUnique=[4;5;5;4;4]
%
%September 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEls=size(vecList,1);
numVecs=size(vecList,2);

if(nargin>1&&~isempty(skipEl))
    vecList(vecList==skipEl)=NaN;
end

%Sorting the array makes all duplicates consecutive and places all NaN
%values at the end of the array (even after Inf).
vecList=sort(vecList,1,'ascend');

numUnique=zeros(numVecs,1);
for curVec=1:numVecs
    curEl=1;
    while(curEl<=numEls)
        %Advance the index past duplicates.
        while(curEl<numEls&&vecList(curEl,curVec)==vecList(curEl+1,curVec))
            curEl=curEl+1;
        end
        
        if(isnan(vecList(curEl,curVec)))
            break; 
        end
        numUnique(curVec)=numUnique(curVec)+1;
        curEl=curEl+1;
    end
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
