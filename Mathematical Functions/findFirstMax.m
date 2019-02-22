function idx=findFirstMax(arr)
%%FINDFIRSTMAX Given an array of sorted vector of values in increasing
%              order, which might be full of duplicates, find the first
%              occurrence of the maximum value.
%
%INPUTS:    arr A vector of values, with possible repeats, sorted in
%               increasing order.
%
%OUTPUTS:   idx The index of the first occurrence of the maximum value
%               (last) element in arr.    
%
%The easiest way to find the first occurrence of the maximum value is to
%scan downward in the array from the end until a value different from the
%maximum value is found. However, that approach has a bad worst-case
%accuracy if all of the elements in the array have the same value. Thus,
%this method first checks for the easiest solutions, and then uses
%something akin to a binary search technique to find the first occurrence
%of the maximum value.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(arr);

%If there is only one element, then the problem is trivial.
if(N==1)
    idx=1;
    return
end

%Check for the simplest solution: The last element is the only maximum.
maxVal=arr(N);
if(arr(N-1)~=maxVal)
    idx=N;
    return
end

%In this case, we must perform a binary search for the second highest
%element.
maxIdx=N-1;
minIdx=1;

while(maxIdx~=minIdx&&maxIdx-minIdx>1)
    midIdx=floor((maxIdx+minIdx)/2);
    if(arr(midIdx)==maxVal)
       maxIdx=midIdx;
    else
       minIdx=midIdx;
    end
end

if(maxIdx==minIdx)
    idx=maxIdx;
else
    if(arr(minIdx)==maxVal)
        idx=minIdx;
    else
        idx=maxIdx;
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
