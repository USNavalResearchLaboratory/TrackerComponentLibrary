function prodVal=numCapacitatedCombos(numPerBin,numBins,totalItems)
%%NUMCAPACITATEDCOMBOS Given totalItems items, find the number of
%               combinations of placing the items into numBins unlabeled
%               bins where each bin can hold numPerBin items. If
%               numPerBin=1, then this is just counting standard
%               combinations of numBins items taken from a total of
%               totalItems. This is the number of ways of putting
%               totalItems labeled balls into numBins unlabeled bins such
%               that each bin holds numPerBin items and the ordering of the
%               items in the bins does not matter.
%
%INPUTS: numPerBin The scalar integer capacity of the bins; numPerBin>=1.
%          numBins The scalar integer number of bins; numBins>=1.
%       totalItems The scalar integer total number of items to place in the
%                  bins; totalitems>=numBins*numPerBin. If omitted or an
%                  empty matrix is passed, then the default of 
%                  totalItems=numBins*numPerBin is used.
%
%OUTPUTS: prodVal The number of ways of placing totalItems items into
%                 numBins unlabeled bins each holding numperBin items,
%                 where the ordering of the items in the bins does not
%                 matter.
%
%If we assume that totalItems=n=numPerBin*numBins, then we can count by
%fixing the first element of each bin to the smallest index unassigned
%to any previous bin in the current combination being constructed. Thus,
%the first bin will always contain item 1. The second bin will always
%contain the smallest unassigned item that is not in the first bin
%and so on. Thus, we have to choose
%binomial(n-(k-1)*numPerBin-1,numPerBin-1) for the kth bin conditioned on
%the previous bins. So, the total number of ways of distributing the items
%across the bins is the product of those binoial terms from k=1 to
%numBins-1 (the last bin only has one way of inserting things, because we
%assume that (totalItems=numPerBin*numBins).
%
%For the case where totalItems>numPerBin*numBins, so totalItems~=n, we take
%the previously found quantity and multiply it by binomial(totalItems,n)
%which is the number of ways of choosing which subset of the totalItems
%items will be put into the bins.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=numPerBin*numBins;

if(nargin<3||isempty(totalItems))
    totalItems=n;
end

if(n>totalItems)
    error('totalItems cannot be less than numBins*numPerBin')
end

%Count how n items can be spread across the bins, given that the bins are
%unlabeled and the ordering in the bins does not matter.
prodVal=1;
for k=1:(numBins-1)
    prodVal=prodVal*binomial(n-(k-1)*numPerBin-1,numPerBin-1);
end

%Now, choose which subset of the n items goes into the bins.
prodVal=binomial(totalItems,n)*prodVal;

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
