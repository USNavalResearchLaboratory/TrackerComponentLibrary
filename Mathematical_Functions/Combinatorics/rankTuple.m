function rank=rankTuple(tuple,maxVals)
%%RANKTUPLE Obtain the order of a particular tuple in the total sequence of
%           tuples given a certain set of maximum values for each digit.
%
%INPUTS: tuple An NX1 or 1XN tuple. The values of the digits in the tuple
%              range from 0 to the value given in maxVals.
%      maxVals An NX1 or 1XN vector of the maximum value that each digit of
%              the tuple can take.
%
%OUTPUTS: rank The rank of the tuple (the order in the sequence of
%              increasing tuples where the first digit is the most
%              significant) starting from 0.
%
%The nDim2Index function actually ranks tuples, but it assumes that the
%digits of the tuples start from 1, not zero. This function just calls
%nDim2Index with the appropriate parameters.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

rank=nDim2Index(maxVals+1,tuple+1)-1;

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
