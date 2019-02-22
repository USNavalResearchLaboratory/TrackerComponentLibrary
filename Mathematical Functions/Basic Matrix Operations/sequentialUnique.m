function [vals,numReps]=sequentialUnique(x)
%%SEQUENTIALUNIQUE Like the unique function, this function removes
%           duplicate values in a vector. However, this function only
%           groups together sequential values that are identical. Thus, if
%           two identical numbers are separated by a different number, then
%           they will not be grouped. This function also returns the number
%           of repeats of each of the values that we combined into a single
%           value.
%
%INPUTS: x A vector of matrix. In matrices, the elements are taken in the
%          order that linear indexation would take them. That is x(i) gives
%          the ith element.
%
%OUTPUTS: vals A numElX1 linear vector giving the sequentially unique
%              values in x. Repetition is determined by exact equality.
%      numReps A numElX1 vector listing how many times each element is
%              repeated.
%
%EXAMPLE:
% x=[2;1;1;1;2;3;4;5;5;6;6;7;7;1;1];
% [vals,numReps]=sequentialUnique(x)
% %One gets vals=[2;1;2;3;4;5;6;7;1];
% %and numReps=[1;3;1;1;1;2;2;2;2];
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numReps=diff(find(diff([-Inf;x(:);Inf])));
vals=x(cumsum(numReps));

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
