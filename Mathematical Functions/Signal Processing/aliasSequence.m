function seqWrap=aliasSequence(seq,period)
%%ALIASSEQUENCE Consider a sequence of length N that should be put into a
%           vector of length period. This function puts the values of the
%           original sequence into an output vector, aliasing values in
%           elements beyond length period back around to the beginning of
%           the output sequence and summing them.
%
%INPUTS: seq An NXnumSeq set of numSeq vectors to alias.
%     period The length of the vector into which seq is aliased
%            (period>=1).
%
%OUTPUTS: seqWrap The NXnumSeq set of aliased sequences. If period>N, then
%                 the final elements are zeros. 
%
%The need for aliasing vectors of data can arise when performing signal
%processing. Issues related to aliasing and the sampling of signals are
%described in Chapter 4.2 of [1]. Aliasing also plays a role in circular
%convolutions, as described in Section 5.4.2 of [1].
%
%EXAMPLE:
%Consider the sequence 1:10 that we want to fit into a spot of length 3.
%The aliasing means that the resulting sequence is the sum of the
%following sequences:
%1  2 3
%4  5 6
%7  8 9
%10 0 0 %Note the zeros padded to the end.
%_______
%Result =22 15 18
% seq=(1:10)';
% result=aliasSequence(seq,3)
%Executing that provides the expected result.
%
%REFERENCES:
%[1] S. K. Mitra, Digital Signal Processing: A Computer-Based Approach,
%    3rd ed. Boston: McGraw Hill, 2006.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(seq,1);
numSeq=size(seq,2);

seqWrap=zeros(period,numSeq);

for k=0:(N-1)
    idx=mod(k,period)+1;
    seqWrap(idx,:)=seqWrap(idx,:)+seq(k+1,:);
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
