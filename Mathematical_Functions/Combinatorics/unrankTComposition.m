function p=unrankTComposition(theRank,t,n,firstElMostSig)
%%UNRANKTCOMPOSITION Return the composition of the given theRank, starting
%                   from zero, of n unlabeled items in t parts. That is, a
%                   method of putting n unlabeled balls into t labeled
%                   slots. Compositions are tuples of t integers >=1 that
%                   sum to n. Compositions are given in lexicographic order
%                   where the first element is the least significant unless
%                   otherwise specified.
%
%INPUTS: theRank The theRank (position in an ordered sequence) of the
%                desired composition starting from zero. There is a total
%                of binomial(n-1,t-1) compositions.
%              t The number of slots that can hold items, >=1.
%              n The number of items that are composed into slots, >=1;
%                n>=t
% firstElMostSig An optional parameter specifying whether the first
%                element is the most or least significant. The default if
%                omitted or an empty matrix is passed is false (the first
%                element is the least significant).
%
%OUTPUTS: p A tX1 vector holding the current composition, whose elements
%           sum to n. Each element is the number of "balls" in that slot.
%           If theRank is larger than the maximum number of compositions,
%           an empty matrix is returned. The values of q can range from 0
%           to n.
%
%This is an implementation of the relation described in Chapter 7.2.1.3 of
%[1] for mapping combinations to compositions. theRank is used to obtain a
%combination in lexicographic order (p(1) is the least significant
%element) and the elements of the combination are mapped to the equivalent
%composition. However, a special case is if the input is one-dimensional,
%in which case the only composition is for theRank=0 and the value is n.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(firstElMostSig))
    firstElMostSig=false;
end

n=n-1;
nCombo=n;
tCombo=t-1;

if(tCombo>0)
    c=unrankCombination(theRank,nCombo,tCombo);

    if(isempty(c))
        p=[];
        return;
    end

    %Transform the combination into a valid composition.
    p=zeros(t,1);
    p(1)=c(1)+1;
    for curIdx=2:tCombo
        p(curIdx)=c(curIdx)-c(curIdx-1);
    end
    p(t)=n-c(tCombo);
else%The t=1 case has to be handled separately.
    if(theRank==0)
        p=n+1;
    else
        p=[]; 
    end
end

if(firstElMostSig)
   p=flipud(p(:)); 
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
