function theRank=rankTComposition(p,firstElMostSig)
%%RANKTCOMPOSITION Obtain the lexicographic order (the rank) of a particular
%                 composition starting at 0. A composition of n into t
%                 is a way of putting n unlabeled balls into t labeled
%                 slots. Compositions are tuples of t integers >0 that
%                 sum to n. Compositions are given in lexicographic order
%                 where the first element is the least significant unless
%                 otherwise specified.
%
%INPUTS: p An tX1 vector holding the current composition, whose elements
%          sum to n>=1.
% firstElMostSig An optional parameter specifying whether the first element
%          is the most or least significant. The default if omitted or an
%          empty matrix is passed is false (the first element is the least
%          significant).
%
%OUTPUTS: theRank The rank of the composition in a lexicographic ordering
%                 of compositions, counting from zero. There is a total of
%                 binomial(n-1,t-1) compositions. Lexicographic ordering is
%                 defined using p(1) as the least significant element.
%
%This is an implementation of the relation described in Chapter 7.2.1.3 of
%[1] for mapping combinations to compositions. However, a special case is
%if the input is one-dimensional. In such an instance, the rank is always
%0.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Special case, empty matrix.
if(isempty(p))
   theRank=[];
   return;
end

%Special case, given a scalar.
if(length(p)==1)
    theRank=0;
   return; 
end

if(nargin<2||isempty(firstElMostSig))
    firstElMostSig=false;
end

if(firstElMostSig)
   p=flipud(p(:)); 
end

m=length(p);

%Turn the composition into the equivalent combination.
mCombo=m-1;

c=zeros(mCombo,1);
c(1)=p(1)-1;

for k=2:mCombo
    c(k)=p(k)+c(k-1); 
end

%get the rank of the combination.
theRank=rankCombination(c);
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
