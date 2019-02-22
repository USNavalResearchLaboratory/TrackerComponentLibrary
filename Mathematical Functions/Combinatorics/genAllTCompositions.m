function tuples=genAllTCompositions(n,t,algorithm)
%%genAllTGENALLTCOMPOSITIONS Find all permutations of t positive (nonzero)
%                   integers that sum to a given constant n (constant-sum
%                   t-tuples). That is, find all compositions of n into t
%                   parts such that no part is ever empty. There is a total
%                   of binomial(n-1,t-1) compositions.
%
%INPUTS: n The integer value that the sum of the elements of the vectors
%          should have.
%        t The integer dimensionality of the vectors of integers to
%          generate.
% algorithm An optional parameter specifying which algorithm to use to
%          generate the compositions. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            relation described in Section 7.2.1.3 of [1] for mapping
%            combinations to compositions.
%          1 Use  Algorithm 8.1.1 of [1], which is also given in Appendix C
%            of [2].
%
%OUTPUTS: tuples A tXnumTuples matrix of column vectors such that the sum
%                across rows for each column equals n. All elements are
%                >=1. If n<t, then an empty matrix is returned. There is a
%                total of numTuples=binomial(n-1,t-1);
%
%Often, one might prefer all partitions of t integers that sum to a
%constant value n INCLUDING zeros. This function can be used to
%produce such a modified set of partitions as
%tuples=genAllTCompositions(n+t,t)-1;
%
%To preallocate the array, one has to know the total number of tuples that
%sum to a particular value. To determine the total number of possible
%tuples:
%Consider a series of n ones. One has to place t-1 commas between
%them. These t-1 commas partition the set into t subsets. The sums of the t
%subsets are the values of an element in a tuple. This means that there are
%consequently binomial(n-1,t-1) possible tuples.
%
%One can observe that for t=3, starting with n=3 and going up, the number
%of tuples is given by triangular numbers. For t=4, increasing n from 4,
%one gets tetrahedral numbers; for t=5, increasing n from 5, one gets
%pentatope numbers. One can see the pattern of figurate numbers arising
%from the number of possible tuples that sum to a given value.
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%[2] T. Gerster, "Sparse grid quadrature methods for computational
%    finance," 2007, Habilitationsschrift, Rheinischen Friedrich-Wilhelms-
%    Universität Bonn (Germany).
%    [Online]. Available: http://wissrech.ins.uni-bonn.de/research/pub/gerstner/gerstner_habil.pdf
%[3] V. Kaarnioja, "Smolyak quadrature," Master's thesis, University of
%    Helsinki, Department of Mathematics and Statistics, Helsinki, Finland,
%    May 2013.
%    [Online]. Available: https://helda.helsinki.fi/bitstream/handle/10138/40159/thesis.pdf
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<t)
    tuples=[];
    return;
end

if(nargin<3||isempty(algorithm))
    algorithm=0; 
end

switch(algorithm)
    case 0
        tuples=genAllTCompositions0(n,t);
    case 1
        tuples=genAllTCompositions1(n,t);
    otherwise
        error('Unknown algorithm specified.')
end
end

function compList=genAllTCompositions0(n,t)
%%GENALLTCOMPOSITIONS Generate all compositions of n unlabeled items in t
%                   parts. That is, a method of putting n unlabeled balls
%                   into t labeled slots. Compositions are tuples of n
%                   integers >=1 that sum to t. 
%
%INPUTS: n The number of items that are composed into slots, >=1; n>=t.
%        t The number of slots that can hold items, >=1.
%              
%OUTPUTS: compList A tXnumComp matrix of the compositions.
%
%This is an implementation of the relation described in Section 7.2.1.3 of
%[1] for mapping combinations to compositions. The total number of
%compositions is binomial(n-1,t-1).
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numCompositions=binomial(n-1,t-1);

n=n-1;
nCombo=n;
tCombo=t-1;

if(tCombo>0)
    compList=zeros(t,numCompositions);
    
    %The first combination that will be mapped to a composition.
    c=(0:(tCombo-1)).';
    curComp=1;
    
    p=zeros(t,1);
    while(~isempty(c))
        %Transform each combination into a valid composition.
        
        p(1)=c(1)+1;
        for curIdx=2:tCombo
            p(curIdx)=c(curIdx)-c(curIdx-1);
        end
        p(t)=n-c(tCombo);

        compList(:,curComp)=p;
        
        c=getNextCombo(c,nCombo,0);
        curComp=curComp+1;
    end
    
else
    compList=n+1; 
end

end

function tuples=genAllTCompositions1(sumVal,d)
%%GENALLTCOMPOSITIONS1 Generate all length d compositions of integers that
%                     sum up to sumVal using a method based on Algorithm
%                     8.1.1 of [1]
%
%The algorithm is based on Algorithm 8.1.1 of [1]. It is also given in
%Appendix C of [2]. Unlike in the references, the special case of d=sumVal
%is treated here.
%
%REFERENCES
%[1] T. Gerster, "Sparse grid quadrature methods for computational
%    finance," 2007, Habilitationsschrift, Rheinischen Friedrich-Wilhelms-
%    Universität Bonn (Germany).
%    [Online]. Available: http://wissrech.ins.uni-bonn.de/research/pub/gerstner/gerstner_habil.pdf
%[2] V. Kaarnioja, "Smolyak quadrature," Master's thesis, University of
%    Helsinki, Department of Mathematics and Statistics, Helsinki, Finland,
%    May 2013.
%    [Online]. Available: https://helda.helsinki.fi/bitstream/handle/10138/40159/thesis.pdf
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Handle the special case of d==sumVal, which is not handled by the
%algorithm.
if(d==sumVal)
   tuples=ones(d,1);
   return;
elseif(d>sumVal)
    tuples=[];
    return;
end

k=ones(d,1);
kHat=(sumVal-d+1)*k;

numTuples=binomial(sumVal-1,d-1);
tuples=zeros(d,numTuples);

q=1;
curTuple=1;
while(k(d)<=sumVal)
    k(q)=k(q)+1;
    if(k(q)>kHat(q))
        if(q~=d)
            k(q)=1;
            q=q+1;
        end
    else
        kHat(1:(q-1))=kHat(q)-k(q)+1;
        k(1)=kHat(1);
        q=1;
        tuples(:,curTuple)=k;
        curTuple=curTuple+1;
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
