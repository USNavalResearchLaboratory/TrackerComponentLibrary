function tuples=getAllCompositions(d,sumVal)
%%GETALLCOMPOSITIONS Find all permutations of d positive (nonzero) integers
%                   that sum to a given constant sumVal (constant-sum
%                   d-tuples). That is, find all compositions of sumVal
%                   into d parts such that no part is ever empty. There is
%                   a total of binomial(sumVal-1,d-1) compositions.
%
%INPUTS:      d The integer dimensionality of the vectors of integers to
%               generate.
%        sumVal The integer value that the sum of the elements of the
%               vectors should have.
%
%OUTPUTS: tuples A dXnumTuples matrix of column vectors such that the sum
%                across rows for each column equals d. All elements are
%                >=1. If sumVal<d, then an empty matrix is returned. The
%                tuples are not in lexicographic order. There is a total of
%                numTuples=binomial(sumVal-1,d-1);
%
%Often, one might prefer all partitions of d integers that sum to a
%constant value sumVal INCLUDING zeros. This function can be used to
%produce such a modified set of partitions as
%tuples=getAllCompositions(d,sumVal+d)-1;
%
%The algorithm is based on Algorithm 8.1.1 of [1]. It is also given in
%Appendix C of [2]. Unlike in the references, the special cases of d=sumVal
%and d<sumVal are treated here. Also, the array of tuples to be returned is
%preallocated so that it does not change size each loop.
%
%To preallocate the array, one has to know the total number of tuples that
%sum to a particular value. To determine the total number of possible
%tuples:
%Consider a series of sumVal ones. One has to place d-1 commas between
%them. These d-1 commas partition the set into d subsets. The sums of the d
%subsets are the values of an element in a tuple. This means that there are
%consequently binomial(sumVal-1,d-1) possible tuples.
%
%One can observe that for d=3, starting with sumVal=3 and going up, the
%number of tuples is given by triangular numbers. For d=4, increasing
%sumVal from 4, one gets tetrahedral numbers; for d=5, increasing sumVal
%from 5, one gets pentatope numbers. One can see the pattern of figurate
%numbers arising from the number of possible tuples that sum to a given
%value.
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
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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
