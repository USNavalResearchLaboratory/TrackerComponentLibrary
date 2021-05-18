function [combosVec,combosInPlace]=genAllMultisetCombos(m,k)
%%GENALLMULTISETCOMBOS Generate all multiset length-k combinations of n
%           unique elements, each which might be repeated in the total set
%           of elements. The output is given in lexicographic order.
%
%INPUTS: m An nX1 vector where each spot represents one unique item in the
%          multiset. m(i) is the number of copies of element i in the
%          multiset (the multiplicity of the ith element). m(i)>=1 for all
%          i.
%        k The number of items to choose from multiset i.
%
%OUTPUS: combosVec An nXnumCombos matrix such that combosVec(m,i) holds the
%                  number of times element m is repeated in combination i.
%                  Combinations are given in lexicographic order.
%    combosInPlace A kXnumCombos matrix such that combosInPlace(m,i) is the
%                  element in spot m of combination i. Numbering of the
%                  elements starts at 1. Combinations are given in
%                  lexicographic order.
%
%The algorithm is just the first procedure mset in Section 2 of [1]. The
%algorithm based on the twisted lexico tree is not used. Unlike the
%presentation in [1], the recursion has been eliminated.
%
%EXAMPLE 1:
%An example is the example given in Section 2 of [1]. Here:
% m=[1;2;2;1;1];
% k=4;
% [combosVec,combosInPlace]=genAllMultisetCombos(m,k);
%One will get 
% combosVec=[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
%             0, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2;
%             2, 1, 2, 2, 0, 1, 1, 2, 1, 2, 2, 0, 1, 1, 2, 0, 0, 1;
%             1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0;
%             1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0];
% combosInPlace=[3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
%                3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2;
%                4, 4, 3, 3, 4, 3, 3, 3, 4, 3, 3, 4, 3, 3, 3, 2, 2, 2;
%                5, 5, 5, 4, 5, 5, 4, 3, 5, 5, 4, 5, 5, 4, 3, 5, 4, 3];
%
%EXAMPLE 2:
%One thing of note is that increasing the number of repetitions for any
%element beyong the number of slots (k) does not change the result one the
%number of repetitions exceeds k. For example.
% m=[2;2;1;1];
% k=2;
% [combosVec0,combosInPlace0]=genAllMultisetCombos(m,k)
% %produces the same results as
% m=[200;200;1;1];
% k=2;
% [combosVec1,combosInPlace1]=genAllMultisetCombos(m,k)
%
%REFERENCES:
%[1] T. Takaoka, "O(1) time generation of adjacent multiset combinations,"
%    arXiv, 28 Feb. 2015. [Online].
%    Available: http://arxiv.org/abs/1503.00067
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(k==0)
    combosVec=[];
    combosInPlace=[];
    return;
end

n=length(m);
a=zeros(n,1);

numCombos=numMultisetCombos(m,k);

combosVec=zeros(n,numCombos);

b(n+1+1)=0;
for i=n:-1:1
   b(i+1)=b(i+1+1)+m(i); 
end

i=1;
kList=zeros(n+1,1);
kList(1)=k;
upperList=zeros(n,1);
jList=zeros(n,1);

goingDown=true;
curCombo=1;
while(i>0)
    if(goingDown)
        if(i<=n)
            k=kList(i);

            lower=max([k-b(i+1+1),0]);
            upperList(i)=min([m(i),k]);
            jList(i)=lower;
            a(i)=jList(i);
            i=i+1;
            kList(i)=k-jList(i-1);
        else%Save the combination and go back up.
           combosVec(:,curCombo)=a;
           curCombo=curCombo+1;
           i=i-1;
           goingDown=false;
        end
    else%If ascending
        jList(i)=jList(i)+1;
        if(jList(i)<=upperList(i))
            k=kList(i);
            a(i)=jList(i);
            
            i=i+1;
            kList(i)=k-jList(i-1);
            goingDown=true;
        else%Keep ascending
            i=i-1;
        end
    end
end

if(nargout>1)
    combosInPlace=zeros(k,numCombos);
    
    for curCombo=1:numCombos
        setIdx=find(combosVec(:,curCombo));
        numIdx=length(setIdx);
        
        curPlace=1;
        for curIdx=1:numIdx
            idx=setIdx(curIdx);
            numRepeats=combosVec(idx,curCombo);
            
            for curRepeat=1:numRepeats
                combosInPlace(curPlace,curCombo)=idx;
                curPlace=curPlace+1;
            end
        end
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
