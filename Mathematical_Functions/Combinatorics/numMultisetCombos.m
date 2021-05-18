function numCombos=numMultisetCombos(m,k)
%%NUMMULTISETCOMBOS Determine the number of multiset combinations for a
%             given multiset choosing k items from the set. A multiset is a
%             set of elements where some of the elements are repeated.
%
%INPUTS: m An nX1 vector where each spot represents one unique item in the
%          multiset. m(i) is the number of copies of element i in the
%          multiset (the multiplicity of the ith element). m(i)>=1 for all
%          i.
%        k The number of items to choose from multiset i.
%
%OUTPUTS: numCombos The integer number of possible multiset combinations.
%
%The algorithm is based on the inclusion-exclusion principle, wbhich is
%discussed in Section 2 of [1], among many other places.
%
%REFERENCES:
%[1] T. Takaoka, "O(1) time generation of adjacent multiset combinations,"
%    arXiv, 28 Feb. 2015. [Online].
%    Available: http://arxiv.org/abs/1503.00067
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(m);

numCombos=binomial(n+k-1,k);

signVal=-1;
for i=1:floor((k+1)/2)
    curCombo=0:(i-1);
    
    sumVal=0;
    while(~isempty(curCombo))
        kCur=k-sum(m(curCombo+1))-i;

        if(kCur>=0)
            sumVal=sumVal+binomial(n+kCur-1,kCur);
        end
        
        curCombo=getNextCombo(curCombo,n);
    end

    numCombos=numCombos+signVal*sumVal;
    
    signVal=-signVal;
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
