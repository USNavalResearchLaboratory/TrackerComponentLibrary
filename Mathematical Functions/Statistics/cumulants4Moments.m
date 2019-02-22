function cumMat=cumulants4Moments(momentMat)
%%CUMULANTS4MOMENTS Given a matrix of non-central univariate or
%           multivariate moments, compute the cumulants corresponding to
%           the moments.
%
%INPUTS: momentMat A matrix taking n indices, where 
%               momentMat(a1,a2,a3...an) corresponds to the coefficient of 
%               a E(x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1)) moment.
%               The order of the moment is a1+a2+..an. It is assumed that
%               moments only up to an order equal to the size of the first
%               dimension of the matrix are available (other entries in
%               momentMat are ignored). All other dimensions must be at
%               least the same size as the first. Also, the zero-th orde
%               moment must be 1 (as the PDF integrates to 1). If a
%               univariate distribution is used, then this is a column
%               vector.
%
%OUTPUTS: cumMat A matrix taking n indices, where cumMat(a1,a2,a3...an)
%                corresponds to the cumulant having the orders
%                corresponding to momentMat(a1,a2,a3...an).
%
%Cumulants can be useful in some applications, such as when using an
%Edgeworth series to approximate a PDF. Also cumulants have some nice
%properties not posessed by moments. For example, the cumulants of the sum
%of two independent random variables is equal to the sum of the cumulants
%of the variables. That is not the case for moments in general. Cumulants
%are often defined in terms of the cumulant generating function, which is
%the natural logarithm of the moment generating function, natural logarithm
%of the moment
%
%Note that obtaining cumulants from sample moments using this function is
%probably biased. However, using an algorithm such as [2] to obtain an
%unbiased result is extremely slow for large number of samples and requires
%so many terms that finite precision errors tend to add up quickly. Thus,
%the approach of obtaining cumulants from sample moments is generally
%preferable.
%
%The recursive formulation used in Equations 10 and 12 of [1] is
%implemented here. Equation 6 and 5 are used to initialize the necessary
%initial batch of cumulants and moments based on the negative of those
%cumulants that are needed for Equation 12 and 10.
%
%REFERENCES:
%[1] P. J. Smith, "A recursive formulation of the old problem of obtaining
%    moments from cumulants and vice versa," The American Statistician,
%    vol. 49, no. 2, pp. 217-218, May 1995.
%[2] I. J. Good, "A new formula for k-statistics," The Annals of
%    Statistics, vol. 5, no. 1, pp. 224-228, Jan. 1977.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDimList=size(momentMat);
maxDeg=numDimList(1)-1;
numIdx=length(numDimList);

%If the last dimension is just a singleton dimension, then shrink
%everything by one. This will speed things up, because Equations 10 and 12
%will not be needed.
if(numIdx==2&&numDimList(2)==1)
    numIdx=1;
end

%Allocate space to hold a table of binomial values that can be quickly
%looked up in the loops rather than having to compute them again and again.
binomTable=makeBinomTable(maxDeg);

%Allocate space for the cumulants and moments that correspond to the
%NEGATIVE of the cumulants (cumulants for moments in Equation 11).
cumMat=zeros(numDimList);
negCumMomMat=zeros(numDimList);

%The zeroth-order moment (from negative cumulants) is always 1.
negCumMomMat(1)=1;
%The zeroth-order cumulant is always zero.
cumMat(1)=0;

%All of the cumulants involving just the first index being nonzero are
%computed using Equation 6. The moments corresponding to the negative
%cumulants are filled in using Equation 5. All of the moments of just the
%first element are the first maxDeg+1 entries in momentMat, so no special
%multi-dimensional indexation is necessary here.
for n=1:maxDeg
    %Equation 6
    cumMat(n+1)=momentMat(n+1);
    for i=1:(n-1)
        cumMat(n+1)=cumMat(n+1)-binomTable(n-1+1,i+1)*cumMat(n-i+1)*momentMat(i+1);
    end
    
    %We also need the moments corresponding to the negative of the
    %cumulants. This is Equation 5.
    for i=0:(n-1)
        negCumMomMat(n+1)=negCumMomMat(n+1)+binomTable(n-1+1,i+1)*(-cumMat(n-i+1))*negCumMomMat(i+1);
    end
end

%We have now computed all of the cumulants for the case where only the
%first index is nonzero. We now must use the recursion of Equation 10 and
%12 to fill in the other cumulants, each time adding a new variable with a
%nonzero power.
for numNonzeroIdx=2:numIdx
    %For a given new index that will be nonzero, we have to go through all
    %of the moments based on negative cumulants and through all of the
    %cumulants for different degrees of moments.
    %We start with the last index being nonzero, because we have
    %already considered all of the instances where it is zero.
    for nonzeroDegVal=1:maxDeg
        %curDeg is the degree of all elements except the last one, which is
        %nonzero and which we go through in the outer loop.
        for curDeg=0:(maxDeg-nonzeroDegVal)
            %Fixing the order of highest index, we must now go through all
            %partitions of possible values that the other elements can
            %take.
            curPart=getNextMPartition(numNonzeroIdx-1+curDeg,numNonzeroIdx-1);
            while(~isempty(curPart))
                %curPart-1 is the (unordered) degrees of the indices before
                %the current index (the last nonzero index).
                %For the given partition, we must go through all possible
                %multiset permutations of the values.
                lowerIdxPerms=genAllMultisetPermutations(curPart);
                numPerms=size(lowerIdxPerms,2);
                
                for curPerm=1:numPerms
                    r=[lowerIdxPerms(:,curPerm)-1;nonzeroDegVal];
                    
                    idxVal=nDim2Index(numDimList(1:numNonzeroIdx),r+1);

                    %Use Equation 12 to update the cumulant terms.
                    cumMat(idxVal)=equation1012Update(r,momentMat,negCumMomMat,1,binomTable);
                    
                    %Use Equation 10 to update the moment terms that
                    %utilize the negative of the cumulants.
                    negCumMomMat(idxVal)=equation1012Update(r,cumMat,negCumMomMat,-1,binomTable);
                end
                curPart=getNextMPartition(curPart);
            end
        end
    end
end

end

function sumVal=equation1012Update(r,mat1,mat2,signMat1,binomTable)
%Equations 10 and 12 in [1] have the same form. The only difference is
%in the terms used in the innermost sums. In Equation 10, the first
%innermost term is cumulants, in 12, it is moments. In 10, the second term
%is moments, in 12, it is a special set of moments based on negative
%cumulants.
%
%We want to evaluate Equation 10 to find the special set of moments based
%on negative cumulants, so that we can then use Equation 12 to find the
%cumulants. When using equation 10 to get moments based on negative
%cumulants, we need to flip the sign of the cumulants. Thus, the coeffMat1
%allows a scalar coefficient to be supplied to multiply the values from
%mat1 that are used.

numIdx=length(r);
maxDim=size(mat2);
maxDim=maxDim(1:numIdx);%Only the first numIdx indices are considered.

%This stores cumulative values of the product of the binomials going down
%each sum. This way, the values do not have to be constantly recomputed.
binomProdVal=ones(numIdx-1,1);

%This vector holds all of the i values. Initially, all sums start at zero.
iVec=zeros(numIdx,1);

%The loop goes through all of the values in the sums.
curLevel=numIdx;
isAscending=false;
sumVal=0;
while(curLevel>=1)
    if(curLevel==numIdx)
        %Go through the innermost sum
        for i=0:(r(numIdx)-1)
            binomVal=binomProdVal(numIdx-1)*binomTable(r(numIdx)-1+1,i+1);
            
            iVec(numIdx)=i;
            kIdxVec=r-iVec;
            
            idx=nDim2Index(maxDim,iVec+1);
            kidx=nDim2Index(maxDim,kIdxVec+1);
            
            sumVal=sumVal+binomVal*(signMat1*mat1(kidx))*mat2(idx);
        end
        curLevel=curLevel-1;
        isAscending=true;
        continue;
    elseif(isAscending)
        %Try incrementing the order at this level...
        iVec(curLevel)=iVec(curLevel)+1;
        if(iVec(curLevel)<=r(curLevel))
            %If the value is valid, then update the binomial product value
            %and descend to the next level.
            if(curLevel~=1)
                binomProdVal(curLevel)=binomProdVal(curLevel-1)*binomTable(r(curLevel)+1,iVec(curLevel)+1);
            else
                binomProdVal(curLevel)=binomTable(r(curLevel)+1,iVec(curLevel)+1);
            end
            isAscending=false;
            curLevel=curLevel+1;
            continue;
        else
            %If the value is invalid, then just keep ascending.
            curLevel=curLevel-1;
            continue;
        end
    else%We are descending in the sums here.
        iVec(curLevel)=0;
        binomProdVal(curLevel)=binomProdVal(curLevel-1);
        curLevel=curLevel+1;
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
