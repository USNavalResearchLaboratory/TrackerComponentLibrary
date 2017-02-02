function quasiMomentMat=moments2QuasiMoments(momentMat,muN,SigmaN)
%%MOMENTS2QUASIMOMENTS Convert a set of (multivariate) noncentral moments
%               to quasi- moments with respect to a particular mean vector
%               and covariance matrix. Quasi-moments are used in
%               generalized Edgeworth series approximations of
%               distributions.
%
%INPUTS: momentMat A matrix taking n indices, where 
%               momentMat(a1,a2,a3...an) corresponds to the coefficient of 
%               a E(x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1)) moment.
%               The order of the moment is a1+a2+..an. It is assumed that
%               moments only up to an order equal to the size of the first
%               dimension of the matrix are available (other entries in
%               momentMat are ignored). All other dimensions must be at
%               least the same size as the first. Also, the zero-th order
%               moment must be 1 (as the PDF integrates to 1). If a
%               univariate distribution is used, then this is a column
%               vector.
%           muN The nX1 mean vector with respect to with the quasi-moments
%               are computed. If this parameter is omitted or an empty
%               matrix is passed, a zero mean vector is used.
%        SigmaN The nXn covariance matrix with respect to which the
%               quasi-moments are computed. If this parameter is omitted or
%               an empty matrix is passed, the identity matrix is used.
%
%OUTPUTS: quasiMomentMat A matrix taking n indices, where 
%               quasiMomentMat(a1,a2,a3...an) corresponds to the quasi-
%               moment whose multivariate order is given by a1-1,a2-1, etc.
%
%Quasi-moments are introduced in [2] and are defined in terms of a certain
%mean vector and covariance matrix. If mean vector and covariance matrix
%coincide with those of the moments associated with the quasi-moments, then
%the quasi-moments are called proper quasi-moments.
%
%The conversion is based on Equation 2.11 in [1]. The equation is
%rearranged to get
%q_\alpha=m_\alpha-\sum_{\beta}^{\alpha-1}C_{\alpha}^\beta m^N_{\alpha-\beta} q_\beta
%where q's represent quasi moments, m^N's represent moments of the Normal
%distribution having covariance matrix Sigma, m's are moments of the
%distribution that are to be converted. The indexation is vector based.
%Thus, beta is actually a vector of beta(1)...beta(n) and alpha goes from
%alpha(1)...alpha(n). The sum is actually multiple sums. The alpha-1
%subtracts 1 from the last nonzero index in the alpha vector. The
%C_{\alpha}^\beta} is equal to the product of binomial(alpha(i),beta(i))
%over i.
%
%REFERENCES:
%[1] I. E. Poloskov, "CAS Mathematica in random studies," in Proceedings
%    of Computational Science - ICCS 2003, Melbourne, Australia and St.
%    Petersburg, Russia, 2-4 Jun. 2003, pp. 781-790.
%[2] P. I. Kuznetsov, R. L. Stratonovich, and V. I. Tikhonov, "Quasi-moment
%    functions in the theory of random processes," Theory of Probability 
%    and its Applications, vol. V, no. 1, pp. 80-97, 1960.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDimList=size(momentMat);
maxDeg=numDimList(1)-1;
numIdx=length(numDimList);

%If the last dimension is just a singleton dimension, then shrink
%everything by one.
if(numIdx==2&&numDimList(2)==1&&(nargin<3||size(SigmaN,2)==1))
    numIdx=1;
end

%Do the conversion with respect to the standard normal distribution if no
%further details are given.
if(nargin<2||isempty(muN))
    muN=zeros(numIdx,1);
end
if(nargin<3||isempty(SigmaN))
    SigmaN=eye(numIdx,numIdx);
end

%Allocate space to hold a table of binomial values that can be quickly
%looked up in the loops rather than having to compute them again and again.
binomTable=makeBinomTable(maxDeg);

%This is a table that will hold values of moments for the normal
%distribution with zero mean and covariance matrix Sigma. These moments
%play a role in the conversion to quasi-moments as they are related to
%hermite polynomials.
momentMatN=zeros(numDimList);

momentMatN(1)=1;%Zeroth-order moment is always 1.
for curDeg=1:maxDeg
    curPart=getNextMPartition(numIdx+curDeg,numIdx);
    
    while(~isempty(curPart))
        %Go through all permutation of the orders of each index.
        powTable=genAllMultisetPermutations(curPart-1);
        numVals=size(powTable,2);
        for curVal=1:numVals
            numDerivs=powTable(:,curVal);
            momentVal=GaussianD.momentGenFun(muN,SigmaN,numDerivs);
            curEl=nDim2Index(numDimList,numDerivs+1);
            momentMatN(curEl)=momentVal;
        end
        
        curPart=getNextMPartition(curPart);
    end
end

%Allocate space for the return values.
quasiMomentMat=zeros(numDimList);
%The zeroth-order moment is always 1.
quasiMomentMat(1)=1;

%First, fill in the values for the case where only the first dimension has
%nonzero exponents. This is just the inverse of Equation 11 for the scalar
%case.
for n=1:maxDeg
    sumVal=0;
    for i=0:(n-1)
        sumVal=sumVal+binomTable(n+1,i+1)*momentMatN(n-i+1)*quasiMomentMat(i+1);
    end
    quasiMomentMat(n+1)=momentMat(n+1)-sumVal;
end

for numNonzeroIdx=2:numIdx
    %For a given new index that will be nonzero. 
    %We start with the last index (not visited) being nonzero, because we
    %have already considered all of the instances where it is zero.
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
                    
                    %Use the inverse of Equation 2.11 to compute the moment
                    %term.
                    quasiMomentMat(idxVal)=equation211InvUpdate(r,momentMat(idxVal),momentMatN,quasiMomentMat,binomTable);
                end
                curPart=getNextMPartition(curPart);
            end
        end
    end
end
end

function sumVal=equation211InvUpdate(alpha,momentVal,momentMatN,quasiMomentMat,binomTable)
%This implements the inverse of Equation 2.11 in [1] for a given alpha
%vector. That is, find q not m.

numIdx=length(alpha);
maxDim=size(quasiMomentMat);
maxDim=maxDim(1:numIdx);%Only the first numIdx indices are considered.

%This stores cumulative values of the product of the binomials going down
%each sum. This way, the values do not have to be constantly recomputed.
binomProdVal=ones(numIdx-1,1);

%This vector holds all of the i values. Initially, all sums start at zero.
iVec=zeros(numIdx,1);

%The loop goes through all of the values in the sums.
curLevel=numIdx;
isAscending=false;
sumVal=momentVal;
while(curLevel>=1)
    if(curLevel==numIdx)
        %Go through the innermost sum
        for i=0:alpha(numIdx)
            binomVal=binomProdVal(numIdx-1)*binomTable(alpha(numIdx)+1,i+1);
            
            iVec(numIdx)=i;
            kIdxVec=alpha-iVec;
            
            idx=nDim2Index(maxDim,iVec+1);
            kidx=nDim2Index(maxDim,kIdxVec+1);
            
            %The last deltaTerm will be subtracted out before returning,
            %because the sum does not include the vary last term, which is
            %for this quasi moment.
            deltaTerm=binomVal*momentMatN(kidx)*quasiMomentMat(idx);
            sumVal=sumVal-deltaTerm;
        end
        curLevel=curLevel-1;
        isAscending=true;
        continue;
    elseif(isAscending)
        %Try incrementing the order at this level...
        iVec(curLevel)=iVec(curLevel)+1;
        if(iVec(curLevel)<=alpha(curLevel))
            %If the value is valid, then update the binomial product value
            %and descend to the next level.
            if(curLevel~=1)
                binomProdVal(curLevel)=binomProdVal(curLevel-1)*binomTable(alpha(curLevel)+1,iVec(curLevel)+1);
            else
                binomProdVal(curLevel)=binomTable(alpha(curLevel)+1,iVec(curLevel)+1);
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

%Subtract out the las term.
sumVal=sumVal-deltaTerm;

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
