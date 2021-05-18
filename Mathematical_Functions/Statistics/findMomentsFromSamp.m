function momentMat=findMomentsFromSamp(xSamp,wSamp,maxDeg,maxDegType,normalize)
%%FINDMOMENTSFROMSAMP Given a set of samples of a univariate or
%               multivariate random variable, determine all of the
%               (non-central) sample moments up to a certain maximum degree
%               for each element of the multivariate vector. To only find a
%               single moment, use the function findMomentFromSamp.
%
%INPUTS: xSamp An nXnumSamp matrix of numSamp samples from which moments
%              are to be computed.
%        wSamp If the samples are weighted, then wSamp is a 1XnumSamp or
%              numSampX1 vector of the weights. If the samples are not
%              weighted, then an empty matrix can be passed.
%       maxDeg A value specifying the nature of the maximum power of each
%              element in the vector to be computed.  A moment has the form
%              E(x1^a1*x2^a1...*xn^pn) where E is the expected
%              value operator and the sum of p1+p2+...pn is the order
%              of the moment. If maxDegType is 0 (the default), then maxDeg
%              is the maximum order of the moment computed. If
%              maxDegType=1, then maxDegType is the maximum power of each
%              component of x for the moments, meaning that the maximum
%              moment computed is order n*maxDeg.
%   maxDegType An optional parameter specifying how many moments are
%              computed. Possible values are
%              0 (The default if omitted or an empty matrix is passed)
%                maxDeg specifies the maximum order of the moments.
%              1 maxDeg specifies the maximum degree of each component of x
%                going into the moments. This is slower than maxDegType 0.
%    normalize An optional boolean value indicating whether the weights in
%              wSamp should be normalized prior to use. The default if
%              this parameter is omitted is false. The weights cannot be
%              normalized if they sum to zero.
%
%OUTPUTS: momentMat A matrix taking n indices, where 
%               momentMat(a1,a2,a3...an) corresponds to the coefficient of 
%               a E(x1^(a1-1)*x2^(a2-1)*x3^(a3-1)...xn^(an-1)) moment.
%               The order of the moment is a1+a2+..an. If maxDegType=0 or
%               is omitted, then moments with orders higher than maxDeg are
%               set to zero.
%
%Noncentral moments are unbiased if one just sums the values of the powers
%of the components of each sample and divides by the number of samples. For
%weighted samples, we multiply each product in the sum by the sample
%weight, and then divide by the sum of the weights.
%
%Note, however, that central samples computed in such a manner ARE biased.
%Hence the reason introductory textbooks use a 1/(numSamp-1) in front of
%variance computations. The bias in central moments (when just dividing by
%numSamp) comes from using the sample mean in place of the true mean to
%centralize the moment.
%
%EXAMPLE:
%This function can be used to verify the accuracy of cubature formula. For
%example, consider 4D, fifth-order cubature points for a Normal 0-I
%distribution:
% [xi, w]=fifthOrderCubPoints(4);
% momentMat=findMomentsFromSamp(xi,w,5)
%One can verify that the zeroth-order moment is 1, the second order
%moments: 
%momentMat(2+1,0+1,0+1,0+1),momentMat(0+1,2+1,0+1,0+1),
%momentMat(0+1,0+1,2+1,0+1), and momentMat(0+1,0+1,0+1,2+1) are all ones
%with the other second-order moments zeros and the fourth-order xi^4
%moments are similarly all 3's with all others numerically zero. Other
%nonzero terms are less than eps(1) and are thus effectively zero.
%Consequently, the function fifthOrderCubPoints correctly matches moments 0
%through 5.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(normalize))
    normalize=false; 
end

if(nargin<4||isempty(maxDegType))
    maxDegType=0;
end

n=size(xSamp,1);
numSamp=size(xSamp,2);

if(isempty(wSamp))
    wSamp=ones(numSamp,1)/numSamp;
end
if(normalize)
    normConst=sum(wSamp);
else
    normConst=1;
end

%Allocate space for the return values; zero the components.
numMaxDims=repmat(maxDeg+1,[1,n]);
totalNumEls=prod(numMaxDims);
if(isscalar(numMaxDims))
    %This prevents Matlab from automatically allocationg a
    %numMaxDimsXnumMaxDims matrix in zeros if numMaxDims is scalar.
    numMaxDims=[numMaxDims,1];
end
momentMat=zeros(numMaxDims);

%Add in the effects of each sample.
for k=1:numSamp
    curSamp=xSamp(:,k);
    
    %We will precompute all of the necessary powers of curSamp.
    sampPow=zeros(n,maxDeg+1);
    sampPow(:,0+1)=ones(n,1);%The zeroth power is always 1.
    for curPow=1:maxDeg
        sampPow(:,curPow+1)=sampPow(:,curPow-1+1).*curSamp;
    end
    
    if(maxDegType==1)
        %Next, we will go through all of the possible permutations of
        %powers of the elements.
        for curEl=1:totalNumEls
            powIdx=index2NDim(numMaxDims,curEl);
            prodVal=1;
            for curDim=1:n
                prodVal=prodVal*sampPow(curDim,powIdx(curDim));
            end
            momentMat(curEl)=momentMat(curEl)+prodVal*wSamp(k);
        end
    else
        %Next, we will go through all of the partitions of the elements for
        %each order and then all of the permutations per patition.
        for curDeg=0:maxDeg
            pRecur=getNextMPartition(n+curDeg,n);

            while(~isempty(pRecur))
                %Get all permutations of pRecur;
                permMat=genAllMultisetPermutations(pRecur);
                numPerm=size(permMat,2);
                for curPerm=1:numPerm
                    prodVal=1;
                    for curDim=1:n
                        prodVal=prodVal*sampPow(curDim,permMat(curDim,curPerm));
                    end
                    idxVal=nDim2Index(numMaxDims,permMat(:,curPerm));

                    momentMat(idxVal)=momentMat(idxVal)+prodVal*wSamp(k);
                end

                pRecur=getNextMPartition(pRecur);
            end
        end
    end
end
momentMat=momentMat/normConst;
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
