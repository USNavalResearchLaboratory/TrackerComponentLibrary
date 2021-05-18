function [thirdMomentMat,idxMat]=computeSampleThirdMoments(x,w,isNormalized,isCentral,algorithm)
%COMPUTESAMPLETHIRDMOMENTS Given possibly weighted samples of a vector
%       random variable, compute the (possibly normalized) third-order
%       moments. This can either be done as an nXn^2 matrix as defined
%       using an expected value and a Kronecker product as used in [1] and
%       [2], or one can just get the unique elements of the matrix.
%       Non-normalized, central third-order moments are given as the
%       expectation:
%       E{kron((x-E{x})*(x-E{x})',(x-E{x}))}
%       Non-normalized noncentral third-order moments are similarly:
%       E{kron(x*x',x)}
%       Normalization, as defined in Equation 30 of [2], involves division
%       by component variances, as described below.
%
%INPUTS: x The xDimXnumSamp set of samples. x and w together could also be,
%          for example, a set of cubature points.
%        w If the samples are weighted, then w is the length numSamp set of
%          all weights that sum to 1. Otherwise, if uniformly weighted, an
%          empty matrix can be passed.
% isNormalized If one desires normalized third moments, then this is true.
%          Otherwise, it is false. The default if omitted or an empty
%          matrix is passed is false.
% isCentral This is true if central moments are used. Otherwise, this is
%          false.
% algorithm An optional parameter indicating which algorithm to use.
%          Possible values are:
%          0 In this instance, thirdMomentMat is the xDimXxDim^2 matrix as
%            defined above, but it is computed in a more efficient manner
%            that using Kronecker products.
%          1 The xDimXxDim^2 matrix from the above definition leads to
%            having a number of repeated entries. Here, we only compute the
%            unique elements and put them in a single
%            (xDim^2+xDim*(xDim-1)*(xDim-2)/6)X1 vector. The output idxMat
%            is a 3X(xDim^2+xDim*(xDim-1)*(xDim-2)/6) matrix where
%            idxMat(:,i) holds the indices of the rows of x involved in
%            that particular third-order moment.
%          2 This is the same as 0, except kronecker products are dorectly
%            used as defined. This is a slow algorithm.
%
%OUTPUTS: thirdMomentMat The third-order moments. If algorithm=0 or
%            algorithm=2, this is an xDimXxDim^2 matrix. If algorithm=1,
%            this is a (xDim^2+xDim*(xDim-1)*(xDim-2)/6)X1 vector.
%     idxMat This 3X(xDim^2+xDim*(xDim-1)*(xDim-2)/6) matrix where
%            idxMat(:,i) holds the indices of the rows of x involved in
%            that particular third-order moment. See below for an
%            explanation of this indexation.
%
%Each of the elements in the expected value matrix
%E{kron((x-E{x})*(x-E{x})',(x-E{x}))} can be expressed as
%M_{i,j,k}=E{((x(i)-E{x(i)})*(x(j)-E{x(j)})*(x(k)-E{x(k)})}
%and similarly without the E{x} terms when using noncentral moments. The
%i,j, k indexation selects the dimensions of x involved in the moment. This
%is the same as the indexation returned by idxMat when algorithm=1. When
%algorithm=0, or algorithm=2, the elements of the xDimXxDim^2 matrix for
%thirdMomentMat can also be related back to these indices. For a row i1 and
%column i2, the corresponding i,j,k value is
%i=i1;
%k=fix((i2-1)/xDim)+1;
%j=i2-(k-1)*xDim;
%
%Normalization is most easily defined in terms of the M_{i,j,k} form of
%the elements. In such an instance, the normalized form is
%M^{norm}_{i,j,k}=M_{i,j,k}/sqrt(P(i,i)*P(j,j)*P(k,k))
%where P is the xDimXxDim covariance matrix of the distribution.
%he normalized third-order sample moments are also referred to as the 
%co-skewness in [2].
%
%EXAMPLE:
%In this example, we show that the different algorithms used here are
%consistent with each other. We generate 4D random samples of the
%exponential distirbutions and compute the normalized, central moments with
%each algorithm. We then show that the matrices from algorithms 0 and 3 are
%the same (within finite precision limits). Also, the values in the full
%matrix are the same as extracting just the unique values with algorithm 1.
%If one were to look at the moments themselves, one would find 4 moments
%near 2, and the others notably smaller (because they should all be 0-- the
%4 exponential distributions are independent).
% numPts=1e6;
% x=ExponentialD.rand([4,numPts],12);
% isNormalized=true;
% isCentral=true;
% tic
% y0=computeSampleThirdMoments(x,[],isNormalized,isCentral,0);
% toc
% tic
% y1=computeSampleThirdMoments(x,[],isNormalized,isCentral,1);
% toc
% tic
% y2=computeSampleThirdMoments(x,[],isNormalized,isCentral,2);
% toc
% 
% %This should be near 0 (difference just finite precision limits.
% max(abs(y0(:)-y2(:)))
% 
% %Show that the values in the full matrix are the same as the list giving
% %only the unique values.
% y0=sort(unique(y0(:)));
% y1=sort(y1(:));
% %This should be 0.
% max(abs(y1-y0))
%
%REFERENCES:
%[1] T. Kollo, "Multivariate skewness and kurtosis measures with
%    application in ICA," Journal of Multivariate Analysis, vol. 99, no.
%    10, pp. 2328-2338, Nov. 2008.
%[2] J. Dunik, O. Straka, B. Noack, J. Steinbring, and U. D. Hanebeck,
%    "Directional splitting of Gaussian density in non-linear random
%    variable transform," IET Signal Processing, vol. 12, no. 9, pp. 1073-
%    1081, Sep. 2018.
%
%February 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);
numPts=size(x,2);

if(nargin<2)
     w=[];
end

if(nargin<3||isempty(isNormalized))
    isNormalized=false;
end

if(nargin<4||isempty(isCentral))
    isCentral=true;
end

if(nargin<5||isempty(algorithm))
    algorithm=0; 
end

if(~isempty(w))
    w=w(:).';
        
    if(isNormalized)
        [xBar,PBar]=calcMixtureMoments(x,w);
        SP=sqrt(diag(PBar));
    else
        xBar=calcMixtureMoments(x,w);
    end

    if(isCentral)
        y=bsxfun(@minus,x,xBar);
    else
        y=x;
    end
else
    if(isNormalized)
        [xBar,PBar]=calcMixtureMoments(x);
        SP=sqrt(diag(PBar));
    else
        xBar=mean(x,2);
    end

    if(isCentral)
        y=bsxfun(@minus,x,xBar);
    else
        y=x;
    end
end

if(algorithm==0)
    %Compute the values without duplicates and put them into the proper
    %spots of the nXn^2 matrix. This computes it from the
    %element-by-element expression of Equation 29 of [2] and not from the
    %expressions directly using Kronecker products in [1] and [2].
    thirdMomentMat=zeros(numDim,numDim*numDim);
    if(~isempty(w))
        %Fill in the ones where all 3 are the different.
        for i1=1:(numDim-1)
            for i2=(i1+1):(numDim-1)
                for i3=(i2+1):numDim
                    momentVal=sum(bsxfun(@times,w,y(i1,:).*y(i2,:).*y(i3,:)));

                    if(isNormalized)
                        momentVal=momentVal/(SP(i1)*SP(i2)*SP(i3));
                    end

                    %For the (i1,i2,i3) order.
                    offset=i2+(i3-1)*numDim;
                    thirdMomentMat(i1,offset)=momentVal;
                    
                    %For the (i1,i3,i2) order.
                    offset=i3+(i2-1)*numDim;
                    thirdMomentMat(i1,offset)=momentVal;

                    %For the (i2,i1,i3) order.
                    offset=i1+(i3-1)*numDim;
                    thirdMomentMat(i2,offset)=momentVal;

                    %For the (i2,i3,i1) order.
                    offset=i3+(i1-1)*numDim;
                    thirdMomentMat(i2,offset)=momentVal;
                    
                    %For the (i3,i1,i2) order.
                    offset=i1+(i2-1)*numDim;
                    thirdMomentMat(i3,offset)=momentVal;
                    
                    %For the (i3,i2,i1) order.
                    offset=i2+(i1-1)*numDim;
                    thirdMomentMat(i3,offset)=momentVal;
                end
            end
        end
    
        %Fill in the elements where at least two are the same.
        for i1=1:numDim
            for i2=1:numDim
                momentVal=sum(bsxfun(@times,w,y(i1,:).*y(i1,:).*y(i2,:)));
                
                if(isNormalized)
                    momentVal=momentVal/(SP(i1)^2*SP(i2));
                end
                
                %For the (i1,i1,i2) order.
                offset=i1+(i2-1)*numDim;
                thirdMomentMat(i1,offset)=momentVal;
                
                %For the (i1,i2,i1) order.
                offset=i2+(i1-1)*numDim;
                thirdMomentMat(i1,offset)=momentVal;
                
                %For the (i2,i1,i1) order.
                offset=i1+(i1-1)*numDim;
                thirdMomentMat(i2,offset)=momentVal;
            end
        end
    else
        %Fill in the ones where all 3 are the different.
        for i1=1:(numDim-1)
            for i2=(i1+1):(numDim-1)
                for i3=(i2+1):numDim 
                    momentVal=mean(y(i1,:).*y(i2,:).*y(i3,:));
                    
                    if(isNormalized)
                        momentVal=momentVal/(SP(i1)*SP(i2)*SP(i3));
                    end
                    
                    %For the (i1,i2,i3) order.
                    offset=i2+(i3-1)*numDim;
                    thirdMomentMat(i1,offset)=momentVal;
                    
                    %For the (i1,i3,i2) order.
                    offset=i3+(i2-1)*numDim;
                    thirdMomentMat(i1,offset)=momentVal;

                    %For the (i2,i1,i3) order.
                    offset=i1+(i3-1)*numDim;
                    thirdMomentMat(i2,offset)=momentVal;

                    %For the (i2,i3,i1) order.
                    offset=i3+(i1-1)*numDim;
                    thirdMomentMat(i2,offset)=momentVal;
                    
                    %For the (i3,i1,i2) order.
                    offset=i1+(i2-1)*numDim;
                    thirdMomentMat(i3,offset)=momentVal;
                    
                    %For the (i3,i2,i1) order.
                    offset=i2+(i1-1)*numDim;
                    thirdMomentMat(i3,offset)=momentVal;
                end
            end
        end
    
        %Fill in the elements where at least two are the same.
        for i1=1:numDim
            for i2=1:numDim
                momentVal=mean(y(i1,:).*y(i1,:).*y(i2,:));
                
                if(isNormalized)
                    momentVal=momentVal/(SP(i1)^2*SP(i2));
                end
                
                %For the (i1,i1,i2) order.
                offset=i1+(i2-1)*numDim;
                thirdMomentMat(i1,offset)=momentVal;
                
                %For the (i1,i2,i1) order.
                offset=i2+(i1-1)*numDim;
                thirdMomentMat(i1,offset)=momentVal;
                
                %For the (i2,i1,i1) order.
                offset=i1+(i1-1)*numDim;
                thirdMomentMat(i2,offset)=momentVal;
            end
        end
    end
    
    idxMat=[];%Not used with this algorithm.
elseif(algorithm==1)
    %Just compute the values without duplicates.
    %The number of elements where all 3 components are different.
    numAll3Diff=numDim*(numDim-1)*(numDim-2)/6;
    
    %The number of elements where 2 components are the same or all 3 are
    %the same.
    numRest=numDim^2;
    thirdMomentMat=zeros(numAll3Diff+numRest,1);
    
    if(~isempty(w))
        %Fill in the ones where all 3 are the different.
        curIdx=1;
        for i1=1:(numDim-1)
            for i2=(i1+1):(numDim-1)
                for i3=(i2+1):numDim 
                    thirdMomentMat(curIdx)=sum(bsxfun(@times,w,y(i1,:).*y(i2,:).*y(i3,:)));
                    
                    if(isNormalized)
                        thirdMomentMat(curIdx)=thirdMomentMat(curIdx)/(SP(i1)*SP(i2)*SP(i3));
                    end
                    
                    curIdx=curIdx+1;
                end
            end
        end
    
        %Fill in the elements where at least two are the same.
        for i1=1:numDim
            for i2=1:numDim
                thirdMomentMat(curIdx)=sum(bsxfun(@times,w,y(i1,:).*y(i1,:).*y(i2,:)));
                
                if(isNormalized)
                    thirdMomentMat(curIdx)=thirdMomentMat(curIdx)/(SP(i1)^2*SP(i2));
                end
                
                curIdx=curIdx+1;
            end
        end
    else
        %Fill in the ones where all 3 are the different.
        curIdx=1;
        for i1=1:(numDim-1)
            for i2=(i1+1):(numDim-1)
                for i3=(i2+1):numDim 
                    thirdMomentMat(curIdx)=mean(y(i1,:).*y(i2,:).*y(i3,:));
                    
                    if(isNormalized)
                        thirdMomentMat(curIdx)=thirdMomentMat(curIdx)/(SP(i1)*SP(i2)*SP(i3));
                    end
                    
                    curIdx=curIdx+1;
                end
            end
        end
    
        %Fill in the elements where at least two are the same.
        for i1=1:numDim
            for i2=1:numDim
                thirdMomentMat(curIdx)=mean(y(i1,:).*y(i1,:).*y(i2,:));
                
                if(isNormalized)
                    thirdMomentMat(curIdx)=thirdMomentMat(curIdx)/(SP(i1)^2*SP(i2));
                end
                
                curIdx=curIdx+1;
            end
        end
    end
    
    if(nargout>1)
        %Fill in a vector of the indices of the values taken.
        idxMat=zeros(3,numAll3Diff+numRest);
        
        curIdx=1;
        for i1=1:(numDim-1)
            for i2=(i1+1):(numDim-1)
                for i3=(i2+1):numDim 
                    idxMat(:,curIdx)=[i1;i2;i3];
                    curIdx=curIdx+1;
                end
            end
        end
    
        %Fill in the elements where at least two are the same.
        for i1=1:numDim
            for i2=1:numDim
                idxMat(:,curIdx)=[i1;i1;i2];
                curIdx=curIdx+1;
            end
        end
    end
else
    %Compute it as an nXn^2 matrix using a kronecker delta as in Equation
    %28 in [2].
    
    thirdMomentMat=zeros(numDim,numDim*numDim);
    if(isempty(w))
        for k=1:numPts
            thirdMomentMat=thirdMomentMat+kron(y(:,k)*y(:,k)',y(:,k)');
        end
        thirdMomentMat=thirdMomentMat/numPts;
    else
        for k=1:numPts
            thirdMomentMat=thirdMomentMat+w(k)*kron(y(:,k)*y(:,k)',y(:,k)');
        end
    end
    
    if(isNormalized)
        thirdMomentMat=bsxfun(@rdivide,thirdMomentMat,SP);
        
        SP=SP';
        sel=1:numDim;
        for k=1:numDim
            thirdMomentMat(:,sel)=thirdMomentMat(:,sel)./(SP*SP(k));
            sel=sel+numDim;
        end
    end
    
    idxMat=[];%Not used with this algorithm.
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
