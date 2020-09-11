function [w,mu,P]=RunnalsGaussMixRed(w,mu,P,K)
%%RUNNALSGAUSSMIXRED Perform Gaussian mixture reduction using the greedy
%                    merging algorithm by Runnals in [1].
%
%INPUTS: w An NX1 or 1XN vector of weights of the components of the
%          original Gaussian mixture.
%       mu An xDimXN matrix of the means of the vector components of the
%          original Gaussian mixture.
%        P An xDim XxDim XN hypermatrix of the covariance matrices for the
%          components of the original Gaussian mixture.
%        K The number of components desired in the mixture after reduction.
%
%OUTPUTS: w The KX1 weights of the mixture after reduction.
%        mu The xDimXK means of the mixture after reduction.
%         P The xDimXxDimXK covariance matrices of the mixture after
%           reduction.
%
%This implements the suboptimal Gaussian mixture reduction algorithm of
%[1]. Note that a real efficient C/C++ implementation would be implemented
%using various data structures to avoid the cost deletion and inseration of
%matrix elements that is used in the reduction algorithm.
%
%REFERENCES:
%[1] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Transactions on Aerospace and Electronic Systems,
%    vol. 43, no. 3, pp. 989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    N=length(w);
    
    %If no reduction is necessary.
    if(N<=K)
        return;
    end

    %We will only be using one triangular of this matrix.
    M=Inf*ones(N,N);%This is the cost matrix.
    
    %We shall fill the cost matrix with the cost of all pairs.
    for cur1=1:(N-1)
        for cur2=(cur1+1):N
            M(cur1,cur2)=BDist(w(cur1),w(cur2),mu(:,cur1),mu(:,cur2),P(:,:,cur1),P(:,:,cur2));
        end
    end
        
    Nr=N;
    for mergeRound=1:(N-K)
        [Mm,minRows]=min(M);%Minimize over the rows
        [~,minCol]=min(Mm);%Minimize over the columns.
        minRow=minRows(minCol);
        
        %Now we know which two hypotheses to merge.
        curClust=[minRow minCol];
        [w,mu,P]=mergeGaussianComp(w,mu,P,curClust);
        
        %Now we must remove the two hypotheses from the cost matrix and add
        %the merged hypothesis to the end of the matrix.
        if(minCol>minRow)
            M(minCol,:)=[];
            M(minRow,:)=[];%Delete the rows and columns.
            M(:,minCol)=[];
            M(:,minRow)=[];
        else
            M(minRow,:)=[];%Delete the rows and columns.
            M(minCol,:)=[];
            M(:,minRow)=[];
            M(:,minCol)=[];
        end
        Nr=Nr-1;
        
        %Now we must add the distances for the new stuff to the end of the
        %matrix.
        M=[M Inf*ones(Nr-1,1);ones(1,Nr)*Inf];
        
        for curRow=1:(Nr-1)
            M(curRow,Nr)=BDist(w(curRow),w(Nr),mu(:,curRow),mu(:,Nr),P(:,:,curRow),P(:,:,Nr));
        end
    end
end

function val=BDist(w1,w2,mu1,mu2,P1,P2)
%%BDIST This is the distance measure given in Equation 21 in Section VI-B
%       in [1].
%
%REFERENCES:
%[1] A. R. Runnalls, "Kullback-Leibler approach to Gaussian mixture
%    reduction," IEEE Trans. Aerosp. Electron. Syst., vol. 43, no. 3, pp.
%    989-999, Jul. 2007.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    diff=mu1-mu2;
    
    w1m=w1/(w1+w2);
    w2m=w2/(w1+w2);
    P12=w1m*P1+w2m*P2+w1m*w2m*(diff*diff');
    
    val=0.5*((w1+w2)*log(det(P12))-w1*log(det(P1))-w2*log(det(P2)));
    
    %Deal with the case where w1 and w2 are both essentially vero.
    if(~isfinite(val))
        val=0;
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
