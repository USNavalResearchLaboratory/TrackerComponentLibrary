function val=calcMSE(xTrue,xEst,type,is3D)
%%CALCMSE Compute the mean-squared error (MSE) of estimates compared to
%         true values. This can be given either as a vector or as a MSE
%         matrix, which is akin to a covariance matrix.
%
%INPUTS: xTrue The truth data. This is either an xDimXNumSamples matrix or
%              an xDimXNXnumSamples matrix. The latter formulation is
%              useful when the MSE over multiple Monte Carlo runs of an
%              entire length-N track is desired. In the first formulation,
%              we take N=1. Alternatively, if the same true value is used
%              for all numSamples, then xTrue can just be an xDimXN matrix.
%              Alternatively, if xTrue is the same for all numSamples and
%              all N, then just an xDimX1 matrix can be passed. N and
%              numSamples are inferred from xEst.
%         xEst An xDimXnumSamples set of estimates or an xDimXNXnumSamples
%              set of estimates (if values at N times are desired).
%         type An optional parameter specifying the type of output to
%              provide. Possible values are
%              0 (The default if omitted or an empty matrix is passed)
%                Return a scalar MSE value.
%              1 Return an xDimX1 MSE vector.
%              2 Return an xDimXxDim MSE matrix. Unlike just the vector,
%                this includes cross terms.
%         is3D An optional indicating that xEst is 3D. This is only used if
%              xEst is a matrix. In such an instance, there is an ambiguity
%              whether xEst is truly 2D, so N=1, or whether numSamples=1
%              and xEst is 3D with the third dimension being 1. If this
%              parameter is omitted or an empty matrix is passed, it is
%              assumed that N=1.
%
%OUTPUTS: val The MSE either given as a 1XN vector of scalar MSEs
%             (algorithm 0), an xDimXN set of MSE vectors, or an
%             xDimXxDimXN collection of MSE matrices
%
%The matrix form can reveal information on the correlation of the errors
%just as the cross terms in a covariance matrix do.
%
%EXAMPLE 1:
%Here, we show that the MSE of a Gaussian random variable is related to the
%covariance matrix.
% R=[28,   4, 10;
%     4,  22, 16;
%     10, 16, 16];%The covariance matrix.
% xTrue=[10;-20;30];
% numRuns=100000;
% xEst=GaussianD.rand(numRuns,xTrue,R);
% valScal=calcMSE(xTrue,xEst,0)
% valVec=calcMSE(xTrue,xEst,1)
% valMat=calcMSE(xTrue,xEst,2)
%One will see that the scalar value is close to the trace of R, the vector
%values are close to the diagonal values of R and the matrix values are
%close to R as a whole.
%
%EXAMPLE 2:
%This is the same as example 1, but here, we show an example where N is not
%1 and where the mean differs for the different dimensions.
% R=[28,   4, 10;
%     4,  22, 16;
%     10, 16, 16];%The covariance matrix.
% xTrue=[[10;-20;30],[20;-8;186]];
% N=size(xTrue,2);
% xDim=size(xTrue,1);
% numRuns=100000;
% xEst=zeros(xDim,N,numRuns);
% for k=1:N
%     xEst(:,k,:)=reshape(GaussianD.rand(numRuns,xTrue(:,k),R),[xDim,1,numRuns]);
% end
% valScal=calcMSE(xTrue,xEst,0)
% valVec=calcMSE(xTrue,xEst,1)
% valMat=calcMSE(xTrue,xEst,2)
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(is3D))
    is3D=false; 
end

if(nargin<3||isempty(type))
    type=0; 
end

xDim=size(xEst,1);
if(ismatrix(xEst)&&is3D==false)
    N=1;
    numSamples=size(xEst,2);
    xEst=reshape(xEst,[xDim,1,numSamples]);
    
    if(size(xTrue,2)==1)
        %If the true values are the same for all samples.
        xTrue=repmat(xTrue,[1,1,numSamples]);
    else
        xTrue=reshape(xTrue,[xDim,1,numSamples]);
    end
else
    N=size(xEst,2);
    numSamples=size(xEst,3);
    
    if(ismatrix(xTrue))
        if(size(xTrue,2)==1)
            %If the true values are the same for all samples and for all N.
            xTrue=repmat(xTrue,[1,N,numSamples]);
        else        
            %If the true values are the same for all samples.
            xTrue=repmat(xTrue,[1,1,numSamples]);
        end
    end
end

switch(type)
    case 0%Scalar
        val=zeros(1,N);
        for k=1:N
            val(k)=sum(sum((xEst(:,k,:)-xTrue(:,k,:)).^2,3))/numSamples;
        end
    case 1%Vector
        val=zeros(xDim,N);
        for k=1:N
            val(:,k)=sum((xEst(:,k,:)-xTrue(:,k,:)).^2,3)/numSamples;
        end
    case 2%Matrix
        val=zeros(xDim,xDim,N);
        
        for k=1:N
            diff=reshape(xEst(:,k,:)-xTrue(:,k,:),[xDim,numSamples]);
            val(:,:,k)=diff*diff'/numSamples;
        end
    otherwise
        error('Unknown type specified.')
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
