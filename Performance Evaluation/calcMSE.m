function val=calcMSE(xTrue,xEst,type)
%%CALCMSE Compute the mean-squared error (MSE) of estimates compared to
%         true values. This can be given either as a vector or as a MSE
%         matrix, which is akin to a covariance matrix.
%
%INPUTS: xTrue The truth data. If the truth is the same for all samples of
%              the estimate, then this is an xDimX1 vector. Otherwise, this
%              is an xDimXN matrix where N is the number of samples in the
%              batch.
%         xEst An xDim X N set of estimates.
%         type An optional parameter specifying the type of output to
%              provide. Possible values are
%              0 (The default if omitted or an empty matrix is passed)
%                Return a scalar MSE value.
%              1 Return an xDimX1 MSE vector.
%              2 Return an xDimXxDim MSE matrix. Unlike just the vector,
%                this includes cross terms.
%
%OUTPUTS: val The MSE either given as a scalar, an xDimX1 vector or an
%             xDimXxDim matrix depending on the type input.
%
%The matrix form can reveal information on the correlation of the errors
%just as the cross terms in a covariance matrix do.
%
%EXAMPLE:
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
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(type))
   type=0; 
end

numEst=size(xEst,2);

if(size(xTrue,2)==1)
    xTrue=repmat(xTrue,[1,numEst]);
end

switch(type)
    case 0%Scalar
        val=sum(sum((xEst-xTrue).^2,2)/numEst);
    case 1%Vector
        val=sum((xEst-xTrue).^2,2)/numEst;
    case 2%Matrix
        diff=(xEst-xTrue);
        val=diff*diff'/numEst;
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
