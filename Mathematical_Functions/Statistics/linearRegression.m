function [ab,PInv]=linearRegression(t,z,W)
%%LINEARREGRESSION Given a dependent variable t and a set of points z,
%           which can be multidimensional, find the parameters of a least-
%           squares fit of a line fo the form z=a*t+b, where z,a, and b can
%           be vectors. a and b are returned as [a,b] and thus polyValVec
%           can be used to evaluate arbitrary points on the line. This
%           function minimizes the quadratic cost (with optional weight W)
%           sum_{i=1}^numPoints (z(:,i)-a*t(i)-b)'*W(:,:,k)*(z(:,i)-a*t(i)-b)
%           
%
%INPUTS: t A length numPoints vector holding the (assumed error free) times
%          of the observations.
%        z A zDimXnumPoints matrix of numPoints observations associated
%          with the times in t.
%        W If W is omitted or an empty matrix is passed, then all of the
%          points and dimensions are not weighted (assigned unity weight).
%          If a length numPoints vector is provided, then the cost function
%          being minimized is
%          sum_{i=1}^numPoints (z(:,i)-a*t(i)-b)'*W(k)*(z(:,i)-a*t(i)-b)
%          If an zDimXzDimXnumPoints hypermatrix is provided, then the cost
%          function 
%          sum_{i=1}^numPoints (z(:,i)-a*t(i)-b)'*W(:,:,k)*(z(:,i)-a*t(i)-b)
%          is optimized. Optimal weights are the inverse covariances of the
%          measurements.
%
%OUTPUTS: ab This is zDimX2 vector where ab(:,1)=a and ab(:,2)=b for the
%            equation of a line. ab can be passed directly to the
%            polyValVec function.
%
%Determining an optimal least squares fit of points to a line is a special
%case of the batch least squares estimation of Chapter 3.4 of [1]. The
%expressions in there are valid for the case where there is no weighting of
%the points as well as assuming that optimal weighting is given, in which
%case W would be inverse covariance matrices of the measurements.
%
%EXAMPLE 1:
%Here, we demonstrate that in the noiseless case, the unweighted solution
%produces the exact result.
% a=[1;2;-3;4];
% b=[-8;0;12;6];
% N=20;
% t=10*randn(1,N);
% z=bsxfun(@plus,bsxfun(@times,a,t),b);
% abFit=linearRegression(t,z);
% AbsErrAB=abFit-[a,b]
%one will see that the absolute errors in the fitted a and b are on the
%order of what one might expect fro finite precision limitations.
%
%EXAMPLE 2:
%Here, we have noisy measurements with one having a very high noise. The
%noise variance is the same for all of the dimensions at one time. We
%show that including the inverse variances as weights produces a better fit
%than exclusing them. Also, we show that the NEES is around 1, which
%indicates that the covariance matrix provided is consistent.
% numMCRuns=1000;
% a=[1;2;-3;4];
% b=[-8;0;12;6];
% numDim=length(a);
% N=100;
% 
% sigma=0.05*ones(N,1);
% sigma(N)=500000;
% w=1./sigma.^2;
% 
% RMSEUnweighted=0;
% RMSEWeighted=0;
% NEES=0;
% for curRun=1:numMCRuns
%     t=10*randn(1,N);
%     zTrue=polyValVec([a,b],t);
%     %Add noise to the measurements
%     z=zeros(numDim,N);
%     for k=1:N
%         z(:,k)=zTrue(:,k)+sigma(k)*randn(numDim,1);
%     end
% 
%     %Fit the lines.
%     ab0=linearRegression(t,z);
%     [ab1,PInv1]=linearRegression(t,z,w);
%     %Evaluate at the fitting points.
%     z0=polyValVec(ab0,t);
%     z1=polyValVec(ab1,t);
%     %Compute RMSE
%     RMSEUnweighted=RMSEUnweighted+calcMSE(z0,zTrue,false);
%     RMSEWeighted=RMSEWeighted+calcMSE(z1,zTrue,false);
%     %Compute NEES
%     diff=ab1(:)-ab(:);
%     NEES=NEES+diff'*PInv1*diff;
% end
% RMSEUnweighted=sqrt(RMSEUnweighted/numMCRuns)
% RMSEWeighted=sqrt(RMSEWeighted/numMCRuns)
% NEES=NEES/(numMCRuns*2*numDim)
%One will observe that the unweighted RMSE is much higher than the weighted
%RMSE.
%
%EXAMPLE 3:
%In this example, Gaussian noise is added to the measurements. The
%covariance matrix has nonzero off-diagonal elements and halfway through,
%the sign of the off-diagonal elements changes. The optimal weighting of
%the samples is by the inverse of the covariance matrix for that sample.
%Here, we compare the unweighted estimate, with an estimate weighted using
%only the diagonal elements of the optimal weighting matrices. Also, we
%look at the NEES to see that the fully weighted estimate is consistent
%with its covariance matrix.
% numMCRuns=5000;
% ab=[[1;2;-3;4],[-8;0;12;6]];
% numDim=4;
% R=[10,  9, -9,  1;
%     9, 10, -9,  1;
%    -9, -9, 10, -1;
%     1,  1, -1, 10];
% SR=chol(R,'lower');
% RInv=inv(R);
% RAlt=[10,  9,  9, -1;
%       9,  10,  9, -1;
%       9,   9, 10,  1;
%      -1,  -1,  1, 10];
% SRAlt=chol(RAlt,'lower');
% RInvAlt=inv(RAlt);
% N=10;
% 
% RMSEUnweighted=0;
% RMSEWeighted=0;
% RMSEMatrixWeighted=0;
% NEES=0;
% for curRun=1:numMCRuns
%     t=10*randn(N,1);
%     zTrue=polyValVec(ab,t);
%     z=zeros(numDim,N);
%     W=zeros(numDim,numDim,N);
%     for k=1:floor(N/2)
%         %Add noise to the measurements.
%         z(:,k)=zTrue(:,k)+SR*randn(numDim,1);
%         W(:,:,k)=RInv;
%     end
%     for k=(floor(N/2)+1):N
%         z(:,k)=zTrue(:,k)+SRAlt*randn(numDim,1);
%         W(:,:,k)=RInvAlt;
%     end
%     WD=zeros(numDim,numDim,N);
%     for k=1:N
%         WD(:,:,k)=diag(diag(W(:,:,k)));
%     end
% 
%     %Fit the lines.
%     ab0=linearRegression(t,z);
%     ab1=linearRegression(t,z,WD);
%     [ab2,PInv2]=linearRegression(t,z,W);
%     %Evaluate at the fitting points.
%     z0=polyValVec(ab0,t);
%     z1=polyValVec(ab1,t);
%     z2=polyValVec(ab2,t);
%     %Compute RMSE
%     RMSEUnweighted=RMSEUnweighted+calcMSE(z0,zTrue);
%     RMSEWeighted=RMSEWeighted+calcMSE(z1,zTrue);
%     RMSEMatrixWeighted=RMSEMatrixWeighted+calcMSE(z2,zTrue);
%     %Compute NEES
%     diff=ab2(:)-ab(:);
%     NEES=NEES+diff'*PInv2*diff;
% end
% RMSEUnweighted=sqrt(RMSEUnweighted/numMCRuns)
% RMSEWeighted=sqrt(RMSEWeighted/numMCRuns)
% RMSEMatrixWeighted=sqrt(RMSEMatrixWeighted/numMCRuns)
% NEES=NEES/(numMCRuns*2*numDim)
%In this instance, one will typically see that the matrix weighted solution
%does the best, the unweighted solution does the worse, and the solution
%that only uses the diagonal of the weights has the worst performance.
%Also, the NEES will be in the vicinity of 1 indicating consistency.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation: Theory, Algorithms and
%    Software. New York: John Wiley and Sons, 2001.
%
%October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=length(t);
if(nargin<3||isempty(W))
    %In the unweighted case, it is fairly simple. Each dimension can be
    %handled independently.

    C=pinv([t(:),ones(numMeas,1)]);
    ab=z*C';
    PInv=[];
elseif(isvector(W))
    %One weight per measurement. The weight is the same for all of the
    %dimensions.
    numDim=size(z,1);
    
    H=[t(:),ones(numMeas,1)];
    R=diag(W);
    
    ab=zeros(numDim,2);
    PInv=zeros(2*numDim,2*numDim);
    for curDim=1:numDim
        PInvCur=H'*R*H;
        abCur=PInvCur\H'*R*(z(curDim,:).');
        PInv(curDim,curDim)=PInvCur(1,1);
        altDim=numDim+curDim;
        PInv(altDim,altDim)=PInvCur(2,2);
        PInv(curDim,altDim)=PInvCur(1,2);
        PInv(altDim,curDim)=PInv(curDim,altDim);
        ab(curDim,1)=abCur(1);
        ab(curDim,2)=abCur(2);
    end
else
    %The W is a bunch of matrices.
    numDim=size(z,1);
    
    Hk=zeros(numDim*numMeas,2*numDim);
    span=1:numDim;
    for k=1:numMeas 
        Hk(span,:)=[t(k)*eye(numDim,numDim),eye(numDim,numDim)];
        span=span+numDim;
    end
    Rk=blkDiagRep(W);
    
    PInv=Hk'*Rk*Hk;
    ab=PInv\Hk'*Rk*z(:);
    ab=reshape(ab,[numDim,2]);
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
