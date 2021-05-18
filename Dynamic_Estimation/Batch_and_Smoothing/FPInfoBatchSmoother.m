function [yEst,PInvEst,xEst,PEst]=FPInfoBatchSmoother(yPred,PInvPred,z,u,H,F,R,Q,kD)
%%FPINFOBATCHSMOOTHER Run the Fraser-Potter smoother for linear dynamic and
%                measurement models using information filters on a batch
%                of data given initial estimates, which could be
%                uninformative. The smoothed result at one time step or
%                along the entire batch are available.
%
%INPUTS: yPred The predicted information state at the time of the initial
%              measurement in z. If no prior information is available, then
%              just pass an empty matrix.
%     PInvPred The inverse covariance matrix associated with the predicted
%              information state at the time of the initial measurement in
%              z. If not prior information is  avaiable, then just pass an
%              empty matrix.
%            z The zDimXN matrix of measurements for the whole batch.
%            u The xDimX(N-1) matrix of control inputs for the whole
%              batch. If there are no control inputs, then set u=[];
%            H The zDimXxDimXN hypermatrix of measurement matrices
%              such that H(:,:,k)*x+w is the measurement at time k, where
%              x is the state and w is zero-mean Gaussian noise with
%              covariance matrix R(:,:,k). Alternatively, if all of the
%              measurement matrices are the same, one can just pass a
%              single zDim X xDim matrix.
%            F The xDimXxDimX(N-1) hypermatrix of state transition
%              matrices. The state at discrete-time k+1 is modeled as
%              F(:,:,k) times the state at time k plus zero-mean
%              Gaussian process noise with covariance matrix Q(:,:,k).
%              Alternatively, if all of the state transition matrices are
%              the same, one can just pass a single xDimXxDim matrix.
%            R The zDimXzDimXN hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a
%              single zDimXzDim matrix.
%            Q The xDimXxDimX(N-1) hypermatrix of process noise
%              covariance matrices. Alternatively, if all of the process
%              noise covariance matrices are the same, one can just pass a
%              single xDimXxDim matrix. This can be singular.
%           kD The discrete time-step at which the smoothed state estimate
%              is desired, where z(:,1) is at discrete time-step 1 (not 0).
%              If kD is omitted ot an empty matrix is passed, then results
%              along the entire batch are obtained.
%
%OUTPUTS: yEst The xDimXN smoothed information state estimates at all steps
%              if kD is not provided or the xDimX1 smoothed information
%              state estimate at step kD if kD is provided.
%      PInvEst The inverse covariance matrices associated with the smoothed
%              information state estimates. This is xDimXxDimXN for the
%              whole batch if kD is not provided and is xDimXxDim if kD is
%              provided.
%         xEst The state estimates corresponding to yEst.
%         PEst The covariance estimates corresponding to PInvEst.
%
%The Kalman smoothing technique utilizing forward and backwards filters is
%described in [1]. The use with information filters is related.
%
%The algorithm works by running an information filter forward just before
%time kD and an information filter backwards to time kD. The results are
%then merged. The measurement update step remains the same in both the
%forward and the reverse algorithms.
%
%REFERENCES:
%[1] D. C. Fraser and J. E. Potter, "The optimal linear smoother as a 
%    combination of two optimum linear filters," IEEE Transactions on
%    Automatic Control, pp. 387-390, Aug. 1969.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(H,2);
N=size(z,2);

if(nargin<9)
   kD=[]; 
end

if(isempty(yPred))
    yPred=zeros(xDim,1); 
end
if(isempty(PInvPred))
    PInvPred=zeros(xDim,xDim);
end
if(isempty(u))
    u=zeros(xDim,N-1); 
end
if(size(H,3)==1)
    H=repmat(H,[1,1,N]);
end
if(size(F,3)==1)
    F=repmat(F,[1,1,N-1]);
end
if(size(R,3)==1)
    R=repmat(R,[1,1,N]); 
end
if(size(Q,3)==1)
    Q=repmat(Q,[1,1,N-1]);
end

%Run the information filter forward until we have the predictions of step
%kD|kD-1 (going forwards).
yFwdPred=zeros(xDim,N);
PInvFwdPred=zeros(xDim,xDim,N);

%The first step uses the priors
yFwdPred(:,1)=yPred;
PInvFwdPred(:,:,1)=PInvPred;

[yFwd,PInvFwd]=infoFilterUpdate(yPred,PInvPred,z(:,1),R(:,:,1),H(:,:,1));
[yFwdPred(:,2),PInvFwdPred(:,:,2)]=infoFilterDiscPred(yFwd,PInvFwd,F(:,:,1),Q(:,:,1),u(:,1));

%The other steps
for curStep=2:(N-1)
    [yFwd,PInvFwd]=infoFilterUpdate(yFwdPred(:,curStep),PInvFwdPred(:,:,curStep),z(:,curStep),R(:,:,curStep),H(:,:,curStep));
    [yFwdPred(:,curStep+1),PInvFwdPred(:,:,curStep+1)]=infoFilterDiscPred(yFwd,PInvFwd,F(:,:,curStep),Q(:,:,curStep),u(:,curStep));
end

%Run the information filter backwards until we have the updated information
%state of step kD|kD (going backwards).
yRev=zeros(xDim,N);
PInvRev=zeros(xDim,xDim,N);

[yRev(:,end),PInvRev(:,:,end)]=infoFilterUpdate(zeros(xDim,1),zeros(xDim,xDim),z(:,N),R(:,:,N),H(:,:,N));
for curStep=(N-1):-1:1
    [yRevPred,PInvRevPred]=infoFilterDiscPredRev(yRev(:,curStep+1),PInvRev(:,:,curStep+1),F(:,:,curStep),Q(:,:,curStep),u(:,curStep));
    [yRev(:,curStep),PInvRev(:,:,curStep)]=infoFilterUpdate(yRevPred,PInvRevPred,z(:,curStep),R(:,:,curStep),H(:,:,curStep));
end

if(~isempty(kD))
    %If the estimate at only one time is desired.
    PInvEst=PInvFwdPred(:,:,kD)+PInvRev(:,:,kD);
    yEst=yFwdPred(:,kD)+yRev(:,kD);
    
    if(nargout>2)
        xEst=PInvEst\yEst;
        PEst=inv(PInvEst);
    end
else
    PInvEst=PInvFwdPred+PInvRev;
    yEst=yFwdPred+yRev;
    
    if(nargout>2)
        xEst=zeros(xDim,kD);
        PEst=zeros(xDim,xDim,kD);
        
        for kD=1:N
            PEst(:,:,kD)=inv(PInvEst(:,:,kD));
            xEst(:,kD)=PInvEst(:,:,kD)\yEst(:,kD);
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
