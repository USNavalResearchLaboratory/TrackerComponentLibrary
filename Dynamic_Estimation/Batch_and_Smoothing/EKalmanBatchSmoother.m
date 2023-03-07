function [xSmooth,PSmooth,xUpd,PUpd]=EKalmanBatchSmoother(xInit,PInit,z,h,HJacob,f,FJacob,R,Q,kD,numIter)
%%EKALMANBATCHSMOOTHER Run the extended forward-backward Kalman smoother
%                for nonlinear dynamic and measurement models on a batch of
%                measurements. The smoothed result at one time step or
%                along the entire batch is available. The initial predicted
%                states cannot be uninformative.
%
%INPUTS: xInit The predicted state at the time of the initial measurement
%              in z.
%        PInit The covariance matrix associated with the predicted state at
%              the time of the initial measurement in z.
%            z The zDim X N matrix of measurements for the whole batch.
%            h A NX1 cell array of function handles for the measurement
%              function that transform the state into the measurement
%              domain at each step. If the same measurement function is
%              used for all steps in the batch, then h can just be the
%              single function handle used.
%       HJacob A NX1 cell array of function handles for the measurement
%              Jacobian matrix that each takes the target state as a
%              parameter. If a single measurement Jacobian matrix is used
%              for all steps of the batch, then HJacob can just be the
%              single function handle used. If an empty matrix is passed or
%              the parameter is omitted, then HJacob will be found using
%              numerical differentiation via the numDiff function with
%              default parameters.
%            f A NX1 cell array of function handles for the state
%              transition function that transform the state into the
%              measurement domain at each step. If the same measurement
%              function is used for all steps in the batch, then h can just
%              be the single function handle used.
%       FJacob A NX1 cell array of function handles for the state
%              transition Jacobian matrix that each takes the target state
%              as a parameter. If a single measurement Jacobian matrix is
%              used for all steps of the batch, then HJacob can just be the
%              single function handle used. If an empty matrix is passed or
%              the parameter is omitted, then HJacob will be found using
%              numerical differentiation via the numDiff function with
%              default parameters.
%            R The zDim X zDim X N hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a single
%              zDim X zDim matrix.
%            Q The xDim X xDim X (N-1) hypermatrix of process noise
%              covariance matrices. Alternatively, if all of the process
%              noise covariance matrices are the same, one can just pass a
%              single xDim X xDim matrix.
%           kD The discrete time-step at which the smoothed state estimate
%              is desired, where z(:,1) is at discrete time-step 1 (not 0).
%              If kD is omitted or an empty matrix is passed, then
%              results along the entire batch are obtained.
%      numIter The number of iterations to perform if an iterated EKF is
%              desired. The default is zero. That is, just use the
%              standard update without any additional iterations.
%
%OUTPUTS: xEst The xDimXN smoothed state estimates at all steps if kD is
%              not provided or the xDimX1 smoothed state estimate at step
%              kD if kD is provided.
%         PEst The covariance matrices associated with the smoothed
%              state estimates. This is xDimXxDimXN for the whole batch
%              if kD is not provided and is xDimXxDim if kD is provided.
%         xUpd The xDimXN state estimates of the forward filter (not
%              smoothed) at all times.
%         PUpd The xDimXxDimXN covariance matrices of the forward filter
%              (not smoothed) at all times.
%
%This function implements an extended Kalman smoother modified from Chapter
%8.6 of [1].
%
%The smoothing iteration algorithm is taken from Section III of [2]. This
%provides an algorithm to iterate the smoothed state, but not the smoothed
%covariance matrix.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] L. A. Johnston and V. Krishnamurthy, "Derivation of a sawtooth
%    iterated extended Kalman smoother via the AECM algorithm," IEEE
%    Transactions on Signal Processing, vol. 49, no. 9, pp. 1899-1909,
%    2001.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if nargin<10
    kD=[];
end
if nargin<11
    numIter=0;
end

xDim=size(xInit,1);
zDim=size(z,1);
N=size(z,2);

if(isempty(FJacob))
    FJacob=@(x)numDiff(x,f,xDim);
end
if(isempty(HJacob))
    HJacob=@(x)numDiff(x,h,zDim);
end

if(isa(HJacob,'function_handle'))
    HJacob=repmat({HJacob},N,1);
end
if(isa(h,'function_handle'))
    h=repmat({h},N,1);
end
if(isa(FJacob,'function_handle'))
    FJacob=repmat({FJacob},N,1);
end
if(isa(f,'function_handle'))
    f=repmat({f},N,1);
end
if(size(R,3)==1)
    R=repmat(R,1,1,N);
end
if(size(Q,3)==1)
    Q=repmat(Q,1,1,N-1);
end

%Run the Kalman filter forward
xPred=zeros(xDim,N);
PPred=zeros(xDim,xDim,N);
xUpd=zeros(xDim,N);
PUpd=zeros(xDim,xDim,N);

%The first step uses the priors
xPred(:,1)=xInit;
PPred(:,:,1)=PInit;

%The rest of the steps
for curStep=1:(N-1)
    [xUpd(:,curStep),PUpd(:,:,curStep)]=EKFUpdate(xPred(:,curStep),PPred(:,:,curStep),z(:,curStep),R(:,:,curStep),h{curStep},HJacob{curStep},numIter);
    [xPred(:,curStep+1),PPred(:,:,curStep+1)]=DiscEKFPred(xUpd(:,curStep),PUpd(:,:,curStep),f{curStep},FJacob{curStep},Q(:,:,curStep));
end
[xUpd(:,end),PUpd(:,:,end)]=EKFUpdate(xPred(:,N),PPred(:,:,N),z(:,N),R(:,:,N),h{N},HJacob{N},numIter);

%Run the backwards Kalman smoother
xSmooth=zeros(xDim,N);
PSmooth=zeros(xDim,xDim,N);

xSmooth(:,N)=xUpd(:,N);
PSmooth(:,:,N)=PUpd(:,:,N);
for curStep=(N-1):-1:1
    F=FJacob{curStep}(xUpd(:,curStep));
    C=PUpd(:,:,curStep)*F'/PPred(:,:,curStep+1);
    xSmooth(:,curStep)=xUpd(:,curStep)+C*(xSmooth(:,curStep+1)-xPred(:,curStep+1));
    PSmooth(:,:,curStep)=PUpd(:,:,curStep)+C*(PSmooth(:,:,curStep+1)-PPred(:,:,curStep+1))*C';
    
    for curIter=1:numIter
        F=FJacob{curStep}(xSmooth(:,curStep));
        H=HJacob{curStep}(xSmooth(:,curStep));
        B=pinv(pinv(PUpd(:,:,curStep))+H'*inv(R(:,:,curStep))*H+F'*inv(Q(:,:,curStep))*F);
        xSmooth(:,curStep)=xUpd(:,curStep)+B*(...
            H'*inv(R(:,:,curStep))*(z(:,curStep)-h{curStep}(xPred(:,curStep))-H*(xPred(:,curStep)-xSmooth(:,curStep)))...
            +F'*inv(Q(:,:,curStep))*(xSmooth(:,curStep+1)-xPred(:,curStep+1)));
    end
end

if(~isempty(kD))
    xSmooth=xSmooth(:,kD);
    PSmooth=PSmooth(:,:,kD);
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
