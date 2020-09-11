function [ySqrtEst, PInvSqrtEst]=sqrtInfoBatchSmoother(ySqrtPred,PInvSqrtPred,z,u,H,F,SR,SQ,Gamma,kD)
%%SQRTINFOBATCHSMOOTHER Run the square root information smoother for linear
%                  dynamic and measurement models on a batch of
%                  measurements. The smoothed result at one time step or
%                  along the entire batch are available. The initial
%                  predicted states can not be uninformative.
%
%INPUTS: ySqrtPred The xDimX1 predicted square root information state at
%                 the time of the initial measurement in z. The predicted
%                 information state is always PInvSqrtPred times the
%                 predicted target state estimate.
%    PInvSqrtPred The inverse square root information matrix associated
%                 with the predicted information state at the time of the
%                 initial measurement in z. If P is the covariance matrix
%                 of a Gaussian state x, then P=PSqrt*PSqrt' and
%                 PInvSqrtPred=inv(PSqrt). This can be either upper
%                 triangular or lower triangular.
%               z The zDim X N matrix of measurements for the whole batch.
%               u The xDim X(N-1) matrix of control inputs for the whole
%                 batch. If there are no control inputs, then set u=[];
%               H The zDim X xDim X N hypermatrix of measurement matrices
%                 such that H(:,:,k)*x+w is the measurement at time k, 
%                 where x is the state and w is zero-mean Gaussian noise 
%                 with covariance matrix R (:,:,k). Alternatively, if all
%                 of the measurement matrices are the same, one can just
%                 pass a single zDim X xDim matrix.
%               F The xDim X xDim X (N-1) hypermatrix of state transition
%                 matrices. The state at discrete-time k+1 is modeled as
%                 F(:,:,k) times the state at time k plus zero-mean
%                 Gaussian process noise with covariance matrix Q(:,:,k).
%                 Alternatively, if all of the state transition matrices
%                 are the same, one can just pass a single xDim X xDim
%                 matrix.
%              SR The zDim X zDim X N hypermatrix of lower-triangular
%                 square-root of the measurement covariance matrices. The
%                 matrices must be invertible. Alternatively, if all of the
%                 measurement covariance matrices are the same, one can
%                 just pass a single zDim X zDim matrix.
%              SQ The xDim X xDim X (N-1) hypermatrix of lower-triangular
%                 square-root of the process noise covariance matrices.
%                 The matrices must be invertible. Alternatively, if all of
%                 the measurement covariance matrices are the same, one can
%                 just pass a single xDim X xDim matrix.
%           Gamma An optional xDim X xDim X (N-1) hypermatrix of matrices
%                 that transform the process noise to the state domain if
%                 the process noise covariance matrix is singular.
%                 Alternatively, if all of the transform matrices are the
%                 same, one can just pass a single xDim X xDim matrix. If
%                 this is omitted entirely, an identity matrix is used
%                 (i.e. there is no Gamma).
%              kD The discrete time-step at which the smoothed state
%                 estimate is desired, where z(:,1) is at discrete
%                 time-step 1 (not 0). If kD is omitted ot an empty matrix
%                 is passed, then results along the entire batch are
%                 obtained.
%
%OUTPUTS: ySqrtEst The xDimXN smoothed square root information state
%                  estimates at all steps if kD is not provided or the
%                  xDimX1 smoothed information state estimate at step kD if
%                  kD is provided.
%      PSqrtInvEst The inverse square root covariance matrices associated
%                  with the smoothed information state estimates. This is
%                  xDimXxDimXN for the whole batch if kD is not provided
%                  and is xDimXxDim if kD is provided.
%
%The algorithm is that of the linear square root information filter and
%smoother that is described in the book [1]. The forward pass steps can be
%found in Chapter IV, and the backwards smoother in Chapter X.
%
%The algorithm works by running a square root information filter forward
%in time and then smoothing backwards. The smoothing step depends on
%information stored during the forward pass. The smoothing step here is
%performed on a state residual rather than the state itself and applies the
%smoothed output to the previous updated measurement. This was based off of
%NASA's Orbit Determination Toolbox.
%
%z*(j+1)=PInvSqrt*(j+1)(x*(j+1)-xPred(j+1)
%A=[Rw(j)+Rwx(j),   Rwx(j)F(j+1,j),         u(j)
%   PInvSqrt*(j+1), PInvSqrt*(j+1)F(j+1,j), z*(j+1)]
%qr(A)=[Rw*(j),   Rwx*(j),      u*(j)
%       0,        PInvSqrt*(j), z*(j)]
%x*(j)=xUpd(j)+PInvSqrt*(j)\z*(j)
%
%where xPred is the forward predicted state, xUpd is the forward updated
%state, Rw and Rwx are outputs from the forward prediction step, u is the
%the control input (often zero), and x* and PInvSqrt* are the smoothed
%states. x*(N)=xUpd(N) and PInvSqrt*(N)=PInvSqrtUpd(N)
%
%REFERENCES:
%[1] G. J. Bierman, "Factorization Methods for Discrete Sequential
%    Estimation. Academic Press, New York, 1977.
%
%March 2015, David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(H,2);
N=size(z,2);

if(nargin<9||isempty(Gamma))
    Gamma=eye(xDim);
end
if(nargin<10||isempty(kD))
    kD=[];
end
if(isempty(ySqrtPred))
    ySqrtPred=zeros(xDim,1); 
end
if(isempty(PInvSqrtPred))
    PInvSqrtPred=zeros(xDim,xDim);
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
if(size(SR,3)==1)
    SR=repmat(SR,[1,1,N]); 
end
if(size(SQ,3)==1)
    SQ=repmat(SQ,[1,1,N-1]);
end
if(size(Gamma,3)==1)
    Gamma=repmat(Gamma,[1,1,N-1]);
end

%Run the SRIF forward until we have the predictions of step kD|kD-1 (going forwards).
ySqrtFwdPred=zeros(xDim,N);
PInvSqrtFwdPred=zeros(xDim,xDim,N);
ySqrtFwdUpd=zeros(xDim,N);
PInvSqrtFwdUpd=zeros(xDim,xDim,N);
Rw=zeros(xDim,xDim,N);
Rwx=zeros(xDim,xDim,N);

%The first step uses the priors
ySqrtFwdPred(:,1)=ySqrtPred;
PInvSqrtFwdPred(:,:,1)=PInvSqrtPred;

%The rest of the steps
for curStep=1:(N-1)
    [ySqrtFwdUpd(:,curStep),PInvSqrtFwdUpd(:,:,curStep)]=sqrtInfoFilterUpdate(ySqrtFwdPred(:,curStep),PInvSqrtFwdPred(:,:,curStep),z(:,curStep),SR(:,:,curStep),H(:,:,curStep));
    [ySqrtFwdPred(:,curStep+1),PInvSqrtFwdPred(:,:,curStep+1),Rw(:,:,curStep+1),Rwx(:,:,curStep+1)]=sqrtInfoFilterDiscPred(ySqrtFwdUpd(:,curStep),PInvSqrtFwdUpd(:,:,curStep),F(:,:,curStep),SQ(:,:,curStep),u(:,curStep),Gamma(:,:,curStep));
end

%Run the SRIS backwards until we have the updated information state of
%step kD|kD (going backwards).
ySqrtEst=zeros(xDim,N);
PInvSqrtEst=zeros(xDim,xDim,N);

%ODTBX method
[ySqrtEst(:,end),PInvSqrtEst(:,:,end)]=sqrtInfoFilterUpdate(ySqrtFwdPred(:,N),PInvSqrtFwdPred(:,:,N),z(:,N),SR(:,:,N),H(:,:,N));
for curStep=(N-1):-1:1
    xPred=PInvSqrtFwdPred(:,:,curStep+1)\ySqrtFwdPred(:,curStep+1); %Forward predicted state_j+1
    xStar=PInvSqrtEst(:,:,curStep+1)\ySqrtEst(:,curStep+1); %Backward smoothed state_j+1
    zStar=PInvSqrtEst(:,:,curStep+1)*(xStar-xPred);
    
    [zSmooth,PInvSqrtEst(:,:,curStep)]=sqrtInfoSmoothStep(zStar,PInvSqrtEst(:,:,curStep+1),F(:,:,curStep),Rw(:,:,curStep),Rwx(:,:,curStep),u(:,curStep),Gamma(:,:,curStep));
    
    xUpd=PInvSqrtFwdUpd(:,:,curStep)\ySqrtFwdUpd(:,curStep); %Forward updated state_j
    ySqrtEst(:,curStep)=PInvSqrtEst(:,:,curStep)*(xUpd+PInvSqrtEst(:,:,curStep)\zSmooth);
end

if(~isempty(kD))
    ySqrtEst=ySqrtEst(:,kD);
    PInvSqrtEst=PInvSqrtEst(:,:,kD);
end
end

function [ySqrtPred, PInvSqrtPred]=sqrtInfoSmoothStep(ySqrtPrev,PInvSqrtPrev,F,Ru,Rux,u,Gamma)
%This is the basic Smoothing step shown in Bierman's book
xDim=size(ySqrtPrev,1);

A=[Ru+Rux*Gamma, Rux*F, u;
    PInvSqrtPrev*Gamma,   PInvSqrtPrev*F,   ySqrtPrev];
[~,T]=qr(A);

ySqrtPred=T((xDim+1):end,end);
PInvSqrtPred=T((xDim+1):end,(end-xDim):(end-1));
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
