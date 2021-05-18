function [xSmooth,PSmooth]=KalmanBatchSmoother(xInit,PInit,z,u,H,F,R,Q,kD,useFP)
%%KALMANBATCHSMOOTHER Run the standard forward-backward Kalman smoother for
%                linear dynamic and measurement models on a batch of
%                measurements. The smoothed result at one time step or
%                along the entire batch are available. The initial
%                predicted states cannot be uninformative. To run the
%                smoother without prior predicted values, use the function
%                FPInfoBatchSmoother.
%
%INPUTS: xInit The predicted state at the time of the initial measurement
%              in z.
%        PInit The covariance matrix associated with the predicted state
%              at the time of the initial measurement in z.
%            z The zDimXN matrix of measurements for the whole batch.
%            u The xDimX(N-1) matrix of control inputs for the whole batch.
%              If there are no control inputs, then set u=[];
%            H The zDimXxDimXN hypermatrix of measurement matrices such
%              that H(:,:,k)*x+w is the measurement at time k, where x is
%              the state and w is zero-mean Gaussian noise with covariance
%              matrix R (:,:,k). Alternatively, if all of the measurement
%              matrices are the same, one can just pass a single
%              zDimXxDim matrix.
%            F The xDimXxDimX(N-1) hypermatrix of state transition
%              matrices. The state at discrete-time k+1 is modeled as
%              F(:,:,k) times the state at time k plus zero-mean
%              Gaussian process noise with covariance matrix Q(:,:,k).
%              Alternatively, if all of the state transition matrices are
%              the same, one can just pass a single xDimXxDim matrix.
%            R The zDimXzDimXN hypermatrix of measurement covariance
%              matrices. Alternatively, if all of the measurement
%              covariance matrices are the same, one can just pass a single
%              zDimXzDim matrix.
%            Q The xDimXxDimX(N-1) hypermatrix of process noise
%              covariance matrices. Alternatively, if all of the process
%              noise covariance matrices are the same, one can just pass a
%              single xDimXxDim matrix.
%           kD The discrete time-step at which the smoothed state estimate
%              is desired, where z(:,1) is at discrete time-step 1 (not 0).
%              If kD is omitted or an empty matrix is passed, then results
%              along the entire batch are obtained.
%        useFP Optional boolean to specify whether to use Fraser-Potter
%              information smoother. If not, standard linear Kalman
%              smoother is used. Default value is true.
%
%OUTPUTS: xEst The xDimXN smoothed state estimates at all steps if kD is
%              not provided or the xDimX1 smoothed information state
%              estimate at step kD if kD is provided.
%         PEst The covariance matrices associated with the smoothed state
%              estimates. This is xDimXxDimXN for the whole batch if kD is
%              not provided and is xDimXxDim if kD is provided.
%
%This function defaults to a wrapper for the function FPInfoSmoother, which
%uses information state and an inverse covariance matrix. This function
%also gives the option to use a non-information filter based Kalman
%smoother as described in Chapter 8.6 of [1].
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%March 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(kD))
    kD=[];
end
if(nargin<10||isempty(useFP))
    useFP=true;
end

if(useFP)
    yInit=PInit\xInit;
    PInvInit=inv(PInit);
    
    [~,~,xSmooth,PSmooth]=FPInfoBatchSmoother(yInit,PInvInit,z,u,H,F,R,Q,kD);
else
    xDim=size(H,2);
    N=size(z,2);
    
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
        [xUpd(:,curStep),PUpd(:,:,curStep)]=KalmanUpdate(xPred(:,curStep),PPred(:,:,curStep),z(:,curStep),R(:,:,curStep),H(:,:,curStep));
        [xPred(:,curStep+1),PPred(:,:,curStep+1)]=discKalPred(xUpd(:,curStep),PUpd(:,:,curStep),F(:,:,curStep),Q(:,:,curStep),u(:,curStep));
    end
    
    %Run the backwards Kalman smoother
    xSmooth=zeros(xDim,N);
    PSmooth=zeros(xDim,xDim,N);
    
    [xSmooth(:,end),PSmooth(:,:,end)]=KalmanUpdate(xPred(:,N),PPred(:,:,N),z(:,N),R(:,:,N),H(:,:,N));
    for curStep=(N-1):-1:1
        C=PUpd(:,:,curStep)*F(:,:,curStep)'/PPred(:,:,curStep+1);
        xSmooth(:,curStep)=xUpd(:,curStep)+C*(xSmooth(:,curStep+1)-xPred(:,curStep+1));
        PSmooth(:,:,curStep)=PUpd(:,:,curStep)+C*(PSmooth(:,:,curStep+1)-PPred(:,:,curStep+1))*C';
    end
    
    if(~isempty(kD))
        xSmooth=xSmooth(:,kD);
        PSmooth=PSmooth(:,:,kD);
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
