function [xEst,PEst]=KalmanFIRSmoother(z,u,H,F,R,Q,kD)
%%KALMANFIRSMOOTHER Run the finite impulse response Kalman smoother on a
%                   batch of data.
%
%INPUTS: z The zDim X N matrix of measurements for the whole batch.
%        u The xDim X(N-1) matrix of control inputs for the whole batch. If
%          there are no control inputs, then set u=[];
%        H The zDim X xDim X N hypermatrix of measurement matrices such
%          that H(:,:,k)*x+w is the measurement at time k, where x is the
%          state and w is zero-mean Gaussian noise with covariance matrix
%          R(:,:,k). Alternatively, if all of the measurement matrices are
%          the same, one can just pass a single zDim X xDim matrix.
%        F The xDim X xDim X (N-1) hypermatrix of state transition
%          matrices. The state at discrete-time k is modeled as F(:,:,k)
%          times the state at time k plus zero-mean Gaussian process noise
%          with covariance matrix Q(:,:,k). Alternatively, if all of the
%          state transition matrices are the same, one can just pass a
%          single xDim X xDim matrix.
%        R The zDim X zDim X N hypermatrix of measurement covariance
%          matrices. Alternatively, if all of the measurement covariance
%          matrices are the same, one can just pass a single zDim X zDim
%          matrix.
%        Q The xDimXxDimX(N-1) hypermatrix of process noise covariance
%          matrices. Alternatively, if all of the process noise covariance
%          matrices are the same, one can just pass a single xDimXxDim
%          matrix.
%       kD The discrete time-step at which the smoothed state estimate is
%          desired, where z(:,1) is at discrete time-step 1 (not 0).
%
%OUTPUTS: xEst The smoothed state estimate at step kD.
%         PEst The covariance matrix associated with the smoothed state
%              estimate at step kD.
%
%The assumed forward-time dynamic equations are
%x(:,k)=F(:,:,k-1)*x(:,k-1)+u(:,k-1)+noise
%z(k)=H(:,:,k)*x(:,k)+noise
%where x is the target state, u is a control input, and z is the
%measurement.
%
%The algorithm is that given in [1]. Note that if H, F, R, Q and kD are
%fixed, the matrices produced by the function KalmanFIRSmootherCoeffs,
%which is called by this function, can be calculated once and used
%repeatedly for different sets of measurements z and control inputs u. The
%matrix G in the paper is omitted, since any control input can be pre-
%multipled by the matrix.
%
%REFERENCES:
%[1] D. F. Crouse, P. Willett, and Y. Bar-Shalom, "A low-complexity
%    sliding-window Kalman FIR smoother for discrete-time models," IEEE
%    Signal Processing Letters, vol. 17, no. 2, pp. 177-180, Feb. 2009.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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

[A, B,PEst]=KalmanFIRSmootherCoeffs(H,F,R,Q,kD);

xEst=zeros(xDim,1);
for idx=1:(N-1)
    xEst=xEst+A(:,:,idx)*z(:,idx)+B(:,:,idx)*u(:,idx);
end

%No control input is used at the final step.
xEst=xEst+A(:,:,end)*z(:,end);
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
