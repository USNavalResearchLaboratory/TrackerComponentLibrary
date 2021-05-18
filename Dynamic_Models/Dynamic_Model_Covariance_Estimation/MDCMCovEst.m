function [R,Q]=MDCMCovEst(z,H,F,u,maxLag)
%%MDCMCOVEST Use the measurement difference correlation method (MDCM) to
%           jointly estimate the measurement covariance matrix and the
%           process noise covariance of a linear time invariant (LTI)
%           system. For observability, this requires that
%           zDim*xDim>=xDim^2. That is, that the measurement dimensionality
%           is >= the state dimensionality. The assumed measurement model
%           if x_k is the target state at time k is:
%           z_k=H*x_k+w_k
%           where H is the measurement matrix and w_k is zero-mean Gaussian
%           measurement noise with covariance matrix R. The dynamic model
%           is:
%           x_{k+1}*F*x_k+v_k+u_k
%           where v_k is zero-mean Gaussian process noise with covariance
%           matrix Q and u_k is an optional control input.
%
%INPUTS: z The zDimXN set of N measurements over time. N should, at a
%          minimum, be large enough for the covariance matrices being
%          estimated to be observable. This goes beyond the observability
%          of the target state. The exact number of measurements necessary
%          depends on F and H.
%        H The zDimXxDim measurement matrix.
%        F The xDimXxDim state transition matrix. This is assumed the same
%          across all discrete time steps.
%        u An optional xDimX(N-1) set of control inputs for state
%          prediction in the filter. An empty matrix can be passed if no
%          control inputs are used.
%   maxLag The maximum delay used for the autocorrelation of the
%          innovations. This must be <=N-1. The default if omitted or an
%          empty matrix is passed is 2. If this parameter is set to be too
%          small or too large, then numerical or observability issues will
%          arise. It should generally be set to be at least the order of
%          the system.
%
%OUTPUTS: R The zDimXzDim estimate of R.If a non-finite term is encountered
%           (usually due to finite-precision or observability issues), then
%           an empty matrix will be returned.
%         Q The xDimXxDim estimate of Q. If a non-finite term is encountered
%           (usually due to finite-precision or observability issues), then
%           an empty matrix will be returned.
%
%The algorithm of [1] is used. Note that in [1], there is also a variant
%for linear time-varying systems. The algorthm is also discussed in less
%detail in [2].
%
%EXAMPLE:
%Here we demonstrate the estimation of both covariance matrices for a
%linear system.
% numRuns=10;%Number of Monte Carlo runs to average.
% numSteps=1000;
% R=[3,   -0.1,   0,  0;
%   -0.1,    2,   0,  0;
%     0,     0,   1,  0.2;
%     0,     0,   0.2,0.08];%True Cartesian measurement covariance matrix.
% SR=chol(R,'lower');
% zDim=size(R,1);
% H=eye(zDim,4);
% 
% xDim=size(H,2);
% T=1/2;
% F=FPolyKal(T,xDim,1);
% q0=1e-3;
% Q=QPolyKal(T,xDim,1,q0);
% SQ=chol(Q,'lower');
% x0=[10e3;-2e3;100;200];%True initial state as [2D position; 2D velocity].
% 
% maxLag=6;
% frobErrR=0;
% frobErrQ=0;
% for curRun=1:numRuns
%     %Generate the sequence of measurements.
%     z=zeros(zDim,numSteps);
%     x=x0;
%     for curStep=1:numSteps
%         z(:,curStep)=H*x+SR*randn(zDim,1);
%         x=F*x+SQ*randn(xDim,1);
%     end
% 
%     [REst,QEst]=MDCMCovEst(z,H,F,[],maxLag);
%     
%     frobErrR=frobErrR+norm(REst-R,'fro');
%     frobErrQ=frobErrQ+norm(QEst-Q,'fro');
% end
% frobErrR=frobErrR/numRuns;
% frobErrQ=frobErrQ/numRuns;
% frobErrRRelative=frobErrR/norm(R,'fro')
% frobErrQRelative=frobErrQ/norm(Q,'fro')
%The relative Frobenius norm error. We find that the error for Q is rather
%bad. One can see that changing q0 above affects the error, with better
%results for larger q0 values. This indicates that in this scenario, the
%algorithm is bad when dealing with very small process noises.
%
%REFERENCES:
%[1] J. Duník, O. Straka, and O. Kost, "Measurement difference
%    autocovariance method for noise covariance matrices estimation," in
%    Proceedings of the IEEE 55th Conference on Decision and Control, Las
%    Vegas, NV, 12-14 Dec. 2016, pp. 365-370.
%[2] J. Duník, O. Straka, O. Kost, and J. Havlík, "Noise covariance
%    matrices in state-space models: A survey and comparison of estimation
%    methods-part I," International Journal of Adaptive Control and Signal
%    Processing, vol. 31, no. 11, pp. 1505-1543, Nov. 2017.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(H,2);
zDim=size(z,1);
N=size(z,2);

if(nargin<5||isempty(maxLag))
    maxLag=2;
end

if(nargin<4||isempty(u))
    u=zeros(xDim,N);
end

if(zDim<xDim)
    error('The input system does not have observability, because zDim<xDim.')
end

%Precompute powers of F.
HFPows=zeros(zDim,xDim,maxLag+1);
FPow=eye(xDim,xDim);
HFPows(:,:,1)=H*FPow;
for k=1:(maxLag+1)
    FPow=F*FPow;
    HFPows(:,:,k+1)=H*FPow;
end

%precompute the pseudoinverse of H.
HInv=pinv(H);

%Estimate the prediction innovations using Equations 3.31 and 3.32 in [2].
%This is equivalent to Equation 32 in [1] with an added term for a control
%input.
innovPred=zeros(zDim,N);
CoeffMat=HFPows(:,:,maxLag+1)*HInv;%precompute
for k=(maxLag+1):N
    sumVal=0;
    for i=0:(maxLag-1)
        sumVal=sumVal+HFPows(:,:,i+1)*u(:,k-i-1);
    end

    innovPred(:,k)=z(:,k)-CoeffMat*z(:,k-maxLag)-sumVal;
end

%Compute the b vector of Equation 29 of [1].
offset=2*maxLag+1;
b=zeros(zDim^2*(maxLag+1),N-offset+1);
for i=(2*maxLag+1):N
    startIdx=1;
    for j=0:maxLag
        temp=innovPred(:,i)*innovPred(:,i-j)';
        b(startIdx:(startIdx+zDim^2-1),i-offset+1)=temp(:);
        startIdx=startIdx+zDim^2;
    end
end
b=mean(b,2);%Average of the autocorrelation matrices.

%Compute the A matrix of Equation 33 of [1].
AMat=zeros(zDim^2*(maxLag+1),xDim^2+zDim^2);
%The first batch of rows is a special case.
A1=zeros(zDim^2,xDim^2);
for k=0:(maxLag-1)
    A1=A1+kron(HFPows(:,:,k+1),HFPows(:,:,k+1));
end
A2=kron(HFPows(:,:,maxLag+1)*HInv,HFPows(:,:,maxLag+1)*HInv)+eye(zDim^2,zDim^2);
AMat(1:zDim^2,:)=[A1,A2];

curStart=1+zDim^2;
for rowBatch=2:maxLag
    A1=zeros(zDim^2,xDim^2);
    for k=0:(maxLag-rowBatch)
        A1=A1+kron(HFPows(:,:,k+1),HFPows(:,:,k+(rowBatch-1)+1));
    end
    AMat(curStart:(curStart+zDim^2-1),:)=[A1,zeros(zDim^2,zDim^2)];
    curStart=curStart+zDim^2;
end
%The last batch of rows is also a special case
A2=kron(eye(zDim,zDim),-HFPows(:,:,maxLag+1)*HInv);
AMat(curStart:(curStart+zDim^2-1),:)=[zeros(zDim^2,xDim^2),A2];

%Least squares solution for Q and R from Equation 21.
QR=linsolve(AMat,b);
Q=reshape(QR(1:xDim^2),xDim,xDim);
Q=(Q+Q')/2;%Force to be symmetric.

if(any(~isfinite(Q(:))))
    %Finite precision (or observability) isues arose.
    Q=[];
    R=[];
    return;
end

%Force it to be a valid covariance value (no negative eigenvalues).
[V,D]=eig(Q);
Q=V*abs(D)*V';

R=reshape(QR((xDim^2+1):(xDim^2+zDim^2)),zDim,zDim);
R=(R+R')/2;%Force to be symmetric.

if(any(~isfinite(R(:))))
    %Finite precision (or observability) isues arose.
    Q=[];
    R=[];
    return;
end

%Force it to be a valid covariance value (no negative eigenvalues).
[V,D]=eig(R);
R=V*abs(D)*V';
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
