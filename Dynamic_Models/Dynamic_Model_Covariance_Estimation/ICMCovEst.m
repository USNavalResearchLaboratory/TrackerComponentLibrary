function [R,Q]=ICMCovEst(z,H,F,u,params2Est,num2Skip,maxLag,xInit,RInit,QInit,numIters)
%%ICMCOVEST Use the indirect correlation method (ICM) to estimate the
%           measurement covariance matrix and/or the process noise
%           covariance matrix of a linear time-invariant (LTI) dynamic
%           system given a sequence of measurements. Only up to xDimXzDim
%           elements of the process noise covariance matrix Q can be
%           uniquely determined, so in many systems, this is not a useful
%           estimator of Q. This algorithm relies on the fact that the
%           autocorrelation of the innovations in a matched Kalman filter
%           Kalman filter should be a zero-mean stochastic process. The
%           assumed measurement model if x_k is the target state at time k
%           is:
%           z_k=H*x_k+w_k
%           where H is the measurement matrix and w_k is zero-mean Gaussian
%           measurement noise with covariance matrix R. The dynamic model
%           is:
%           x_{k+1}*F*x_k+v_k+u_k
%           where v_k is zero-mean Gaussian process noise with covariance
%           matrix Q and u_k is an optional control input. The algorithm of
%           [1] first runs a mismatched constant filter with a constant
%           gain to obtain the innovations used by a subsequent estimation
%           step.
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
% params2Est A parameter indicating whether R and/or Q should be estimated.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Just
%            estimate the measurement covariance matrix R.
%          1 Just estimate the process noise covariance matrix Q.
%          2 Estimate both R and Q.
% num2Skip The results in [1] are for asymptotic results. If one wishes to
%          skip a few of the initial innovation values so that the
%          innovation sequence used is more characteristic of the
%          asymptotic behaviour, this parameter can be set to something
%          other than 0. Again, one must retain enough values such that the
%          parameters to be estimated remain observable. The default if
%          omitted or an empty matrix is passed is 0.
%   maxLag The maximum delay used for the autocorrelation of the
%          innovations. This must be <=N-1. The default if omitted or an
%          empty matrix is passed is 2. If this parameter is set to be too
%          small or too large, then numerical or observability issues will
%          arise. It should generally be set to be at least the order of
%          the system.
% xInit, RInit, QInit The value of the target state to use as the predicted
%          target state at the first time instant as well as nominal values
%          of the R and Q matrices (or actual values if one or the other is
%          known). After the first iteration, xInit is updated. The R, Q
%          parameters begin estimated, as specified by the params2Est
%          parameter, are updated at the end fo each iteration. The
%          defaults if omitted or empty matrices are passed as
%          zeros(xDim,1), eye(zDim,zDim), and eye(xDim,xDim).
% numIters The number of iterations of the algorithm to perform. This value
%          must be >=1. The default if omitted or an empty matrix is passed
%          is 3.
%
%OUTPUTS: R The zDimXzDim estimate of R. If params2Est=1, then this is just
%           RInit. If a non-finite term is encountered (usually due to
%           finite-precision or observability issues), then an empty matrix
%           will be returned.
%         Q The xDimXxDim estimate of Q. If params2Est==0 then this is just
%           QInit. If a non-finite term is encountered (usually due to
%           finite-precision or observability issues), then an empty matrix
%           will be returned.
%
%This function implements the algorithm of [1], which is also discussed in
%[2]. The design of the algorithm requires a stable Kalman filter
%predictor gain as a design parameter as well as an initial estimate.
%Setting too terrible an initial estimate can result in very bad estimates.
%To avoid this and to always get a stable gain, a nominal R and Q are used
%to get the optimal Kalman filter gain for a specific filter. The filter
%gain depends on the asymptotic state prediction covariance matrix and is
%obatained by solving the Riccatti equation. To get rid of the reliance of
%the estimates on the initial estimate x, an option has been made to
%iterate the estimator: After estimating R and Q, a Fisher-Potter
%information filter smoother is used to get an improved initial estimate of
%x and a new gain is computed using the estimates of R and Q. Also, since
%the results of [1] are only asymptotically valid, an option num2Skip is
%given so that the first few innovations in the sequence can be
%disregarded.
%
%EXAMPLE:
%Here we demonstrate the estimation of a covariance matrix. The process
%noise covariance matrix is unobservable and is thus provided, not
%estimated. However, it is worth noting that the estimates of R do not
%depend much on Q, so passing a very bad value of Q will not affect the
%estimate of R much here.
% numRuns=10;%Number of Monte Carlo runs to average.
% numSteps=1000;
% R=[3,   -0.1;
%   -0.1, 2];%True Cartesian measurement covariance matrix.
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
% params2Est=0;%Estimate the measurement covariance matrix.
% num2Skip=0;
% maxLag=4;
% numIters=3;
% 
% frobErr=0;
% for curRun=1:numRuns
%     %Generate the sequence of measurements.
%     z=zeros(zDim,numSteps);
%     x=x0;
%     for curStep=1:numSteps
%         z(:,curStep)=H*x+SR*randn(zDim,1);
%         x=F*x+SQ*randn(xDim,1);
%     end
%     
%     REst=ICMCovEst(z,H,F,[],params2Est,num2Skip,maxLag,[],[],Q,numIters);
%     
%     frobErr=frobErr+norm(REst-R,'fro');
% end
% frobErr=frobErr/numRuns;
% frobErrRelative=frobErr/norm(R,'fro')
%We display the relative Frobenius norm of the estimate of R.
%
%REFERENCES:
%[1] R. K. Mehra, "On the identification of variances and adaptive Kalman
%    filtering," IEEE Transactions on Automatic Control, vol. 15, no. 2,
%    pp. 175-184, Apr. 1972.
%[2] J. Duník, O. Straka, O. Kost, and J. Havlík, "Noise covariance
%    matrices in state-space models: A survey and comparison of estimation
%    methods-part I," International Journal of Adaptive Control and Signal
%    Processing, vol. 31, no. 11, pp. 1505-1543, Nov. 2017.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Given z and u.
xDim=size(F,1);
zDim=size(H,1);
N=size(z,2);%The number of measurements.

if(nargin<4||isempty(u))
    u=zeros(xDim,N-1);
end

if(nargin<5||isempty(params2Est))
    params2Est=0;%Estimate just R.
end

if(nargin<6||isempty(num2Skip))
    num2Skip=0;
end

if(nargin<7||isempty(maxLag))
    maxLag=2;
end

switch(params2Est)
    case 0
        estR=true;
        estQ=false;
    case 1
        estR=false;
        estQ=true;
    case 2
        estR=true;
        estQ=true;
    otherwise
        error('Unnown value for params2Est specified.') 
end

if(estQ&&zDim^2<xDim*zDim)
    warning('Estimation of Q is requested, but the measurement dimensionality is too small to ensure observability.')
end

if(nargin<8||isempty(xInit))
    xInit=zeros(xDim,1);
end

if(nargin<9||isempty(RInit))
    RInit=eye(zDim,zDim);
end

if(nargin<10||isempty(QInit))
    QInit=eye(xDim,xDim);
end

if(nargin<11||isempty(numIters))
    numIters=3;
end

%Compute powers of F outside of the loop so that they do not have to be
%recomputed.
FPow=zeros(xDim,xDim,maxLag+1);
FPow(:,:,1)=eye(xDim,xDim);
for k=1:maxLag
    FPow(:,:,k+1)=F*FPow(:,:,k);
end

R=RInit;
Q=QInit;
x0=xInit;
for curIter=1:numIters
    %Given  a new R and Q estimate, we can try to improve xInit. This wil
    %bias the estimator, but any initial estimate will bias the estimator.
    if(curIter>1)
        [~,~,x0,~]=FPInfoBatchSmoother([],[],z,u,H,F,R,Q,1);
    end

    %Get the asymptotic prediction covariance given the current R and Q
    %estimates.
    PPred=RiccatiPredNoClutter(H,F,R,Q);

    %Get the Kalman filter gain corresponding to the asymptotic prediction
    %covariance matrix.
    [~,PzPred,otherInfo]=KalmanMeasPred([],PPred,H);
    K=calcKalmanGain(R,PzPred,otherInfo);

    innov=getInnovSequence(x0,z,H,F,K,u,num2Skip);

    %Compute all of the autocorrelation matrices from the innovation.
    %Equation 3.6 in [2].
    C=autocorrMats(innov,maxLag);

    %Defined after Equation 3.1b.
    FBar=F-F*K*H;

    KC0=K*C(:,:,1);
    %Solve Equation 3.4 in [2] for P*H'. This requires building a system of
    %Equations. KC0 is xDimXzDim. H*FBar^(j-1)*F*K*C0 is zDimXzDim.
    %P*H' is xDimXzDim.

    A=zeros(zDim*N,xDim);
    B=zeros(zDim*N,zDim);

    FBarPow=eye(xDim,xDim);
    startIdx=1;
    for j=1:maxLag
        coeffMat=H*FBarPow*F;

        A(startIdx:(startIdx+zDim-1),:)=coeffMat;
        B(startIdx:(startIdx+zDim-1),:)=C(:,:,j+1)+coeffMat*KC0;

        FBarPow=FBarPow*FBar;
        startIdx=startIdx+zDim;
    end
    %Get the least-squares solution of P*H' from Equation 3.4.
    PHT=linsolve(A,B);

    if(estR)
        %Solve for R in Equation 3.4
        R=C(:,:,1)-H*PHT;
        %Force symmetry.
        R=(R+R')/2;
        
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

    if(estQ)
        %Equation 29 in [1].
        Omega=F*(-K*PHT'-PHT*K'+K*C(:,:,1)*K')*F';

        %This will hold all of the right-hand-side equations from Equation
        %28.
        RhsMatEq28=zeros(maxLag*zDim^2,1); 
        startIdx=1;
        for k=1:maxLag
            %rhsMat is the right-hand side of Equation 28 in [1].
            rhsMat=PHT'*(FPow(:,:,k)'\H')-H*FPow(:,:,k+1)*PHT;
            for j=1:k
                rhsMat=rhsMat-H*FPow(:,:,j-1+1)*Omega/FPow(:,:,j+1)'*H';
            end
            RhsMatEq28(startIdx:(startIdx+zDim^2-1))=rhsMat(:);
            startIdx=startIdx+zDim^2;
        end
        
        %Hold the left-hand side of Equation 28 in [1] without Q
        %included (Q is being estimated),
        LHSMatEq28=zeros(maxLag*zDim^2,xDim^2);
        startIdx=1;
        for k=1:maxLag
            lhsMat=zeros(zDim^2,xDim^2);
            for j=0:k-1
                lhsMat=lhsMat+kron(H/FPow(:,:,k-j+1),H*FPow(:,:,j+1));
            end
            LHSMatEq28(startIdx:(startIdx+zDim^2-1),:)=lhsMat;
            startIdx=startIdx+zDim^2;
        end

        Q=reshape(linsolve(LHSMatEq28,RhsMatEq28),xDim,xDim);
        Q=(Q+Q')/2;%Force symmetry.

        if(any(~isfinite(Q(:))))
            %Finite precision (or observability) isues arose.
            Q=[];
            R=[];
            return;
        end

        %Force it to be a valid covariance value (no negative eigenvalues).
        [V,D]=eig(Q);
        Q=V*abs(D)*V';
    end
end
end

function innov=getInnovSequence(x0,z,H,F,K,u,num2Skip)
%%GETINNOVSEQUENCE This function just goes through Equations 3.1a and 3.1b
%               in [1] to generate a sequence of innovation terms over
%               time givena fixed Kalman gain.
%
%REFERENCES:
%[2] J. Duník, O. Straka, O. Kost, and J. Havlík, "Noise covariance
%    matrices in state-space models: A survey and comparison of estimation
%    methods-part I," International Journal of Adaptive Control and Signal
%    Processing, vol. 31, no. 11, pp. 1505-1543, Nov. 2017.
%
%July 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

zDim=size(z,1);
N=size(z,2);

%Nominal initial estimate.
x=x0;

%Get the innovation values through the filter.
innov=zeros(zDim,num2Skip);
for k=1:N
    %The measurement prediction after Eq. 3.3 in [2].
    zPred=H*x;
    curInnov=z(:,k)-zPred;
    %Only start saving the innovations after some initial number of steps.
    if(k>num2Skip)
        innov(:,k-num2Skip)=curInnov;
    end

    %State update in Equation 3.1a.
    x=x+K*curInnov;
    if(k<N)
        %State prediction in Equation 3.1b.
        x=F*x+u(:,k);
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
