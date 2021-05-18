function [x,t]=randLinTrajConstEnds(x0,x2,deltaT,q0,order,N,numTraj,FFun,QFun)
%%RANDLINTRAJCONSTENDS Generate a random trajectory based on a linear
%               dynamic model with additive noise, as is typically used in
%               tracking systems, given known locations of the start and
%               end points. A polynomial dynamic model is used by default,
%               but one can specify other models. This is essentially a
%               generalization of a Brownian bridge.
%
%INPUTS: x0, x2 The xDimX1 target states at the beginning and end of the
%               region over which paths are to be simulated. These are
%               treated as known values.
%        deltaT The time between when state x0 occurred and when state x2
%               occurred. This is a positive value.
%            q0 The power spectral density of the process noise. This is
%               passed to QFun, as described below and by default would
%               corresponds to the same-named input in QPolyKal. If a
%               custom QFun is used and this parameter is not needed, then
%               an empty matrix can be passed.
%         order The order >=0 of the filter. If order=1, then it is
%               constant velocity, 2 means constant acceleration, 3 means
%               constant jerk, etc. This is passed to FFun and QFun. If
%               custom FFun and QFun functions are used and this is not
%               needed, then an empty matrix can be passed. The default if
%               omitted or an empty matrix is passed is 1.
%             N The number of points in the random path to generate between
%               x0 and x2. The default if omitted or an empty matrix is
%               passed is 100; N>=0.
%       numTraj The number of random trajectories to generate. The default
%               if omitted or an empty matrix is passed is 1.
%    FFun, QFun Optional function handles to functions that respectively
%               provide the state transition matrix and the process noise
%               covariance matrix as a function of the time duration over
%               which the prediction is taken. The defaults if omitted or
%               empty matrices are passed are handles to FPolyKal and
%               QPolyKal. The functions are called as FFun(T,xDim,order)
%               and QFun(T,xDim,order,q0), where T is the duration over
%               which the prediction is taken.
%
%OUTPUTS: x An xDimX(N+2)XnumTraj set of random trajectories generated
%           according to the specified model. For the ith trajectory,
%           x(:,1,i) is always x0 and x(:,N+2,i) is always x2.
%         t The (N+2)X1 set of time offsets (starting at 0) where the
%           values in x are generated.
%
%Consider a process at points 0, 1, and 2, where we are given the value at
%time x0. The transition models are:
%x1=F0*x0+w0
%x2=F1*x1+w1
%Where x0 is given (deterministic) and w0 and w1 are independent zero-mean
%Gaussian random variables with covariance matrices Q0 and Q1.
%x1 and x2 are jointly normal with mean mu and covariance matrix Sigma:
%mu=[x1Pred; = [F0*x0
%    x2Pred]    F1*F0*x0];
%Sigma=[   Q0,        Q0*F1'; =[Pxx, Pxz;
%       F1*Q0,  F1*Q0*F1'+Q1]   Pxz',Pzz];
%When given the value of x2, the conditional distribution of x1|x0,x2 is
%the same as the Kalman filter update taking the measurement covariance
%matrix to be Q1, the measurement matrix H to be F1, the prior state to be
%F0*x0, and the covariance of the prior state to be Q0. Thus, as derived,
%for example, in [Ch. 5,1], [2], the posterior distribution of x1 is
%Gaussian with mean and covariance matrix
%x1Bar=x1Pred+Pxz/Pzz*(x2-x2Pred)  Px1=Pxx-Pxz*inv(Pzz)*Pxz';
%Given the similarity to the Kalman filter, the KalmanFilterUpdate function
%is used, which puts the covariance update in Joseph form. Given the
%conditional distribution, one can generate a random sample of a point
%between x1 and x2. Then, that point can become a new endpoint and another
%sample can be generated, etc.
%
%EXAMPLE:
%This plots a number of random trajectories. All of the trajectories have
%to eventually turn around, since the final point is not in a direction
%close to the original velocity direction.
% q2=14;
% order=1;
% x0=[1e3;0;30;20];
% x2=[709;
%     527;
%     -39;
%       4];
% N=500;
% numTraj=25;
% deltaT=100;
% x=randLinTrajConstEnds(x0,x2,deltaT,q2,order,N,numTraj);
% figure(1)
% clf
% hold on
% for curTraj=1:numTraj
%     plot(x(1,:,curTraj),x(2,:,curTraj))
% end
% h1=xlabel('x');
% h2=ylabel('y');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[2] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(numTraj))
    numTraj=1;
end

if(nargin<6||isempty(N))
    N=100;
end

if(nargin<5||isempty(order))
    order=1;
end

xDim=size(x0,1);

if(nargin<9||isempty(QFun))
    Q=@(T)QPolyKal(T,xDim,order,q0);
else
    Q=@(T)QFun(T,xDim,order,q);
end

if(nargin<8||isempty(FFun))
    F=@(T)FPolyKal(T,xDim,order);
else
    F=@(T)FFun(T,xDim,order);
end

x=zeros(xDim,N+2,numTraj);

%The increment between the sample points.
deltaTInc=deltaT/(N+1);

x(:,1,:)=repmat(x0,[1,1,numTraj]);
x(:,N+2,:)=repmat(x2,[1,1,numTraj]);

F0=F(deltaTInc);
Q1=Q(deltaTInc);
for curTraj=1:numTraj
    for curPoint=1:N
        deltaT0=deltaTInc*(curPoint-1);
        deltaT2=deltaT-deltaT0-deltaTInc;

        F1=F(deltaT2);
        Q2=Q(deltaT2);

        [x1Mean,P1Mean]=KalmanUpdate(F0*x(:,curPoint,curTraj),Q1,x2,Q2,F1);

        x(:,curPoint+1,curTraj)=x1Mean+chol(P1Mean,'lower')*randn(xDim,1);
    end
end

if(nargout>1)
    t=[0;deltaTInc*(1:(N+1)).'];
    t(end)=deltaT;%Make it exact despite finite precision limitations.
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
