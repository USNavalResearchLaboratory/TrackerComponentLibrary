function [x,t]=buildBrownianBridge(x1,x2,deltaT,N)
%%BUILDBROWNIANBRIDGE Given the values of a Wiener process at times
%           separated by deltaT, generate N samples of a Wiener process
%           that goes between those points. Wiener processes are also known
%           as Bronian motion. The distribution of a Wiener process
%           conditioned on two endpoints is known as a Brownian bridge.
%
%INPUTS: x1 The numPointsX1 set of start values of the numPoints
%           independent Wiener processes being considered.
%        x2 The numPointsX1 set of end values of the Wiener processes.
%    deltaT The time between the Wiener processes. Wiener processes are
%           defined such that the distribution of w_2-w_1, which are two
%           values of the Wiener process sampled by deltaT, is normal with
%           zero mean and variance sqrt(deltaT).
%         N The number of points to generate between the endpoints; N>=0.
%           If omitted or an empty matrix is passed, then the default of
%           100 is used.
%
%OUTPUTS: x The numPointsX(N+2) points of the Brownian bridges for the
%           numPoints independent processes. The first point is always x1
%           and the last is always x2.
%         t The (N+2)X1 set of time offsets (starting at 0) where the
%           values in x are generated.
%
%A scalar Wiener process is a Markov process such that
%W_{k+1}=W_{k}+sqrt(deltaT)*n, where n is a normal random variable with
%mean zero and variance 1. A vector process is just a bunch of parallel
%scalar processes. The same applied when considering 3 points:
%W_{k+1}=W_{k}+sqrt(deltaT1)*n
%W_{k+2}=W_{k+1}+sqrt(deltaT2)*n
%where deltaT1 and deltaT2 are the time steps going from W_k to W_{k+1} and
%from W_{k+1} to W_{k+2}. So, here, we are given the values of W_{k} and
%W_{k+2} and building the Brownian bridge involves finding the conditional
%distribution of W_{k+1} and sampling it. Thus, assuming that the value at
%time 0 is deterministic, W_{k}, W_{k+1} and W_{k+2} are jointly Gaussian
%with mean 0 and covariance matrix Sigma:
%Sigma=[T0,T0,T0;
%       T0,T1,T1;
%       T0,T2,T2]
%The distribution of W_{k+1} given W_{k} and W_{k+2} can thus be shown to
%be Gaussian with mean mu and variance sigma2 given by:
%mu=((T3-T2)*W_{k}+(T2-T1)*W_{k+2})/(T3-T1)
%sigma2=(T3-T2)*(T2-T1)/(T3-T1)
%Note that the actual values of T1,T2, and T3, do not matter, only the
%relative differences between them. Thus, adapting it for the formulation
%at hand, we have
%mu=(deltaT2*W_{k}+deltaT1*W_{k+2})/(deltaT1+deltaT2);
%sigma2=deltaT1*deltaT2/(deltaT1+deltaT2)
%With deltaT=deltaT1+deltaT2, we thus know the distribution to sample from
%to generate intermediate points. After each point is generated, it becomes
%a new endpoint for the next sample over a smaller region.
%
%EXAMPLE:
%Here, we display 100 Brownian bridges that connect 4 and 10 over a
%deltaT=10.
% numRuns=100;
% x1=4*ones(numRuns,1);
% x2=10*ones(numRuns,1);
% deltaT=10;
% N=500;
% [x,t]=buildBrownianBridge(x1,x2,deltaT,N);
% 
% figure(1)
% clf
% plot(t,x)
% h1=xlabel('t');
% h2=ylabel('x');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(N))
    N=100;
end

xDim=size(x1,1);
x=zeros(xDim,N+2);

%The increment between the sample points.
deltaTInc=deltaT/(N+1);

x(:,1)=x1;
x(:,N+2)=x2;

for curPoint=1:N
    deltaT0=deltaTInc*(curPoint-1);
    deltaTSum=deltaT-deltaT0;
    deltaT2=deltaTSum-deltaTInc;

    mu=(deltaT2*x(:,curPoint)+deltaTInc*x2)/deltaTSum;
    sigma2=deltaTInc*deltaT2/deltaTSum;
    
    x(:,curPoint+1)=mu+sqrt(sigma2)*randn(xDim,1);
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
