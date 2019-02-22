function [xPred,PPred]=BlackScholesPredGaussPrior(xPrior,PPrior,a,D,deltaT)
%%BLACKSCHOLESPREDGAUSSPRIOR A multivariate geometric Brownian motion model
%           (a Black-Scholes model) is given by the stochastic differential
%           equation
%           dx(i)=a(i)*x(i) dt+sum_{j=1}^m D(i,j)x(i) dW(j)
%           where we are considering the ith component of x and dW is the
%           differential of a Wiener process. This function provides the
%           mean and covariance matrix of the prediction of the function a
%           time-step of deltaT into the future. when the prior is
%           Gaussian.
%
%INPUTS: xPrior The dX1 mean of the distribution prior to propagation.
%        PPrior The dXd covariance matrix of the distribution prior to
%               propagation. This can be singular.
%        a The dX1 drift coefficient.
%        D The dXm diffusion coefficient matrix.
%   deltaT The time duration of the step; the step size.
%
%OUTPUTS: xPred The dX1 mean of the process deltaT in the future.
%         PPred The dXd covariance matrix of the process deltaT in the
%               future.
%
%The solution for propagation given a deterministic prior is given in
%Equations 108 and 109 of [1]. Given a Gaussian prior, the mean is trivial
%to find. The covariance is obtained using the Law of Total Covariance.
%That is that Cov(X)=E{Cov(X|Y)}+Cov(E{X|Y}) Here, X is the predicted value
%and Y is the prior value.
%
%EXAMPLE:
%We demonstrate that the value obtained here is consistent with computing
%the mean of BlackScholesPred and taking the total covariance of those
%random runs. This can be a bit slow.
% numMC=1e6;%Number Monte Carlo Runs
% deltaT=1/5;
% xPrev=[1;2;3;1];
% SPrev=diag([4;1;2;0.5]);
% PPrev=SPrev*SPrev';
% a=[0.9;
%    1.7;
%    1.3;
%    0.1];
% D=[1.6154,  0.0284;
%    0.1034,  0.4361;
%    0.9386,  0.0641;
%    1.1955,  0.4186];%d=4,m=2.
% d=size(D,1);
% 
% xPredMC=zeros(d,numMC);
% PPredMC=zeros(d,d);
% for k=1:numMC
%     x=xPrev+SPrev*randn(d,1);
%     [xPred,PPred]=BlackScholesPred(x,a,D,deltaT);
%     xPredMC(:,k)=xPred;
%     PPredMC=PPredMC+PPred;
% end
% [xPredMC,PPred0]=calcMixtureMoments(xPredMC);
% %Law of total covariance.
% PPredMC=PPred0+PPredMC/numMC;
% 
% [xPred,PPred]=BlackScholesPredGaussPrior(xPrev,PPrev,a,D,deltaT);
% max(max(abs((xPred-xPredMC)./xPredMC)))%Maximum relative mean error
% max(max(abs((PPred-PPredMC)./PPredMC)))%Maximum relative covariance error
%The mean error is typically better than 1%. The covariance error is
%typically better than 5%.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, pp. 4-41, Feb. 2015.
%
%January 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    d=size(a,1);

    mu=(a-(1/2)*diag(D*D'))*deltaT;
    Sigma=D*D'*deltaT;
    
    const=diag(exp(mu+(1/2)*diag(Sigma)));
    xPred=const*xPrior;
    
    if(nargout>1)
        PPred=zeros(d,d);
        for i=1:d
            for j=i:d 
                %The expected value E{x(i)*x(j)}
                mij=PPrior(i,j)+xPrior(i)*xPrior(j);

                PPred(i,j)=mij*exp(mu(i)+mu(j)+(1/2)*(Sigma(i,i)+Sigma(j,j)))*(exp(Sigma(i,j))-1);
                PPred(j,i)=PPred(i,j);
            end
        end
        
        %Law of Total Covariance
        PPred=PPred+const*PPrior*const';
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
