function [xPred,PPred]=BlackScholesPred(x,a,D,deltaT)
%%BLACKSCHOLESPRED A multivariate geometric Brownian motion model (a Black-
%           Scholes model) is given by the stochastic differential equation
%           dx(i)=a(i)*x(i) dt+sum_{j=1}^m D(i,j)x(i) dW(j)
%           where we are considering the ith component of x and dW is the
%           differential of a Wiener process. This function provides the
%           mean and covariance matrix of the prediction of the function a
%           time-step of deltaT into the future when given a deterministic
%           initival value.
%
%INPUTS: x The dX1 initial state.
%        a The dX1 drift coefficient.
%        D The dXm diffusion coefficient matrix.
%   deltaT The time duration of the step; the step size.
%
%OUTPUTS: xPred The dX1 mean of the process deltaT in the future.
%         PPred The dXd covariance matrix of the process deltaT in the
%               future.
%
%The solution is given in Equations 108 and 109 of [1]. However, the claim
%in [1] that D must be square is incorrect.
%
%EXAMPLE:
%Here we compare the prediction given by this function to the mean and
%covariance matrix obtained by randomly simulations using the function
%BlackScholesStep.
% x0=[1;2;3;1];
% a=[-0.9365;
%    -1.7956;
%     1.3305;
%    -0.1772];
% D=[1.6154, -0.0284, -0.1684;
%    0.1034,  0.4361, -1.3050;
%    0.9386,  0.0641,  0.7215;
%   -1.1955,  0.4186,  1.0875];
% d=size(D,1);
% deltaT=1/8;
% [xPred,PPred]=BlackScholesPred(x0,a,D,deltaT)
% 
% numRuns=1e5;
% x=zeros(d,numRuns);
% for curRun=1:numRuns
%     x(:,curRun)=BlackScholesStep(x0,a,D,deltaT);
% end
% [xSim,PSim]=calcMixtureMoments(x)
%One will typically see that the results agree to about one or two decimal
%places.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, pp. 4-41, Feb. 2015.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    d=size(a,1);

    mu=(a-(1/2)*diag(D*D'))*deltaT;
    Sigma=D*D'*deltaT;
    
    xPred=diag(x)*exp(mu+(1/2)*diag(Sigma));
    
    if(nargout>1)
        PPred=zeros(d,d);
        for i=1:d
            for j=i:d 
                PPred(i,j)=x(i)*x(j)*exp(mu(i)+mu(j)+(1/2)*(Sigma(i,i)+Sigma(j,j)))*(exp(Sigma(i,j))-1);
                PPred(j,i)=PPred(i,j);
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
