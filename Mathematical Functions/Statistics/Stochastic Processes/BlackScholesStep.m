function xT=BlackScholesStep(x,a,D,deltaT)
%%BLACKSCHOLESSTEP A multivariate geometric Brownian motion model (a Black-
%           Scholes model) is given by the stochastic differential equation
%           dx(i)=a(i)*x(i) dt+sum_{j=1}^m D(i,j)x(i) dw(j)
%           where we are considering the ith component of x and dw is the
%           differential of a Wiener process. This function simulates a
%           random step forward a time increment of deltaT when a and B are
%           constants. The result is random, because the Wiener process is
%           random. 
%
%INPUTS: x The dX1 initial state.
%        a The dX1 drift coefficient.
%        D The dXm diffusion coefficient matrix.
%   deltaT The time duration of the step; the step size.
%
%OUTPUTS: xT A dX1 random sample of x at time deltaT forward.
%
%The solution is given in Equation 99 of [1], and is taken from Chapter 2.4
%of [2] and is also in Chapter 3.2.3 of [3]. Note that Equation 98 in [1]
%contains a typo in the sum for the noise term of the stochastic
%differential equation. The sum in the description above is correct.
%
%EXAMPLE:
%Here, we compare the result of this explicit solution to integrating the
%stochastic differential equation using the Euler-Maruyama method and a
%small step size.
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
% deltaT=1/3;
% 
% numRuns=1e3;
% x=zeros(d,numRuns);
% xMC=zeros(d,numRuns);
% numSimSteps=500;
% deltaTSim=deltaT/numSimSteps;
% aFun=@(x,t)aGeoBrownian(x,a);
% BFun=@(x,t)DGeoBrownian(x,D);
% parfor curRun=1:numRuns
%     x(:,curRun)=BlackScholesStep(x0,a,D,deltaT);
%     
%     xCur=x0;
%     for curStep=1:numSimSteps
%         xCur=strongRungeKStep(xCur,0,aFun,BFun,deltaTSim);
%     end
%     xMC(:,curRun)=xCur;
% end
% [xStep,PStep]=calcMixtureMoments(x)
% [xMC,PMC]=calcMixtureMoments(xMC)
%There will typically be a digit or more of agreement between the values in
%xStep, PStep and those in xMC, PMC. Better agreement is obtained with more
%Monte Carlo runs.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear continuous-time dynamic
%    models," IEEE Aerospace and Electronic Systems Magazine, vol. 30, no.
%    2, pp. 4-41, Feb. 2015.
%[2] E. Platen and N. Bruti-Liberati, Numerical Solution of Stochastic
%    Differential Equations with Jumps in Finance. Berlin: Springer-Verlag,
%    2010.
%[3] P. Glasserman, Monte Carlo Methods in Financial Engineering. New
%   York: Springer, 2004.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(D,2);

deltaW=sqrt(deltaT)*randn(m,1);
xT=diag(x)*exp((a-(1/2)*diag(D*D'))*deltaT+D*deltaW);

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
