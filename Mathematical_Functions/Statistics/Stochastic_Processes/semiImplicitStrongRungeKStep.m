function [yNext,didConverge]=semiImplicitStrongRungeKStep(y,t,a,B,deltaT,algorithm,aCur,BCur,useNewton,maxIter,theta,RelTol,AbsTol,p,ItoAlg)
%%SEMIIMPLICITSTRINGRUNGEKSTEP Perform a single step of a semi-implicit
%           strong stochastic Runge-Kutta method under Itô calculus. This
%           integrates a d-dimensional stochastic differential equation of
%           the form
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process. As noted in Chapter 9.8 of [1], fully implicit strong
%           Itô-Taylor expansions are unstable. Thus, fully
%           implicit strong stochastic Runge Kutta methods would also be
%           unstable. Thus, semi-implicit techniques are typically used.
%           Strong methods converge to an optimal path as the stepsize
%           decreases, not considering finite precision errors. Runge-Kutta
%           methods use finite differencing instead of explicit
%           derivatives.
%
%INPUTS: y The dX1 initial value of the random process.
%        t The scalar initial time of the random process. If an empty
%          matrix is passed, t=0 is used.
%        a A function handle to the drift function. This is called as
%          a(y,t) and returns a dX1 vector. If one wishes to use Newton's
%          method for the implicit iteration, then the calling format is
%          [aVal,papy]=a(y,t), where papy is the dXd matrix of partial
%          derivatives of a with respect to the elements of y papy(:,i) is
%          the derivative with respect to the ith component of y.
%        B A function handle to the diffusion matrix function. This is
%          called as B(y,t) and returns a dXm matrix.
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. Possible values
%          are:
%          -1 Use the semi-implicit Euler-Maruyama method (order 0.5) of
%             Equation 2.2 of Chapter 12.2 of [1]. This is the semi-
%             implicit version of Equation 2.1 in Chapter 10.2 of [1].
%           0 Semi-implicit order 1.0 method for scalar noise (generalizing
%             Equation 3.1 in Chapter 12.3 of [1] to multiple dimensions).
%             This is the semi-implicit version of Equation 1.5 in Chapter
%             11.1 of [1]. It has been generalized in the same manner as
%             Equation 3.2 of Chapter 12.3 of [1] to allow for theta values
%             other than 1. This requires that m=1.
%           1 (The default if omitted or an empty matrix is passed) Use the
%             semi-implicit order 1.0 method for general noise from
%             Equation 3.2 of Chapter 12.3 of [1]. This is the semi-
%             implicit version of Equation 1.7 in Chapter 11.1 of [1].
%           2 Use the semi-implicit order 1.0 method for diagonal noise.
%             This is a simplification of Equation 3.2 of Chapter 12.3 of
%             [1] in the same manner that Equation 1.9 of Chapter 11.1 for
%             the explicit diagonal form is a simplification of Equation
%             1.7 in Chapter 11.1.
%           3 Use the semi-implicit order 1.5 method for autonomous scalar
%             problems from Equation 3.4  in Chapter 12.3 of [1]. This is
%             the semi-implicit version of Equation 2.1 in Chapter 11.2 of
%             [1]. This algorithm does not use the input theta; the
%             equivalent of theta is fixed at 1/2.
%           4 Use the semi-implicit order 1.5 method for non-autonomous
%             additive noise from Equation 3.6 of Chapter 12.3 of [1]. This
%             is the semi-implicit version of Equation 2.19 in Chapter 11.2
%             of [1]. This algorithm does not use the input theta; the
%             equivalent of theta is fixed at 1/2.
%           5 Use the semi-implicit order 1.5 method for autonomous
%             additive noise from Equation 3.6 of Chapter 12.3 of [1],
%             removing the relevant B terms. This is the semi-implicit
%             version of Equation 2.19 in Chapter 11.2 of [1] removing the
%             relevant B terms.
% aCur, BCur Often one might already have the values a(x,t) and B(x,t). If
%          so, then they should be provided as the dX1 and dXm aCur and
%          BCur to avoid recalculation. If unavailable, these values can be
%          omitted or empty matrices passed.
% useNewton Indicates whether the implicit iteration should be performed
%          using Newton's method or fixed-point iteration. The default if
%          omitted or an empty matrix is passed is false (fixed-point
%          iteration). If Newton's method is used, then papy must be
%          returned by the function a.
%  maxIter The maximum number of iterations to perform. The default if
%          omitted or an empty matrix is passed is 2.
%    theta A dX1 parameter that affects the degree of impliciteness of the
%          method in each dimension. This can be scalar if the value is the
%          same in all dimensions. This varies from 0 to 1. The default if
%          omitted or an empty matrix is passed is 1/2. 0 means an explicit
%          step and 1 is the most implicit. This parameter is ignored for
%          algorithm=3,4,5, as it is fixed at 1/2 for those algorithms.
% RelTol, AbsTol The relative and absolute tolerances on the iterations
%          before declaring convergence. If these are set to 0 (the default
%          if omitted or empty matrices are passed), then the algorithm
%          will just iterate for the maximum number of iterations. the
%          tolerances apply to each element of x. Convergence is declared
%          if all(diff<=AbsTol)||all(diff<=RelTol*abs(yNext)).
% p, ItoAlg These values are only used if  algorithm=1. These values
%          correspond to the inputs p and algorithm of the function
%          sampleItoIntegrals, which is called by this function for
%          necessary terms. Defaults if omitted or an empty matrix is
%          passed are 15 and 0.
%
%OUTPUTS: yNext The estimated value of the process after taking a step of
%               deltaT. This is a random value.
%   didConverge If RelTol and/or AbsTol are not zero and maxIter>0, then
%               this indicates whether the iterations converged to the
%               desired accuracy. Otherwise, this is just an empty matrix.
%
%EXAMPLE 1:
%Here, we compare the implicit strong Taylor method to the explicit method
%on a linear model. This is done with the same noise driving both
%processes.
% algorithm=5;
% useNewton=false;
% numMC=1e4;
% deltaT=1/3;
% numSteps=5;
% y0=[1/4;-12];
% A=[1.1,0.1;
%    -0.2,2.2];
% D=[1.5,-0.4;
%    0.1,1];
% d=size(D,1);
% [F,Q]=linDynMod2Disc(deltaT,A,D);
% mu=F*y0;
% P=Q;
% aDrift=@(x,t)(A*x);
% BDiff=@(x,t)(D);
% 
% aFun=@(y,t)dealRobust(aDrift(y,t),A);
% 
% valsRKImp=zeros(d,numMC);
% valsRK=zeros(d,numMC);
% for curMC=1:numMC
%     yExpl=y0;
%     y=y0;
%     t=0;
%     for curStep=1:numSteps
%         s=rng();%Record the state prior to generating the random variables.
%         y=semiImplicitStrongRungeKStep(y,t,aFun,BDiff,deltaT/numSteps,algorithm,[],[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         yExpl=strongRungeKStep(yExpl,t,aDrift,BDiff,deltaT/numSteps,algorithm);
%         t=t+deltaT/numSteps;
%     end
%     valsRKImp(:,curMC)=y;
%     valsRK(:,curMC)=yExpl;
% end
% [muRKImp,PRKImp]=calcMixtureMoments(valsRKImp);
% [muRK,PRK]=calcMixtureMoments(valsRK);
% norm(muRKImp-mu)./norm(mu)%Relative mean error, implicit.
% norm(PRKImp-P,'fro')./norm(P,'fro')%Relative variance error,implicit.
% norm(muRK-mu)./norm(mu)%Relative mean error, explicit.
% norm(PRK-P,'fro')./norm(P,'fro')%Relative variance error,explicit.
%One will typically see that the relative errors of the implicit method are
%notably better than the explicit method.
%
%EXAMPLE 2:
%This is an example of what is considered a "stiff" problem in Section 12.2
%of [1] with non-additive noise. The problem has an explicit solution in
%terms of the Wiener process W. We use quadrature integration to obtain the
%mean of that solution, which is used as the true moments for comparison.
%We compare this implicit method with an explicit step. This is a 2X1
%process with scalar noise.
% rng(0)%Make the exact run repeatable.
% algorithm=1;
% useNewton=true;
% numMC=1e3;
% a=25;
% b=2;
% deltaT=1/10;
% numSteps=10;
% A=[-a,a;
%     a,-a];
% B=[b,0;
%    0,b];
% aDrift=@(y,t)A*y;
% BDiff=@(y,t)B*y;
% d=2;
% papy=A;
% aFun=@(y,t)dealRobust(aDrift(y,t),papy);
% %The explicit solution; equation 2.5 in Chapter 12.2 of [1].
% expSol=@(W,deltaT,y0)(expm((A-(1/2)*B^2)*deltaT+B*W)*y0);
% y0=[1;0.1];
% 
% %Take the expected value of the explicit solution using quadrature
% %integration.
% [xi,w]=quadraturePoints1D(6);%2*6-1=11th order.
% numPts=length(w);
% xi=sqrt(deltaT)*xi;
% muCub=zeros(d,1);
% for k=1:numPts
%     muCub=muCub+w(k)*expSol(xi(:,k),deltaT,y0);
% end
% 
% valsRK=zeros(d,numMC);
% valsRKExplicit=zeros(d,numMC);
% for curMC=1:numMC
%     y=y0;
%     yExp=y0;
%     t=0;
%     for curStep=1:numSteps
%         s=rng();%Record the state prior to generating the random variables.
%         y=semiImplicitStrongRungeKStep(y,t,aFun,BDiff,deltaT/numSteps,algorithm,[],[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         yExp=strongRungeKStep(yExp,t,aDrift,BDiff,deltaT/numSteps,algorithm);
%         t=t+deltaT/numSteps;
%     end
%     valsRK(:,curMC)=y;
%     valsRKExplicit(:,curMC)=yExp;
% end
% muRK=mean(valsRK,2);
% muRKExp=mean(valsRKExplicit,2);
% norm((muRK-muCub)./norm(muCub))%Relative mean error, implicit.
% norm((muRKExp-muCub)./norm(muCub))%Relative mean error, explicit.
%One will get a relative implicit error of about 0.0167 and a relative
%explicit error of about 0.0173. Getting rid of the random seed, the
%implicit solution is typically more accurate than the explicit. 
%
%EXAMPLE 3:
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its solution
%are from Equation 4.40. We compare the performance of the explicit
%solution (taking the expected value using cubature integration) and the
%implicit solution.
% rng(1)%Make exact run repeatable.
% numMC=1e3;
% numSteps=1;
% algorithm=0;
% useNewton=true;
% deltaT=1.1;
% y0=0.5;
% aDrift=@(y,t)((1/3)*y^(1/3));
% BDiff=@(y,t)(y^(2/3));
% papy=@(y)1/(9*y^(2/3));
% aFun=@(y,t)dealRobust(aDrift(y,t),papy(y));
% explSim=@(W)(y0^(1/3)+(1/3)*W)^3;
% 
% %Take the expected value of the explicit solution using quadrature
% %integration.
% [xi,w]=quadraturePoints1D(6);%2*6-1=11th order.
% numPts=length(w);
% xi=sqrt(deltaT)*xi;
% muCub=0;
% for k=1:numPts
%     muCub=muCub+w(k)*explSim(xi(:,k));
% end
% 
% valsRKImp=zeros(1,numMC);
% valsRK=zeros(1,numMC);
% for curMC=1:numMC
%     y=y0;
%     yExp=y0;
%     t=0;
%     for curStep=1:numSteps
%         s=rng();%Record the state prior to generating the random variables.
%         y=semiImplicitStrongRungeKStep(y,t,aFun,BDiff,deltaT/numSteps,algorithm,[],[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         yExp=strongRungeKStep(yExp,t,aDrift,BDiff,deltaT/numSteps,algorithm);
%         t=t+deltaT/numSteps;
%     end
%     valsRKImp(curMC)=y;
%     valsRK(curMC)=yExp;
% end
% muRKImp=mean(valsRKImp);
% muRK=mean(valsRK);
% abs((muRKImp-muCub)./muCub)%Relative mean error, implicit.
% abs((muRK-muCub)./muCub)%Relative mean error, explicit.
%The implicit error will be about 0.0072 and the explicit error will be
%about 0.0127. Thus, the algorithm improves the explicit method in this
%instance.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<15||isempty(ItoAlg))
    ItoAlg=0;
end

if(nargin<14||isempty(p))
    p=15;
end

if(nargin<13||isempty(AbsTol))
    AbsTol=0;
end

if(nargin<12||isempty(RelTol))
    RelTol=0; 
end

if(nargin<11||isempty(theta))
    theta=1/2;
end

if(nargin<10||isempty(maxIter))
    maxIter=2;
end

if(nargin<9||isempty(useNewton))
    useNewton=false;
end

if(nargin<8||isempty(BCur))
    BCur=B(y,t); 
end

if(nargin<7||isempty(aCur))
    aCur=a(y,t); 
end

if(nargin<6||isempty(algorithm))
    algorithm=1;
end

if(isempty(t))
    t=0; 
end

if(algorithm>=-1&&algorithm<=2)%0, 1 or 2
    [yNext,endTerm]=strongRungeKStep(y,t,a,B,deltaT,algorithm,aCur,BCur,p,ItoAlg);

    fixedTerm=(1-theta).*aCur*deltaT+endTerm;
    tNext=t+deltaT;
    if(useNewton)
        [yNext,didConverge]=implicitNewtonIter(y,a,deltaT*theta,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
    else
        [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT*theta,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
    end
elseif(algorithm<=5)
    [yNext,fixedTerm]=strongRungeKStep(y,t,a,B,deltaT,algorithm,aCur,BCur,p,ItoAlg);
    tNext=t+deltaT;
    if(useNewton)
        [yNext,didConverge]=implicitNewtonIter(y,a,deltaT*(1/2),yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
    else
        [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT*(1/2),yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol); 
    end
else
    error('Unknown algorithm specified.')
end
end

function [yNext,didConverge]=implicitFixedPointIter(yOld,a,coeff,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol)
%IMPLICITFIXEDITER This iterates the function
%                  yNext=a(yNext,tNext)*coeff+fixedTerm
%                  Iterations continue until the RelTol and or AbsTol
%                  conditions are met or maxIter iterations have elapsed.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(RelTol==0&&AbsTol==0)
        %Fixed point iteration for the maximum number of steps.
        for curIter=1:maxIter
            yNext=a(yNext,tNext)*coeff+fixedTerm;
        end
        didConverge=[];
    else%Iterate until meeting the convergence bound.
        didConverge=false;
        for curIter=1:maxIter
            yNext=a(yNext,tNext)*coeff+fixedTerm;

            diff=abs(yOld-yNext);
            if(all(diff<=AbsTol)||all(diff<=RelTol*abs(yNext)))
                didConverge=true;
                break; 
            end

            yOld=yNext;
        end
    end
end

function [yNext,didConverge]=implicitNewtonIter(yOld,a,coeff,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol)
%%IMPLICITNEWTONITER This function uses Newton's method to try to solve
%                    a(yNext,tNext)*coeff+fixedTerm-yNext=0
%                    for yNext.Iterations continue until the RelTol and or
%                    AbsTol conditions are met or maxIter iterations have
%                    elapsed.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    d=size(yNext,1);
    
    I=eye(d,d);
    if(RelTol==0&&AbsTol==0)
        %Newton's iteration for the maximum number of steps.
        for curIter=1:maxIter
            [aVal,papy]=a(yNext,tNext);
            F=yNext-fixedTerm-coeff*aVal;
            dF=I-coeff*papy;
            yNext=yNext-dF\F;
        end
        didConverge=[];
    else%Iterate until meeting the convergence bound.
        didConverge=false;
        for curIter=1:maxIter
            [aVal,papy]=a(yNext,tNext);
            F=aVal-fixedTerm-coeff*aVal;
            dF=I-coeff*papy;
            yNext=yNext-dF\F;

            diff=abs(yOld-yNext);
            if(all(diff<=AbsTol)||all(diff<=RelTol*abs(yNext)))
                didConverge=true;
                break; 
            end

            yOld=yNext;
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
