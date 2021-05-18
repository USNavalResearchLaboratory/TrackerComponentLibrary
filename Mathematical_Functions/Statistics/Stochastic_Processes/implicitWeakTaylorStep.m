function [yNext,didConverge]=implicitWeakTaylorStep(y,curT,a,B,deltaT,algorithm,givenVals,useNewton,maxIter,thetaVals,RelTol,AbsTol)
%%IMPLICITWEAKTAYLORSTEP Perform a single step of an implicit or
%           semi-implicit weak -Itô-Taylor expansion to integrate a d-
%           dimensional stochastic differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process and t is time. As noted in Chapter 9.8 of [1], fully
%           implicit strong Itô-Taylor expansion expansions are unstable,
%           because particular values of the noise process can make them
%           diverge. However, this issue can often be avoided with weak
%           Itô-Taylor expansions. As the stepsize used decreases, weak
%           methods converge such that integrals with the random process
%           are a measure are correct (for example, to determine moments).
%           However, they do not converge to the optimal path, unlike
%           strong methods.
%
%INPUTS: y The dX1 initial value of the random process.
%     curT The current time associated with the random process.
%        a A function handle to the drift function. This is called as
%          [aVal,papy,papt,p2apypy]=a(y,t), where aVal is the dX1 value of
%          the drift function, papy and p2apypy are its dXd and dXdXd
%          matrices of first and second partial derivatives with respect to
%          the elements of y; papt is its dX1 partial derivative with
%          respect to time. Not all algorithms require all of the outputs;
%          if an algorithm doesn't use a particular outputs, then an empty
%          matrix can be returned for that output. papt(:,i) is the first
%          partial derivative with respect to y(i) and p2apypy(:,i,j) is
%          the second partial derivatives with respect to y(i) and y(j).
%        B A function handle to the diffusion matrix function. This is
%          called as [BVal,pBpy,pBpt,p2Bpypy]=B(y,t), where BVal is a dXm
%          matrix and pBpy,pBpt,p2Bpypy are first and second partial
%          derivatives analogous to those returned b a. Not all algorithms
%          require all outputs.
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. All algorithms use
%          the first outputs of a and B. Additionally, all methods that
%          support a Newton iteration use papy. Possible values of
%          algorithm are:
%          0 Use an implicit weak-Euler-Maruyama method (order 1.0). This
%            takes theta1 and theta2 as inputs. This is Equation 4.6 in
%            Chapter 15.4 of [1], which is the same as Equation 5.4 in
%            Chapter 15.5. By default theta1=1/2 and theta2=0, which
%            corresponds to the trapezoidal method. theta2 adjusts the
%            implicitness of the diffusion and theta1 adjusts the
%            implicitness of the drift. If theta2 is not zero, then B must
%            be able to return pBpy. Newton's method is not supported if
%            theta2 is not zero.
%          1 Use the implicit order 2.0 weak Taylor scheme for autonomous
%            scalar problems; Equation 5.7 of Chapter 15.5 of [1], which
%            is just Equation 4.11 of Chapter 15.5. The predictor is the
%            explicit order 2.0 weak Taylor scheme in Equation 2.2 of
%            Chapter 14.2, which is just Equation 5.8 in Chapter 15.5. This
%            requires papy,p2apypy,pBpy,p2Bpypy, that a and B not
%            depend on time and that d=m=1. This does not use the theta
%            values.
%          2 Use the implicit order 2.0 weak Taylor scheme for general
%            problems; Equation 5.10 of Chapter 15.5 of [1], which is the
%            same as Equation 4.11 of Chapter 15.4 of [1]. The predictor is
%            the general explicit order 2.0 weak Taylor scheme in Equation
%            2.7 of Chapter 14.2 of [1], which is just Equation 5.11 of
%            Chapter 15.5 of [1]. This requires papy,papt,p2apypy,pBpy,pBpt,
%            and p2Bpypy. This does not use the theta values.
% givenVals If the outputs of a(y,curT) and B(y,curT) (as described above)
%          are available, they can be passed as members of this
%          structure, each having the names used above. For example, papy
%          for the first derivatives a with respect to the y. If omitted or
%          an empty matrix is passed, then the functions a with b will be
%          called for the initial values, as needed. All initially required
%          values must be provided for a or b if any for a or b are
%          provided. For example, providing papy but not aVal will result
%          in an error. However, providing neither will just result in an
%          extra call to a. One can pass values for just a or for just b
%          or for both.
% useNewton Indicates whether the implicit iteration should be performed
%          using Newton's method or fixed-point iteration. The default if
%          omitted or an empty matrix is passed is false (fixed-point
%          iteration). If Newton's method is used, then papy must be
%          returned by the function a.
%  maxIter The maximum number of iterations to perform. The default if
%          omitted or an empty matrix is passed is 2.
% thetaVals Algorithm 0 has the option for parameters theta1 and theta2
%          that affect the degree of impliciteness in the drift and
%          diffusion terms. A value of 1 is full impliciteness and 0 is
%          explicit. thetaVals is a 2X1 or 1X2 vector holding the values
%          of theta1 and theta2  The default if omitted or an empty matrix
%          matrix is passed is [1/2;0]. If a scalar is passed, it is
%          assumed that one only wishes to change theta1.
% RelTol, AbsTol The relative and absolute tolerances on the iterations
%          before declaring convergence. If these are set to 0 (the default
%          if omitted or empty matrices are passed), then the algorithm
%          will just iterate for the maximum number of iterations. The
%          tolerances apply to each element of y. Convergence is declared
%          if all(diff<=AbsTol)||all(diff<=RelTol*abs(yNext)), where diff
%          is the change in y between subsequent steps.
%
%OUTPUTS: yNext The estimated value of the process after taking a step of
%               deltaT. This is a random value.
%   didConverge If RelTol and/or AbsTol are not zero and maxIter>0, then
%               this indicates whether the iterations converged to the
%               desired accuracy. Otherwise, this is just an empty matrix.
%
%EXAMPLE 1:
%This is an example of what is considered a "stiff" problem in Section 12.2
%of [1] with non-additive noise. The problem has an explicit solution in
%terms of the Wiener process W. We use quadrature integration to obtain the
%mean of that solution, which is used as the true moments for comparison.
%We compare this implicit method with an explicit step. This is a 2X1
%process with scalar noise.
% rng(0)%Make the exact run repeatable.
% algorithm=0;
% useNewton=false;
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
% m=1;
% papy=A;
% papt=zeros(d,1);
% p2apypy=zeros(d,d,d);
% pBpy=reshape(B,[d,m,d]);
% pBpt=zeros(d,m);
% p2Bpypy=zeros(d,m,d,d);
% aFun=@(y,t)dealRobust(aDrift(y,t),papy,papt,p2apypy);
% BFun=@(y,t)dealRobust(BDiff(y,t),pBpy,pBpt,p2Bpypy);
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
%         y=implicitWeakTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExp,t);
%         BCur=BDiff(yExp,t);
%         simplified=2;
%         yExp=weakStochTaylorStep(yExp,aCur,BCur,deltaT/numSteps,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt);
%         t=t+deltaT/numSteps;
%     end
%     valsRK(:,curMC)=y;
%     valsRKExplicit(:,curMC)=yExp;
% end
% muRK=mean(valsRK,2);
% muRKExp=mean(valsRKExplicit,2);
% norm((muRK-muCub)./norm(muCub))%Relative mean error, implicit.
% norm((muRKExp-muCub)./norm(muCub))%Relative mean error, explicit.
%One will get a relative implicit error of about 0.0069 and a relative
%explicit error of about 0.0082. Getting rid of the random seed, the
%implicit solution is typically more accurate than the explicit. The same
%relation between explicit and implicit solutions should hold when using
%algorithm 2. Note that we are comparing the norms, not the individual
%components.
%
%EXAMPLE 2:
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
% simplified=2;
% useNewton=true;
% deltaT=1.1;
% y0=0.5;
% aDrift=@(y,t)((1/3)*y^(1/3));
% BDiff=@(y,t)(y^(2/3));
% papy=@(y)1/(9*y^(2/3));
% papt=0;
% p2apypy=@(y)-(2/(27*y^(5/3)));
% pBpy=@(y)2/(3*y^(1/3));
% pBpt=0;
% p2Bpypy=@(y)-(2/(9*y^(4/3)));
% aFun=@(y,t)dealRobust(aDrift(y,t),papy(y),papt,p2apypy(y));
% BFun=@(y,t)dealRobust(BDiff(y,t),pBpy(y),pBpt,p2Bpypy(y));
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
%         y=implicitWeakTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExp);
%         BCur=BDiff(yExp);
%         pBpyCur=pBpy(yExp);
%         p2BpypyCur=p2Bpypy(yExp);
%         papyCur=papy(yExp);
%         p2apypyCur=p2apypy(yExp);
%         yExp=weakStochTaylorStep(yExp,aCur,BCur,deltaT/numSteps,algorithm,simplified,pBpyCur,p2BpypyCur,papyCur,p2apypyCur,papt,pBpt);
%         t=t+deltaT/numSteps;
%     end
%     valsRKImp(curMC)=y;
%     valsRK(curMC)=yExp;
% end
% muRKImp=mean(valsRKImp);
% muRK=mean(valsRK);
% abs((muRKImp-muCub)./muCub)%Relative mean error, implicit.
% abs((muRK-muCub)./muCub)%Relative mean error, explicit.
%The implicit error will be about 0.004 and the explicit error will be
%about 0.010. Thus, the algorithm improves the Euler-Maruyama method.
%
%EXAMPLE 3:
%Here, we compare the implicit Weak Taylor method to the explicit method on
%a linear model. This is done with the same noise driving both processes.
% algorithm=0;
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
% m=size(D,2);
% [F,Q]=linDynMod2Disc(deltaT,A,D);
% mu=F*y0;
% P=Q;
% aDrift=@(x,t)(A*x);
% BDiff=@(x,t)(D);
% 
% aFun=@(y,t)dealRobust(aDrift(y,t),A,zeros(d,1),zeros(d,d,d));
% BFun=@(y,t)dealRobust(BDiff(y,t),zeros(d,m,d),zeros(d,m),zeros(d,m,d,d));
% 
% valsRKImp=zeros(d,numMC);
% valsRK=zeros(d,numMC);
% for curMC=1:numMC
%     yExpl=y0;
%     y=y0;
%     t=0;
%     for curStep=1:numSteps
%         s=rng();%Record the state prior to generating the random variables.
%         y=implicitWeakTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExpl);
%         BCur=BDiff(yExpl);
%         simplified=2;
%         yExpl=weakStochTaylorStep(yExpl,aCur,BCur,deltaT/numSteps,algorithm,simplified,[],[],A);
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
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<12||isempty(AbsTol))
    AbsTol=0;
end

if(nargin<11||isempty(RelTol))
    RelTol=0;
end

if(nargin<10||isempty(thetaVals))
    theta1=1/2;
    theta2=0;
else
    if(length(thetaVals)==2)
        theta1=thetaVals(1);
        theta2=thetaVals(2);
    else
        theta1=thetaVals;
        theta2=0;
    end
end

if(nargin<9||isempty(maxIter))
    maxIter=2; 
end

if(nargin<8||isempty(useNewton))
    useNewton=false; 
end

aVal=[];
papy=[];
papt=[];
p2apypy=[];
BVal=[];
pBpy=[];
pBpt=[];
p2Bpypy=[];

if(nargin>6&&~isempty(givenVals))
    if(isfield(givenVals,'aVal'))
        aVal=givenVals.aVal;
    end
    if(isfield(givenVals,'papy'))
        papy=givenVals.papy;
    end
    if(isfield(givenVals,'papt'))
        papt=givenVals.papt;
    end
    if(isfield(givenVals,'p2apypy'))
        p2apypy=givenVals.p2apypy;
    end

    if(isfield(givenVals,'BVal'))
        BVal=givenVals.BVal;
    end
    if(isfield(givenVals,'pBpy'))
        pBpy=givenVals.pBpy;
    end
    if(isfield(givenVals,'pBpt'))
        pBpt=givenVals.pBpt;
    end
    if(isfield(givenVals,'p2Bpypy'))
        p2Bpypy=givenVals.p2Bpypy;
    end
end

if(nargin<6||isempty(algorithm))
    algorithm=0;
end

if(algorithm==0)
    %Weak Euler prediction (Equation 1.2 in Chapter 14.1) and
    %Trapezoid method correction by default, but allowing the family of
    %implicit schemes.
    if(isempty(aVal))
        aVal=a(y,curT); 
    end

    if(isempty(BVal))
       if(theta2==0)
           BVal=B(y,curT);
       else
           [BVal,pBpy]=B(y,curT);
       end
    end
    d=size(BVal,1);
    m=size(BVal,2);

    sqrtDeltaT=sqrt(deltaT);
    
    deltaW=rand(m,1);
    sel=(deltaW<=0.5);
    deltaW(sel)=-sqrtDeltaT;
    deltaW(~sel)=sqrtDeltaT;

    %Prediction.
    yNext=y+aVal*deltaT+BVal*deltaW;

    tNext=curT+deltaT;
    if(theta2==0)
        fixedTerm=(1-theta1).*aVal*deltaT+BVal*deltaW+y;

        if(useNewton)
            [yNext,didConverge]=implicitNewtonIter(y,a,theta1*deltaT,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        else
            [yNext,didConverge]=implicitFixedPointIter(y,a,theta1*deltaT,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        end
    else
        aBar=aVal;
        for j=1:m
            for k=1:d
                aBar=aBar-theta2.*BVal(k,j)*pBpy(:,j,k);
            end
        end

        fixedTerm=y+(1-theta1).*aBar*deltaT+(1-theta2).*BVal*deltaW;
        f=@(y,t)a1Mod(y,t,a,B,theta1,theta2,deltaT,deltaW);

        if(useNewton)
            error('Iteration with Newton''s method is not supported for this algorithm unless theta2=0.')
        else
            [yNext,didConverge]=implicitFixedPointIter(y,f,deltaT,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        end
    end
elseif(algorithm==1||algorithm==2)
    if(algorithm==1)
        if(isempty(aVal))
            [aVal,papy,~,p2apypy]=a(y,curT);
        end
        
        if(isempty(BVal))
            [BVal,pBpy,~,p2Bpypy]=B(y,curT);
        end
        [yNext,fixedTerm]=weakStochTaylorStep(y,aVal,BVal,deltaT,algorithm,2,pBpy,p2Bpypy,papy,p2apypy);
    else
        if(isempty(aVal))
            [aVal,papy,papt,p2apypy]=a(y,curT);
        end
        
        if(isempty(BVal))
            [BVal,pBpy,pBpt,p2Bpypy]=B(y,curT);
        end
        [yNext,fixedTerm]=weakStochTaylorStep(y,aVal,BVal,deltaT,algorithm,2,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt);
    end
    tNext=curT+deltaT;
    fixedTerm=fixedTerm+(1/2)*aVal*deltaT;

    if(useNewton)
        [yNext,didConverge]=implicitNewtonIter(y,a,deltaT/2,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
    else
        [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT/2,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
    end
else
    error('Unknown algorithm specified.')
end
end

function val=a1Mod(y,t,a,B,theta1,theta2,deltaT,deltaW)
%A1MOD This implements the implicit terms in Equation 5.4 of Chapter 15.5
%      of [1].
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

[BVal,pBpy]=B(y,t);
aVal=a(y,t);

d=size(BVal,1);
m=size(BVal,2);

aBar=aVal;
for j=1:m
    for k=1:d
        aBar=aBar-theta2*BVal(k,j)*pBpy(:,j,k);
    end
end

val=theta1*aBar*deltaT+theta2*BVal*deltaW;

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
