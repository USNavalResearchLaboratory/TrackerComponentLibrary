function [yNext,didConverge]=semiImplicitStrongTaylorStep(y,curT,a,B,deltaT,algorithm,givenVals,useNewton,maxIter,thetaVals,RelTol,AbsTol,p,ItoAlg)
%%SEMIIMPLICITSTRONGTAYLORSTEP Perform a single step of a semi-implicit
%           strong Itô-Taylor expansion to integrate a d-dimensional
%           stochastic differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process and t is time. Semi-implicit techniques are slower, but
%           typically have better stability compared to explicit
%           algorithms. As noted in Chapter 9.8 of [1], fully implicit
%           strong Itô-Taylor expansion expansions are unstable, because
%           particular values of the noise process can make them diverge.
%           Thus, one typically uses semi-implicit methods, as discussed in
%           Chapter 12 of [1]. Strong methods converge to an optimal path
%           as the stepsize decreases, not considering finite precision
%           errors.
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
%          the first output of a and all but 0.1 use the first output of B.
%          Additionally, all methods that support a Newton iteration use
%          papy. The methods with fractional numbers are typically not that
%          great. Possible values of algorithm are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            semi-implicit Euler-Maruyama method (order 0.5) of Equation
%            2.2 of Chapter 12.2 of [1]. This is the semi-implicit version
%            of Equation 2.1 in Chapter 10.2 of [1].
%          0.1 Use split-step backward Euler-Maruyama algorithm (order 0.5)
%            of [2].
%          0.2 Perform the error corrected Euler-Maruyama algorithm (order
%            0.5) of [3]. This is not a true implicit method, because there
%            is no iteration and the useNewton input is ignored. This
%            algorithm does not use theta1 and theta2.
%          1 Use the semi-implicit Milstein scheme (order 1.0) for scalar
%            problems; Equation 2.9 in Chapter 12.2 of [1], but generalized
%            to a multivariate state with scalar noise. This is the semi-
%            implicit version of Equation 3.2 in Chapter 10.3. The equation
%            has been generalized to support theta1 values other than one,
%            analogously to what was done in Equation 2.13 in Chapter 12.2.
%            of [1]. This required pBpy and that m=1.
%          1.1 The improved semi-implicit Milstein scheme (order 1.0) for
%            scalar problems; Equation 1.2 in [4]. This requires papy and
%            pBpy.
%          2 Use the semi-implicit Milstein scheme for general multivariate
%            noise (order 1.0); Equation 2.13 of Chapter 12.2 of [1]. This
%            is the semi-implicit form of Equation 3.3 of Chapter 10.3 of
%            [1]. This requires pBpy.
%          2.1 The semi-implicit Milstein scheme (order 1.0) for general
%            multivariate noise; Equation 2.17 in [4]. This requires papy
%            and pBpy.
%          3 Use the semi-implicit Milstein scheme for diagonal noise. This
%            relates to explicit method Equation 3.12 of Chapter 10.3 in
%            the same manner that Equation 2.13 of Chapter 12.2 of [1] is
%            the semi-implicit form of Equation 3.3 of Chapter 10.3 of [1].
%            This requires pBpy and that d=m.
%          4 Use the semi-implicit Milstein scheme for commutative noise;
%            Equation 2.11 in Chapter 12.2. of [1]. This is the semi-
%            implicit form of Equation 3.16 of Chapter 10.3 of [1]
%            (modified Itô rather than Stratonovich). The equation has
%            been generalized to support theta1 values other than one,
%            analogously to Equation 2.13 in Chapter 12.2. of [1]. This
%            requires pBpy and that d=m.
%          5 Use the semi-implicit strong order 1.5 Taylor method for
%            autonomous scalar problems; Equation 2.15 in Chapter 12.2 of
%            [1]. This is the semi-implicit form of Equation 4.1 of Chapter
%            10.4 of [1]. This requires pBpy,p2Bpypy,papy, p2apypy, that
%            d=m=1, and that a and B don't depend on t. This algorithm does
%            not use theta1; rather the expression used assumes
%            theta1=theta2=1/2.
%          6 Use the semi-implicit strong order 1.5 Taylor method for
%            additive noise; Equation 2.17 of Chapter 12.2 of [1]. This is
%            the semi-implicit form of Equation 4.10 of Chapter 10.4 of
%            [1]. This requires papy, p2apypy, papt, and pBpt and that B
%            not depend on x. Except when theta1=theta2, this algorithm
%            only supports useNewton=false.
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
% thetaVals All of the algorithms except 0.2, 1.1 and 2.1 have the option
%          for a parameter, 0<=theta1<=1, that influences how implicit the
%          step is in with respect to the drift term. Algorithm 6 has an
%          additional value 0<=theta2<=1 that effects the level of
%          implicitness of the derivatives. thetaVals=[theta1,theta2] or
%          just thetaVals=theta1 if theta2 is not needed or one wishes for
%          it to just be 1. A value of 1 is full implicitness. theta1 and
%          theta2 can be dX1 vectors if
%          different weightings for different dimensions are desired, or
%          they can just be scalars. The default if omitted or an empty
%          matrix is passed is [1/2,1/2].
% RelTol, AbsTol The relative and absolute tolerances on the iterations
%          before declaring convergence. If these are set to 0 (the default
%          if omitted or empty matrices are passed), then the algorithm
%          will just iterate for the maximum number of iterations. The
%          tolerances apply to each element of y. Convergence is declared
%          if all(diff<=AbsTol)||all(diff<=RelTol*abs(yNext)), where diff
%          is the change in y between subsequent steps.
% p, ItoAlg These values are only used if  algorithm=2. These values
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
% algorithm=0;
% useNewton=true;
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
%         y=semiImplicitStrongTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExpl);
%         BCur=BDiff(yExpl);
%         yExpl=strongStochTaylorStep(yExpl,aCur,BCur,deltaT/numSteps,algorithm,[],[],A);
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
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its solution
%are from Equation 4.40. We compare the performance of the explicit
%solution (taking the expected value using cubature integration) and the
%implicit solution.
% rng(1)%Make exact run repeatable.
% numMC=1e3;
% numSteps=1;
% algorithm=5;
% 
% useNewton=false;
% deltaT=1/4;
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
%         y=semiImplicitStrongTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExp);
%         BCur=BDiff(yExp);
%         pBpyCur=pBpy(yExp);
%         p2BpypyCur=p2Bpypy(yExp);
%         papyCur=papy(yExp);
%         p2apypyCur=p2apypy(yExp);
%         yExp=strongStochTaylorStep(yExp,aCur,BCur,deltaT/numSteps,algorithm,pBpyCur,p2BpypyCur,papyCur,p2apypyCur,papt,pBpt);
%         t=t+deltaT/numSteps;
%     end
%     valsRKImp(curMC)=y;
%     valsRK(curMC)=yExp;
% end
% muRKImp=mean(valsRKImp);
% muRK=mean(valsRK);
% abs((muRKImp-muCub)./muCub)%Relative mean error, implicit.
% abs((muRK-muCub)./muCub)%Relative mean error, explicit.
%The implicit error will be about 0.0067 and the explicit error will be
%about 0.0176. Thus, the algorithm improves the order 1.5 Taylor method.
%
%EXAMPLE 3:
%This is an example of what is considered a "stiff" problem in Section 12.2
%of [1] with non-additive noise. The problem has an explicit solution in
%terms of the Wiener process W. We use quadrature integration to obtain the
%mean of that solution, which is used as the true moments for comparison.
%We compare this implicit method with an explicit step. This is a 2X1
%process with scalar noise.
% rng(0)%Make the exact run repeatable.
% algorithm=1;
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
%         y=semiImplicitStrongTaylorStep(y,t,aFun,BFun,deltaT/numSteps,algorithm,[],useNewton);
% 
%         rng(s)%Drive the explicit step with the same random process as the
%               %implicit one.
%         aCur=aDrift(yExp,t);
%         BCur=BDiff(yExp,t);
%         yExp=strongStochTaylorStep(yExp,aCur,BCur,deltaT/numSteps,algorithm,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt);
%         t=t+deltaT/numSteps;
%     end
%     valsRK(:,curMC)=y;
%     valsRKExplicit(:,curMC)=yExp;
% end
% muRK=mean(valsRK,2);
% muRKExp=mean(valsRKExplicit,2);
% norm((muRK-muCub)./norm(muCub))%Relative mean error, implicit.
% norm((muRKExp-muCub)./norm(muCub))%Relative mean error, explicit.
%The implicit solution is typically more accurate than the explicit. The
%same relation between explicit and implicit solutions should hold when
%using algorithm 0 or 2. Note that we are comparing the norms, not the
%individual components.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] D. J. Higham, X. Mao, and A. M. Stuart, "Strong convergence of Euler-
%    type methods for nonlinear stochastic differential equations," SIAM
%    Journal on Numerical Analysis, vol. 40, no. 3, pp. 1041-1063, 2002.
%[3] Z. Yin and S. Gan, "An error-corrected Euler-Maruyama method for stiff
%    stochastic differential equations," Applied Mathematics and
%    Computation, vol. 25, pp. 630-641, 1 Apr. 2015.
%[4] Z. Yin and S. Gan, "An improved Milstein method for stiff stochastic
%    differential equations," Advances in Difference Equations, vol. 369, 1
%    Dec. 2015.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

d=size(y,1);

if(nargin<14||isempty(ItoAlg))
    ItoAlg=0;
end

if(nargin<13||isempty(p))
    p=15;
end

if(nargin<12||isempty(AbsTol))
    AbsTol=0;
end

if(nargin<11||isempty(RelTol))
    RelTol=0; 
end

if(nargin<10||isempty(thetaVals))
    theta1=1/2;
    theta2=1/2;
else
    if(length(thetaVals)==2)
        theta1=thetaVals(1);
        theta2=thetaVals(2);
    else
        theta1=thetaVals;
        theta2=1/2;
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

if(algorithm>=0&&algorithm<=4)
    %If it is an order 0.5 or an order 1.0 method.
    if(algorithm<1&&algorithm~=fix(algorithm))
        switch(algorithm)
            case 0.1%If we are to perform the split-step backward Euler-
                  %Maruyama algorithm of [2].
                if(isempty(aVal))
                    if(useNewton)
                        [aVal,papy]=a(y,curT);
                        initVals.fCur=aVal;
                        initVals.df=papy;
                    else
                        aVal=a(y,curT);
                        initVals.fCur=aVal;
                    end
                else
                    initVals=[];
                    initVals.fCur=aVal;
                    initVals.df=papy;
                end

                [yNext,didConverge]=implicitEulerStep(y,curT,a,deltaT,initVals,useNewton,maxIter,theta1,RelTol,AbsTol);
                BVal=B(yNext,curT+deltaT);
                m=size(BVal,2);
                
                deltaW=sqrt(deltaT)*randn(m,1);

                yNext=yNext+BVal*deltaW;
            case 0.2
                if(isempty(aVal))
                    aVal=a(y,curT);
                end

                if(isempty(BVal))
                    BVal=B(y,curT); 
                end
                
                %The explicit step of the Euler-Maruyama method.
                yNext=strongStochTaylorStep(y,aVal,BVal,deltaT,0);
                
                %If we are to perform the error corrected Euler-Maruyama
                %algorithm of [3].
                [aNext,daNext]=a(yNext,curT+deltaT);

                I=eye(d,d);
                yNext=yNext+deltaT*pinv(I-deltaT*daNext)*(aNext-aVal);
                didConverge=[];
            otherwise
                error('Unknown variant of Euler''s algorithm selected.')
        end
    elseif(algorithm==1.1)%Improved Milstein scheme of [1] for scalar
                          %problems; Equation 2.15 in [4].
        if(isempty(aVal))
            aVal=a(y,curT);
        end
        if(isempty(BVal))
            [BVal,pBpy]=B(y,curT);
        end
        m=size(BVal,2);
        
        if(m~=1||d~=1)
            error('This algorithm is only for scalar problems.')
        end
        
        deltaW=sqrt(deltaT)*randn(m,1);
        
        %Equation 2.15
        yBar=y+aVal*deltaT+BVal*deltaW+(1/2)*pBpy*BVal*(deltaW^2-deltaT);
        [aBar,paBarpt]=a(yBar,curT+deltaT);
        
        %If 1-deltaT*paNextpt=0, then pinv will return 0. 
        yNext=yBar+pinv(1-deltaT*paBarpt)*(aBar-aVal)*deltaT;
        
        didConverge=[];
    elseif(algorithm==2.1)%Improved Milstein scheme of [1] for general
                          %problems.
        if(isempty(aVal))
            aVal=a(y,curT);
        end
        if(isempty(BVal))
            [BVal,pBpy]=B(y,curT);
        end
        m=size(BVal,2);
        %Stochastic integrals are needed.
        [I1Idx,I2Idx]=sampleItoIntegrals(m,deltaT,p,ItoAlg);
        
        deltaW=I1Idx.Ij;
        Ij1j2=I2Idx.Ij1j2;
        LjVal=LjOperator(BVal,pBpy);
        
        yBar=y+aVal*deltaT+BVal*deltaW+sum(sum(bsxfun(@times,LjVal,reshape(Ij1j2,[1,m,m])),2),3);
        [aBar,paBarpt]=a(yBar,curT+deltaT);
        
        yNext=yBar+pinv(eye(d,d)-paBarpt*deltaT)*(aBar-aVal)*deltaT;

        didConverge=[];               
    else%Implicit Euler or Milstein algorithm variant.
        if(isempty(aVal))
            aVal=a(y,curT);
        end

        if(isempty(BVal))
            if(algorithm>0)
                [BVal,pBpy]=B(y,curT);
            else
                BVal=B(y,curT);
                pBpy=[];
            end
        end
        
        [yNext,endTerm]=strongStochTaylorStep(y,aVal,BVal,deltaT,algorithm,pBpy,[],[],[],[],[],p,ItoAlg);

        fixedTerm=(1-theta1).*aVal*deltaT+endTerm;
        tNext=curT+deltaT;
        if(useNewton)
            [yNext,didConverge]=implicitNewtonIter(y,a,deltaT*theta1,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        else
            [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT*theta1,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        end
    end
elseif(algorithm<=6)%Implicit order 1.5 algorithm variant.
    if(algorithm==5)
        if(isempty(aVal))
            [aVal,papy,~,p2apypy]=a(y,curT);
        end
        
        if(isempty(BVal))
            [BVal,pBpy,~,p2Bpypy]=B(y,curT);
        end
        
        %The strong order 1.5 Taylor method for autonomous scalar
        %problems; Equation 2.15 in Chapter 12.2 of [1]. This is
        %an implicit form of Equation 4.1 of Chapter 10.4 of [1].
        [yNext,endTerm,randTerms]=strongStochTaylorStep(y,aVal,BVal,deltaT,algorithm,pBpy,p2Bpypy,papy,p2apypy,[],[],p,ItoAlg);
        deltaW=randTerms.deltaW;
        deltaZ=randTerms.deltaZ;

        fixedTerm=(1/2)*aVal*deltaT+papy*BVal*(deltaZ-(1/2)*deltaW*deltaT)+endTerm;
        tNext=curT+deltaT;
        if(useNewton)
            [yNext,didConverge]=implicitNewtonIter(y,a,deltaT*theta1,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        else
            [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT*theta1,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
        end
    else%Algorithm==6.
        if(isempty(aVal))
            [aVal,papy,papt,p2apypy]=a(y,curT);
        end
        
        if(isempty(BVal))
            [BVal,~,pBpt]=B(y,curT);
        end

        [yNext,endTerm,randTerms]=strongStochTaylorStep(y,aVal,BVal,deltaT,algorithm,[],[],papy,p2apypy,papt,pBpt,p,ItoAlg);
        deltaW=randTerms.deltaW;
        deltaZ=randTerms.deltaZ;

        t1=endTerm.val;
        L0a=endTerm.L0a;
        Lja=endTerm.Lja;

        fixedTerm=(1-theta1).*aVal*deltaT+(1/2-theta1).*(1-theta2).*L0a*deltaT^2+Lja*(deltaZ-theta1*deltaW*deltaT)+t1;
        tNext=curT+deltaT;
        if(all(theta1==1/2))
            if(useNewton)
                [yNext,didConverge]=implicitNewtonIter(y,a,deltaT*(1/2),yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
            else
                [yNext,didConverge]=implicitFixedPointIter(y,a,deltaT*(1/2),yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
            end
        else
            aFun=@(x,t)a15Mod(x,t,a,BVal,deltaT,theta1,theta2);

            if(useNewton)
                error('Iteration with Newton''s method is not supported for this algorithm unless theta1=1/2.')
            else
                [yNext,didConverge]=implicitFixedPointIter(y,aFun,deltaT,yNext,tNext,fixedTerm,maxIter,RelTol,AbsTol);
            end
        end
    end
else
    error('Unknown algorithm specified.')
end
end

function val=a15Mod(y,t,a,BVal,deltaT,theta1,theta2)
%%A15MOD This implements the implicit terms in Equation 2.17 of Chapter
%        12.2 of [1]. This involves a and the L0 operator applied to a.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

[aNext,papy,papt,p2apypy]=a(y,t);
L0a=L0Operator(aNext,BVal,papt,papy,p2apypy);

val=theta1.*aNext*deltaT+((1/2)-theta1).*theta2*L0a;
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
