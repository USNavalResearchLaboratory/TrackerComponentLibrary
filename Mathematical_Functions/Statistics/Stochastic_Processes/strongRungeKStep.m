function [y,fixedTerm]=strongRungeKStep(y,t,a,B,deltaT,algorithm,aCur,BCur,p,ItoAlg)
%%STRONGRUNGEKSTEP Perform a single step of an explicit strong stochastic
%           Runge-Kutta method under It� calculus. This integrates a d-
%           dimensional stochastic differential equation of the form
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process. Strong methods converge to an optimal path as the
%           stepsize decreases, not considering finite precision errors.
%           Runge-Kutta methods use finite differencing instead of
%           explicit derivatives.
%
%INPUTS: y The dX1 initial value of the random process.
%        t The scalar initial time of the random process. If an empty
%          matrix is passed, t=0 is used.
%        a A function handle to the drift function. This is called as
%          a(y,t) and returns a dX1 vector.
%        B A function handle to the diffusion matrix function. This is
%          called as B(y,t) and returns a dXm matrix.
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. Possible values
%          are:
%          -1 (The default if omitted or an empty matrix is passed) Use the
%             Euler-Maruyama method from Equation 2.4 in Chapter 10.2 of
%             [1]. This is an order 0.5 method. The other algorithms will
%             typically outperform the Euler-Maruyama method.
%           0 Use the explicit order 1.0 method for scalar noise from
%             Equation 1.5 in Chapter 11.1 of [1]. This requires that m=1.
%           1 Use the explicit order 1.0 method for general noise from
%             Equation 1.7 in Chapter 11.1 of [1].
%           2 Use the explicit order 1.0 method for diagonal noise from
%             Equation 1.9 in Chapter 11.1 of [1]. This requires that d=m.
%           3 Use the explicit order 1.5 method for autonomous scalar
%             problems from Equation 2.1 in Chapter 11.2 of [1]. This
%             requires that d=m=1 and that a and B are not functions of t.
%           4 Use the explicit order 1.5 method for non-autonomous additive
%             noise from Equation 2.19 in Chapter 11.2 of [1]. This
%             requires that B is not a function of x.
%           5 Use the explicit order 1.5 method for autonomous additive
%             noise from Equation 2.7 in Chapter 11.2 of [1]. This requires
%             that a and B are not functions of t and that B is not a
%             function of x.
% aCur, BCur Often one might already have the values a(x,t) and B(x,t). If
%          so, then they should be provided as the dX1 and dXm aCur and
%          BCur to avoid recalculation. If unavailable, these values can be
%          omitted or empty matrices passed.
% p, ItoAlg These values are only used if algorithm=1. These values
%          correspond to the inputs p and algorithm of the function
%          sampleItoIntegrals, which is called by this function for
%          necessary terms. Defaults if omitted or an empty matrix is
%          passed are 15 and 0.
%
%OUTPUTS: y The estimated value of the process after taking a step of
%           deltaT. This is a random value.
% fixedTerm This output is used by the semiImplicitStrongRungeKStep
%           function. It is the fixed term in the iterations when using the
%           implicit version of the algorithm.
%
%EXAMPLE 1:
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its implicit
%solution are from Equation 4.40.
% numMC=1e4;
% numSteps=3;
% algorithm=3;
% deltaT=1.1;
% x0=0.5;
% aDrift1=@(x,t)((1/3)*x^(1/3));
% BDiff1=@(x,t)(x^(2/3));
% explSim1=@(W)(x0^(1/3)+(1/3)*W)^3;
% 
% valsSim=zeros(1,numMC);
% valsRK=zeros(1,numMC);
% for curMC=1:numMC
%     W=sqrt(deltaT)*randn(1);
%     valsSim(curMC)=explSim1(W);
% 
%     y=x0;
%     for curStep=1:numSteps
%         y=strongRungeKStep(y,[],aDrift1,BDiff1,deltaT/numSteps,algorithm);
%     end
%     valsRK(curMC)=y;
% end
% [muSim,PSim]=calcMixtureMoments(valsSim);
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-muSim)./muSim)%Relative mean error.
% abs((PRK-PSim)./PSim)%Relative variance error.
%The error in the estimate of the mean and the variance should be less than
%10%. However, if one were to switch to algorithm=-1, the Euler-Maruyama
%method, then the error will often be over 10%.
%
%EXAMPLE 2:
%Here, we consider a multivariate linear dynamic model with multivariate
%noise. The mean and covariance matrix of the sample paths are compared to
%the exact solution.
% algorithm=5;
% numMC=1e4;
% deltaT=1/3;
% numSteps=3;
% x0=[1/4;-12];
% A=[1.1,0.1;
%    -0.2,2.2];
% D=[1.5,0.4;
%     0.0,1];
% [F,Q]=linDynMod2Disc(deltaT,A,D);
% mu=F*x0;
% P=Q;
% aDrift=@(x,t)(A*x);
% BDiff=@(x,t)(D);
% 
% valsRK=zeros(2,numMC);
% for curMC=1:numMC
%     y=x0;
%     for curStep=1:numSteps
%         y=strongRungeKStep(y,[],aDrift,BDiff,deltaT/numSteps,algorithm);
%     end
%     valsRK(:,curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-mu)./mu)%Relative mean error.
% abs((PRK-P)./P)%Relative variance error.
%The errors in the estimate of the mean and the variance elements should be
%less than 10%. However, if one were to switch to algorithm=-1, the
%Euler-Maruyama method, then the errors will often be over 10%.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(t))
    t=0; 
end

if(nargin<10||isempty(ItoAlg))
    ItoAlg=0;
end

if(nargin<9||isempty(p))
    p=15; 
end

if(nargin<8||isempty(aCur))
    aCur=a(y,t);
end

if(nargin<7||isempty(BCur))
    BCur=B(y,t);
end

if(nargin<6||isempty(algorithm))
    algorithm=-1;
end

d=size(BCur,1);
m=size(BCur,2);
    
switch(algorithm)
     case -1%The Euler-Maruyama method; Equation 2.4 of [1].
        sqrtDeltaT=sqrt(deltaT);
        deltaW=sqrtDeltaT*randn(m,1);
        
        fixedTerm=y+BCur*deltaW;
        y=aCur*deltaT+fixedTerm;
    case 0%The explicit order 1.0 method for scalar noise; Equation 1.5 in
          %Chapter 11.1 of [1].
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
        sqrtDeltaT=sqrt(deltaT);
        deltaW=sqrtDeltaT*randn(1);
        
        yBar=y+aCur*deltaT+BCur*sqrtDeltaT;
        BBar=B(yBar,t);
        
        fixedTerm=y+BCur*deltaW+(1/(2*sqrtDeltaT))*(BBar-BCur)*(deltaW.^2-deltaT);
        
        y=aCur*deltaT+fixedTerm;
    case 1%The explicit order 1.0 method for general noise; Equation 1.7 in
          %Chapter 11.1 of [1].

        sqrtDeltaT=sqrt(deltaT);
        %Stochastic integrals are needed.
        [I1Idx,I2Idx]=sampleItoIntegrals(m,deltaT,p,ItoAlg);
        deltaW=I1Idx.Ij;
        Ij1j2=I2Idx.Ij1j2;
        
        sumVal=zeros(d,1);
        for j1=1:m
            yBar=y+aCur*deltaT+BCur(:,j1)*sqrtDeltaT;
            BDiff=B(yBar,t)-BCur;
            
            for j2=1:m
                sumVal=sumVal+BDiff(:,j2)*Ij1j2(j1,j2);
            end
        end
        
        fixedTerm=y+BCur*deltaW+(1/sqrtDeltaT)*sumVal;
        
        y=aCur*deltaT+fixedTerm;
    case 2%The explicit order 1.0 method for diagonal noise; Equation 1.9
          %in Chapter 11.1 of [1].
        if(m~=d)
           error('This algorithm is only valid for d=m.');
        end
          
        sqrtDeltaT=sqrt(deltaT);
        deltaW=sqrtDeltaT*randn(m,1);
        
        BCur=diag(BCur);
        yBar=y+aCur*deltaT+BCur*sqrtDeltaT;
        BBar=diag(B(yBar,t));

        fixedTerm=y+BCur.*deltaW+(1/(2*sqrtDeltaT))*(BBar-BCur).*(deltaW.^2-deltaT);
        
        y=aCur*deltaT+fixedTerm;
    case 3%The explicit order 1.5 method for autonomous scalar problems;
          %Equation 2.1 in Chapter 11.2 of [1].
        sqrtDeltaT=sqrt(deltaT);

        if(m~=1||d~=1)
            error('This algorithm is only for scalar problems.')
        end

        U1=randn(1);
        U2=randn(1);
        deltaW=U1*sqrtDeltaT;
        deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);

        yBarPlus=y+aCur*deltaT+BCur*sqrtDeltaT;
        yBarMinus=y+aCur*deltaT-BCur*sqrtDeltaT;
        aBarPlus=a(yBarPlus,t);
        aBarMinus=a(yBarMinus,t);
        BBarPlus=B(yBarPlus,t);
        BBarMinus=B(yBarMinus,t);

        PhiBarPlus=yBarPlus+BBarPlus*sqrtDeltaT;
        PhiBarMinus=yBarPlus-BBarPlus*sqrtDeltaT;
        BPhiPlus=B(PhiBarPlus,t);
        BPhiMinus=B(PhiBarMinus,t);
        
        fixedTerm=y+BCur*deltaW...
                   +(1/(4*sqrtDeltaT))*(BBarPlus-BBarMinus)*(deltaW.^2-deltaT)...
                   +(1/(2*deltaT))*(BBarPlus-2*BCur+BBarMinus)*(deltaW*deltaT-deltaZ)...
                   +(1/(4*deltaT))*(BPhiPlus-BPhiMinus-BBarPlus+BBarMinus)*((1/3)*deltaW.^2-deltaT).*deltaW;

        temp=(1/(2*sqrtDeltaT))*(aBarPlus-aBarMinus);
               
        y=temp*deltaZ...
          +(1/4)*(aBarPlus+2*aCur+aBarMinus)*deltaT+fixedTerm;
      
        if(nargout>1)
            fixedTerm=fixedTerm+temp*(deltaZ-(1/2)*deltaW*deltaT)+(1/2)*aCur*deltaT;
        end
    case 4%The explicit order 1.5 method for non-autonomous additive noise;
          %Equation 2.19 in Chapter 11.2 of [1].
        sqrtDeltaT=sqrt(deltaT);

        U1=randn(m,1);
        U2=randn(m,1);
        deltaW=U1*sqrtDeltaT;
        deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);

        tNext=t+deltaT;
        deltaB=B(y,tNext)-BCur;

        sum1=zeros(d,1);
        sum2=zeros(d,1);
        sum3=zeros(d,1);
        if(nargout>1)
            %If the fixed term for the implicit algorithm is desired.
            sum2Imp=zeros(d,1);
            for j=1:m
                yBarPlus=y+(1/m)*aCur*deltaT+BCur(:,j)*sqrtDeltaT;
                yBarMinus=y+(1/m)*aCur*deltaT-BCur(:,j)*sqrtDeltaT;
                aBarPlus=a(yBarPlus,tNext);
                aBarMinus=a(yBarMinus,tNext);

                sum1=sum1+(aBarPlus-2*aCur+aBarMinus)*deltaT;

                diffA=aBarPlus-aBarMinus;
                sum2=sum2+diffA*deltaZ(j);
                sum2Imp=sum2Imp+diffA*(deltaZ(j)-(1/2)*deltaW(j)*deltaT);
                sum3=sum3+deltaB(:,j)*(deltaW(j)*deltaT-deltaZ(j));
            end
        else
            for j=1:m
                yBarPlus=y+(1/m)*aCur*deltaT+BCur(:,j)*sqrtDeltaT;
                yBarMinus=y+(1/m)*aCur*deltaT-BCur(:,j)*sqrtDeltaT;
                aBarPlus=a(yBarPlus,tNext);
                aBarMinus=a(yBarMinus,tNext);

                sum1=sum1+(aBarPlus-2*aCur+aBarMinus)*deltaT;
                sum2=sum2+(aBarPlus-aBarMinus)*deltaZ(j);
                sum3=sum3+deltaB(:,j)*(deltaW(j)*deltaT-deltaZ(j));
            end
        end

        fixedTerm=y+BCur*deltaW+(1/deltaT)*sum3;
        
        y=aCur*deltaT+(1/4)*sum1+(1/(2*sqrtDeltaT))*sum2+fixedTerm;
        
        if(nargout>1)
            fixedTerm=fixedTerm+(1/2)*aCur*deltaT+sum2Imp/(2*sqrtDeltaT);
        end
    case 5%The explicit order 1.5 method for autonomous additive noise;
          %Equation 2.7 in Chapter 11.2 of [1].
        sqrtDeltaT=sqrt(deltaT);

        U1=randn(m,1);
        U2=randn(m,1);
        deltaW=U1*sqrtDeltaT;
        deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);

        sum1=zeros(d,1);
        sum2=zeros(d,1);
        if(nargout>1)
            %If the fixed term for the implicit algorithm is desired.
            sum1Imp=zeros(d,1);
            for j=1:m
                yBarPlus=y+(1/m)*aCur*deltaT+BCur(:,j)*sqrtDeltaT;
                yBarMinus=y+(1/m)*aCur*deltaT-BCur(:,j)*sqrtDeltaT;
                aBarPlus=a(yBarPlus,t);
                aBarMinus=a(yBarMinus,t);

                diffA=aBarPlus-aBarMinus;
                sum1=sum1+diffA*deltaZ(j);
                sum1Imp=sum1Imp+diffA*(deltaZ(j)-(1/2)*deltaW(j)*deltaT);

                sum2=sum2+aBarPlus-((2*(m-2))/m)*aCur+aBarMinus;
            end
        else
            for j=1:m
                yBarPlus=y+(1/m)*aCur*deltaT+BCur(:,j)*sqrtDeltaT;
                yBarMinus=y+(1/m)*aCur*deltaT-BCur(:,j)*sqrtDeltaT;
                aBarPlus=a(yBarPlus,t);
                aBarMinus=a(yBarMinus,t);

                sum1=sum1+(aBarPlus-aBarMinus)*deltaZ(j);
                sum2=sum2+aBarPlus-((2*(m-2))/m)*aCur+aBarMinus;
            end
        end

        fixedTerm=y+BCur*deltaW;
        y=(1/(2*sqrtDeltaT))*sum1+(deltaT/4)*sum2+fixedTerm;
        
        if(nargout>1)
            fixedTerm=fixedTerm+(1/2)*aCur*deltaT+sum1Imp/(2*sqrtDeltaT);
        end
    otherwise
        error('Unknown algorithm specified.')
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
