function [y,fixedTerm]=weakRungeKStep(y,t,a,B,deltaT,algorithm,aCur,BCur,useGaussian)
%%WEAKRUNGEKSTEP Perform a single step of an explicit weak stochastic
%           Runge-Kutta method under It� calculus. This integrates a d-
%           dimensional stochastic differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process. As the stepsize used decreases, weak methods converge
%           such that integrals with the random process are a measure are
%           correct (for example, to determine moments). However, they do
%           not converge to the optimal path, unlike strong methods.
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
%          0 Use the explicit order 2.0 scheme for scalar noise from
%            Equation 1.1 in Chapter 15.1 of [1]. This requires that m=1.
%          1 Use the autonomous explicit order 2.0 weak scheme from
%            Equation 1.3 in Chapter 15.1 of [1]. This requires that a and
%            B not depend on t.
%          2 Use the autonomous explicit order 2.0 weak scheme for additive
%            noise from Equation 1.4 in Chapter 15.1 of [1]. This requires
%            that a and B not depend on t and that B does not depend on x.
% aCur, BCur Often one might already have the values a(x,t) and B(x,t). If
%           so, then they should be provided as the dX1 and dXm aCur and
%           BCur to avoid recalculation. If unavailable, these values can
%           be omitted or empty matrices passed.
% useGaussian Algorithms 0 and 1 have a choice of how the random component
%           is generated. If useGaussian=true, then Gaussian random
%           variables will be used. Otherwise, simpler random variables
%           having the same desired moment properties will be used. the
%           default if omitted or an empty matrix is passed is true.
%
%OUTPUTS: y The estimated value of the process after taking a step of
%           deltaT. This is a random value.
% fixedTerm This output is used by the impicitWeakRungeKStep function. It
%           is the fixed term in the iterations when using the implicit
%           version of the algorithm.
%
%EXAMPLE 1:
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its solution
%are from Equation 4.41.
% numMC=1e4;
% numSteps=1;
% algorithm=0;
% deltaT=0.9;
% x0=0.5;
% aDrift=@(x,t)((1/3)*x^(1/3));
% BDiff=@(x,t)(x^(2/3));
% explSim=@(W)(x0^(1/3)+(1/3)*W)^3;
% 
% valsSim=zeros(1,numMC);
% valsRK=zeros(1,numMC);
% for curMC=1:numMC
%     W=sqrt(deltaT)*randn(1);
%     valsSim(curMC)=explSim(W);
% 
%     y=x0;
%     for curStep=1:numSteps
%         y=weakRungeKStep(y,[],aDrift,BDiff,deltaT/numSteps,algorithm);
%     end
%     valsRK(curMC)=y;
% end
% [muSim,PSim]=calcMixtureMoments(valsSim);
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-muSim)./muSim)%Relative mean error.
% abs((PRK-PSim)./PSim)%Relative variance error.
%The error in the estimate of the mean and the variance should be less than
%4%. One can verify that this is better than what one would get with
%the Euler-Maruyama method in strongRungeKStep. 
%
%EXAMPLE 2:
%Here, we consider a multivariate linear dynamic model with multivariate
%noise. The mean and covariance matrix of the sample paths are compared to
%the exact solution.
% algorithm=2;
% numMC=1e4;
% deltaT=1/3;
% numSteps=5;
% x0=[1/4;-12];
% A=[1.1,0.1;
%    -0.2,2.2];
% D=[1.5;
%     0.1];
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
%         y=weakRungeKStep(y,[],aDrift,BDiff,deltaT/numSteps,algorithm);
%     end
%     valsRK(:,curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-mu)./mu)%Relative mean error.
% abs((PRK-P)./P)%Relative variance error.
%The error in the estimate of the mean and the variance should be less than
%4%. One can verify that this is better than what one would get with
%the Euler-Maruyama method in strongRungeKStep. 
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9||isempty(useGaussian))
    useGaussian=true;
end

if(nargin<8||isempty(aCur))
    aCur=a(y,t);
end

if(nargin<7||isempty(BCur))
    BCur=B(y,t);
end

if(isempty(t))
    t=0; 
end

d=size(BCur,1);
m=size(BCur,2);

switch(algorithm)
    case 0%The explicit order 2.0 weak scheme for scalar noise; Equation
          %1.1 in Chapter 15.1 of [1].
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
        
        sqrtDeltaT=sqrt(deltaT);
        if(useGaussian)
            deltaW=sqrtDeltaT*randn(1);
        else
            U=rand();
            if(U<=1/6)
                deltaW=sqrt(3*deltaT);
            elseif(U<=2/6)
                deltaW=-sqrt(3*deltaT);
            else
                deltaW=0;
            end
        end
            
        yBar=y+aCur*deltaT+BCur*deltaW;
        yBarPlus=y+aCur*deltaT+BCur*sqrtDeltaT;
        yBarMinus=y+aCur*deltaT-BCur*sqrtDeltaT;
        
        aBar=a(yBar,t);
        BBarPlus=B(yBarPlus,t);
        BBarMinus=B(yBarMinus,t);
        
        fixedTerm=y+(1/2)*aCur*deltaT...
                   +(1/4)*(BBarPlus+BBarMinus+2*BCur)*deltaW...
                   +(1/4)*(BBarPlus-BBarMinus)*(deltaW.^2-deltaT)/sqrtDeltaT;
        
        y=(1/2)*aBar*deltaT+fixedTerm;
    case 1%The autonomous explicit order 2.0 weak scheme; Equation 1.3 in
          %Chapter 15.1 of [1].
        sqrtDeltaT=sqrt(deltaT);
        if(useGaussian)
            deltaW=sqrtDeltaT*randn(m,1);
        else
            deltaW=zeros(m,1);

            for k1=1:m
                U=rand();
                if(U<=1/6)
                    deltaW(k1)=-sqrt(3*deltaT);
                elseif(U<=2/6)
                    deltaW(k1)=sqrt(3*deltaT);
                else
                    deltaW(k1)=0;
                end
            end
            %The other elements of deltaW remain 0.
        end

        V=zeros(m,m);
        for k1=1:m
            V(k1,k1)=-deltaT;
            for k2=1:(k1-1)
                if(rand()<=1/2)
                    V(k1,k2)=-deltaT;
                else
                    V(k1,k2)=deltaT;
                end
                V(k2,k1)=-V(k1,k2);
            end
        end
            
        yBar=y+aCur*deltaT+BCur*deltaW;
        aBar=a(yBar,t);
        
        UBarPlus=bsxfun(@plus,y,BCur*sqrtDeltaT);
        UBarMinus=bsxfun(@minus,y,BCur*sqrtDeltaT);
        
        BBarUPlus=zeros(d,m,m);
        BBarUMinus=zeros(d,m,m);
        for j=1:m
            BBarUPlus(:,:,j)=B(UBarPlus(:,j),t);
            BBarUMinus(:,:,j)=B(UBarMinus(:,j),t);
        end
        
        BTermSum=zeros(d,1);
        for j=1:m
            RBarPlus=aCur*deltaT+UBarPlus(:,j);
            RBarMinus=aCur*deltaT+UBarMinus(:,j);
            
            BPlusR=B(RBarPlus,t);
            BMinusR=B(RBarMinus,t);
            
            BTermSum=BTermSum+(1/4)*(BPlusR(:,j)+BMinusR(:,j)+2*BCur(:,j))*deltaW(j)...
                               +(1/4)*(BPlusR(:,j)-BMinusR(:,j))*(deltaW(j)^2-deltaT);
                           
            BUTermSum=zeros(d,1);
            for r=1:m
                if(r==j)
                    continue;
                end
                BUTermSum=BUTermSum+(BBarUPlus(:,j,r)+BBarUPlus(:,j,r)-2*BCur(:,j))*deltaW(j)...
                                   +(BBarUPlus(:,j,r)-BBarUPlus(:,j,r))*(deltaW(j)*deltaW(r)+V(r,j))/sqrtDeltaT;
            end
            BTermSum=BTermSum+(1/4)*BUTermSum;
        end
        
        fixedTerm=y+(1/2)*aCur*deltaT+BTermSum;
        y=(1/2)*aBar*deltaT+fixedTerm;
    case 2%The autonomous explicit order 2.0 weak scheme for additive
          %noise; Equation 1.4 in Chapter 15.1 of [1].
        if(useGaussian)
            deltaW=sqrt(deltaT)*randn(m,1);
        else
            deltaW=zeros(m,1);
            U=rand(m,1);
            sel1=(U<=1/6);
            sel2=~sel1&(U<=2/6);
            deltaW(sel1)=sqrt(3*deltaT);
            deltaW(sel2)=-sqrt(3*deltaT);
            %The other elements of deltaW remain 0.
        end
        
        BTerm=BCur*deltaW;
        
        aBar=a(y+aCur*deltaT+BTerm,t);
        
        fixedTerm=y+(1/2)*aCur*deltaT+BTerm;
        y=(1/2)*aBar*deltaT+fixedTerm;
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
