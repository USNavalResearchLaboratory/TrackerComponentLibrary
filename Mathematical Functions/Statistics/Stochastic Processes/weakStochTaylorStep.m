function [y,fixedTerm]=weakStochTaylorStep(y,aCur,BCur,deltaT,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt,p,ItoAlg)
%%WEAKSTOCHTAYLORSTEP Perform a single step of an explicit weak Itô-
%           Taylor expansion. This integrates a d-dimensional stochastic
%           differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process and t is time. As the stepsize used decreases, weak
%           methods converge such that integrals with the random process
%           are a measure are correct (for example, to determine moments).
%           However, they do not converge to the optimal path, unlike
%           strong methods.
%
%INPUTS: y The dX1 initial value of the random process.
%     aCur The value of the drift function a at a(y,t), where t is the
%          current time.
%     BCur The value of the diffusion matrix function B at B(y,t).
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            (weak order 1.0) Euler scheme. For simplified=0, this is just
%            the Euler-Maruyama method from Equation 2.4 in Chapter 10.2 of
%            [1]. For simplified=1, this is Equation 1.2 in Chapter 14.1 of
%            [1]. The other algorithms will typically outperform the Euler
%            scheme.
%          1 Use the order 2.0 weak Taylor scheme for autonomous scalar
%            problems from Equation 2.1 in Chapter 14.2 of [1] if
%            simplified=0 or Equation 2.2 for the simplified scheme if
%            simplfied=1 or 2. If simplified=1, then Gaussians are used for
%            the random variables; otherwise a 3-point distribution is
%            used. This requires pBpy, p2Bpypy, papy, p2apypy and that
%            d=m=1.
%          2 Use the order 2.0 weak Taylor scheme for general problems from
%            Equation 2.6 in Chapter 14.2 of [1] if simplified=0 or
%            Equation 2.7 for the simplified scheme if simplified=1 or 2.
%            If simplified=1, then Gaussians are used for the random
%            variables; otherwise a 3-point distribution is used. This
%            requires papy, papt, and pBpt.
%          3 Use the order 2.0 weak Taylor scheme for problems with scalar
%            noise from Equation 2.5 in Chapter 14.2 of [1]. If
%            simplified=0, then the noise parameters are correlated
%            Gaussians. If simplified=1, then the deltaW noise parameter is
%            Gaussian and deltaZ=(1/2)*deltaW*deltaT, as described before
%            Equation 2.2 in Chapter 14.2 of [1]. If simplified=2, then
%            deltaW is a 3-point distribution and
%            deltaZ=(1/2)*deltaW*deltaT. This requires pBpy, p2Bpypy, papy,
%            p2apypy, papt, pBpt and that m=1.
% simplified A parameter indicating whether the simplfiied form of the
%          algorithms should be used. The values it takes depends on
%          algorithm (see above). The default if omitted or an empty matrix
%          is passed is 1.
%     pBpy A dXmXd matrix of the partial derivatives of B with respect to
%          the d elements of y. pBpy(:,:,i) is the derivative with respect
%          to the ith element of y. If this input is omitted or an empty
%          matrix is passed, then a matrix of zeros will be used.
%  p2Bpypy A dXmXdXd matrix of second partial derivatives of B with respect
%          to the elements of y. p2Bpypy(:,:,i,j) is the derivative with
%          respect to the ith and jth elements of y. If this input is
%          omitted or an empty matrix is passed, then a matrix of zeros
%          will be used.
%     papy A dXd matrix of the partial derivatives of B with respect to the
%          d elements of y. papy(:,i) is the derivative with respect to the
%          ith element of y. If this input is omitted or an empty matrix is
%          passed, then a matrix of zeros will be used.
%  p2apypy A dXdXd matrix of second partial derivatives of B with respect
%          to the elements of y. p2apypy(:,i,j) is the derivative with
%          respect to the ith and jth elements of y. If this input is
%          omitted or an empty matrix is passed, then a matrix of zeros
%          will be used.
%     papt The dX1 partial derivative vector of the drift function a with
%          respect to time t. If this input is omitted or an empty matrix
%          is passed, then a matrix of zeros will be used.
%     pBpt The dXm partial derivative vector of the diffusion function a
%          with respect to time t. If this input is omitted or an empty
%          matrix is passed, then a matrix of zeros will be used.    
% p, ItoAlg These values are only used if algorithm=2 and simplified=0.
%           These values correspond to the inputs p and algorithm of the
%           function sampleItoIntegrals, which is called by this function
%           for necessary terms. Defaults if omitted or an empty matrix is
%           passed are 15 and 0.
%
%OUTPUTS: y The estimated value of the process after taking a step of
%           deltaT. This is a random value.
% fixedTerm This output is used by the implicitWeakTaylorStep function. It
%           is the fixed term in the iterations when using the implicit
%           version of the algorithm. It is an empty matrix for the
%           simplifications that do not map to algorithms in
%           implicitWeakTaylorStep.
%
%EXAMPLE 1:
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its solution
%are from Equation 4.40.
% numMC=1e5;
% numSteps=1;
% algorithm=1;
% simplified=2;
% deltaT=1.1;
% y0=0.5;
% aDrift=@(x,t)((1/3)*x^(1/3));
% BDiff=@(x,t)(x^(2/3));
% explSim=@(W)(y0^(1/3)+(1/3)*W)^3;
% 
% %Take the mean and covariance matrix of the explicit solution using
% %quadrature integration.
% [xi,w]=quadraturePoints1D(6);%2*6-1=11th order.
% numPts=length(w);
% xi=sqrt(deltaT)*xi;
% for k=1:numPts
%     xi(:,k)=explSim(xi(:,k));
% end
% [muCub,PCub]=calcMixtureMoments(xi,w);
% 
% valsRK=zeros(1,numMC);
% for curMC=1:numMC
%     y=y0;
%     for curStep=1:numSteps
%         aCur=aDrift(y);
%         BCur=BDiff(y);
%         
%         pBpy=2/(3*y^(1/3));
%         p2Bpypy=-(2/(9*y^(4/3)));
%         papy=1/(9*y^(2/3));
%         p2apypy=-(2/(27*y^(5/3)));
% 
%         y=weakStochTaylorStep(y,aCur,BCur,deltaT/numSteps,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy);
%     end
%     valsRK(curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-muCub)./muCub)%Relative mean error.
% abs((PRK-PCub)./PCub)%Relative variance error.
%One will see that the error in the mean and variance estimates are
%typically below 5%.
%
%%EXAMPLE 2:
%Here, we consider a multivariate linear dynamic model with multivariate
%noise. The mean and covariance matrix of the sample paths are compared to
%the exact solution.
% algorithm=2;
% simplified=1;
% numMC=1e4;
% deltaT=1/3;
% numSteps=5;
% x0=[1/4;-12];
% A=[1.1,0.1;
%    -0.2,2.2];
% D=[1.5,-0.4;
%    0.1,1];
% d=size(D,1);
% [F,Q]=linDynMod2Disc(deltaT,A,D);
% mu=F*x0;
% P=Q;
% aDrift=@(x,t)(A*x);
% BDiff=@(x,t)(D);
% 
% papy=A;
% valsRK=zeros(d,numMC);
% for curMC=1:numMC
%     y=x0;
%     for curStep=1:numSteps
%         aCur=aDrift(y);
%         BCur=BDiff(y);
%         y=weakStochTaylorStep(y,aCur,BCur,deltaT/numSteps,algorithm,simplified,[],[],papy);
%     end
%     valsRK(:,curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-mu)./mu)%Relative mean error.
% abs((PRK-P)./P)%Relative variance error.
%One will see that the error in the mean and variance terms are typically
%below 10%.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<14||isempty(ItoAlg))
    ItoAlg=0;
end

if(nargin<13||isempty(p))
    p=15;
end

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

if(nargin<6||isempty(simplified))
    simplified=1;
end

d=size(aCur,1);%Dimensionality of the state
m=size(BCur,2);%Dimensionality of the noise.

if(nargin<7||isempty(pBpy))
    pBpy=zeros(d,m,d);
end

if(nargin<8||isempty(p2Bpypy))
    p2Bpypy=zeros(d,m,d,d);
end

if(nargin<9||isempty(papy))
    papy=zeros(d,d);
end

if(nargin<10||isempty(p2apypy))
    p2apypy=zeros(d,d,d);
end

if(nargin<11||isempty(papt))
    papt=zeros(d,1);
end

if(nargin<12||isempty(pBpt))
    pBpt=zeros(d,m);
end

switch(algorithm)
    case 0%The order 1.0 simplified weak Euler scheme; Equation 1.2 in
          %Chapter 14.1 of [1] for the simplified scheme; the unsimplified
          %scheme is just the Euler-Maruyama method.
        sqrtDeltaT=sqrt(deltaT);
        if(simplified==0)%Euler-Maruyama
            deltaW=sqrtDeltaT*randn(m,1);
        else%Simplified weak Euler scheme.
            deltaW=rand(m,1);
            sel=(deltaW<=0.5);
            deltaW(sel)=-sqrtDeltaT;
            deltaW(~sel)=sqrtDeltaT;
        end
        
        fixedTerm=BCur*deltaW+y;
        y=aCur*deltaT+fixedTerm;
    case 1%The order 2.0 weak Taylor scheme for scalar problems; Equation
          %2.1 in Chapter 14.2 of [1] or Equation 2.2 for the simplified
          %scheme.
        if(m~=1||d~=1)
            error('This algorithm is only for scalar problems.')
        end
          
        if(simplified==0)%Equation 2.1
            U1=randn(1);
            U2=randn(1);
            deltaW=U1*sqrt(deltaT);
            deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);
            
            y=y+aCur*deltaT+BCur*deltaW+(1/2)*BCur*pBpy*(deltaW.^2-deltaT)...
                +papy*BCur*deltaZ+(1/2)*(aCur*papy+(1/2)*p2apypy*BCur.^2)*(deltaT^2)...
                +(aCur*pBpy+(1/2)*p2Bpypy*BCur^2)*(deltaW*deltaT-deltaZ);
            fixedTerm=[];
        else
            if(simplified==1)
                deltaW=sqrt(deltaT)*randn(1);
            else
                U=rand();
                if(U<=1/6)
                    deltaW=-sqrt(3*deltaT);
                elseif(U<=2/6)
                    deltaW=sqrt(3*deltaT);
                else
                    deltaW=0;
                end
            end
            
            fixedTerm=y+BCur*deltaW+(1/2)*BCur*pBpy*(deltaW.^2-deltaT)...
                       +(1/2)*(aCur*pBpy+(1/2)*p2Bpypy*BCur^2)*deltaW*deltaT;
            
            y=aCur*deltaT+...
                +(1/2)*(papy*BCur)*deltaW*deltaT...
                +(1/2)*(aCur*papy+(1/2)*p2apypy*BCur.^2)*deltaT^2+fixedTerm;
        end
    case 2%The order 2.0 weak Taylor scheme for general problems; Equation
          %2.6 in Chapter 14.2 of [1] or Equation 2.7 for the simplified
          %scheme.
        L0a=L0Operator(aCur,BCur,papt,papy,p2apypy);
        L0B=L0Operator(aCur,BCur,pBpt,pBpy,p2Bpypy);
        Lja=LjOperator(BCur,papy);
        LjB=LjOperator(BCur,pBpy);
          
        if(simplified==0)%Equation 2.6
            %Stochastic integrals are needed.
            [I1Idx,I2Idx]=sampleItoIntegrals(m,deltaT,p,ItoAlg);
            deltaW=I1Idx.Ij;
            I0j=I2Idx.I0j;
            Ij0=I2Idx.Ij0;
            Ij1j2=I2Idx.Ij1j2;

            endBSum=zeros(d,1);
            for j1=1:m
                for j2=1:m
                    endBSum=endBSum+LjB(:,j2,j1)*Ij1j2(j1,j2);
                end
            end

            fixedTerm=[];
            y=aCur*deltaT+(1/2)*L0a*(deltaT^2)...
                +Lja*Ij0+y+endBSum+BCur*deltaW+L0B*I0j;
        else
            if(simplified==1)
                deltaW=sqrt(deltaT)*randn(m,1);
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

            endBSum=zeros(d,1);
            for j1=1:m
                for j2=1:m
                    endBSum=endBSum+LjB(:,j2,j1).*(deltaW(j1)*deltaW(j2)+V(j1,j2));
                end
            end
            endBSum=endBSum/2;
            
            fixedTerm=y+endBSum+BCur*deltaW+(1/2)*L0B*deltaW*deltaT;
            
            y=aCur*deltaT+(1/2)*L0a*(deltaT^2)...
                +(1/2)*Lja*deltaW*deltaT+fixedTerm;
        end
    case 3%The order 2.0 weak Taylor scheme for scalar noise; Equation
          %2.5 in Chapter 14.2 of [1]. The simplified version uses the same
          %approximation for deltaW and deltaZ that were used in Equation
          %2.2.
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
          
        L0a=L0Operator(aCur,BCur,papt,papy,p2apypy);
        L0B=L0Operator(aCur,BCur,pBpt,pBpy,p2Bpypy);
        Lja=LjOperator(BCur,papy);
        LjB=LjOperator(BCur,pBpy);
        if(simplified==0)
            U1=randn(1);
            U2=randn(1);
            deltaW=U1*sqrt(deltaT);
            deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);
        else
            if(simplified==1)
                deltaW=sqrt(deltaT)*randn(1);
            else
                U=rand();
                if(U<=1/6)
                    deltaW=-sqrt(3*deltaT);
                elseif(U<=2/6)
                    deltaW=sqrt(3*deltaT);
                else
                    deltaW=0;
                end
            end

            deltaZ=(1/2)*deltaW*deltaT;
        end

        y=y+aCur*deltaT+BCur*deltaW+(1/2)*LjB*(deltaW^2-deltaT)...
            +(1/2)*L0a*deltaT^2+L0B*(deltaW*deltaT-deltaZ)...
            +Lja*deltaZ;

        fixedTerm=[];
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
