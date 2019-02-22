function [y,fixedTerm,randTerms]=strongStochTaylorStep(y,aCur,BCur,deltaT,algorithm,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt,p,ItoAlg)
%%STRONGSTOCHTAYLORSTEP Perform a single step of an explicit strong Itô-
%           Taylor expansion. This integrates a d-dimensional stochastic
%           differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process and t is time. Strong methods converge to an optimal
%           path as the stepsize decreases, not considering finite
%           precision errors.
%
%INPUTS: y The dX1 initial value of the random process.
%     aCur The value of the drift function a at a(y,t), where t is the
%          current time.
%     BCur The value of the diffusion matrix function B at B(y,t).
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Euler-Maruyama method from Equation 2.4 in Chapter 10.2 of
%            [1]. This is an order 0.5 method. The other algorithms will
%            typically outperform the Euler-Maruyama method.
%          1 Use the Milstein scheme (order 1.0) for scalar noise from
%            Equation 3.2 of Chapter 10.3 of [1]. This requires pBpy and
%            that m=1.
%          2 Use the Milstein scheme (order 1.0) for general multivariate
%            noise from Equation 3.3 of Chapter 10.3 of [1]. This requires
%            pBpy.
%          3 Use the Milstein scheme (order 1.0) for diagonal noise from
%            Equation 3.12 of Chapter 10.3 of [1]. This requires
%            pBpy and that d=m.
%          4 Use the Milstein scheme (order 1.0) for commutative noise.
%            Equation 3.16 of Chapter 10.3 of [1] provides the solution
%            under Stratonovich calculus. However, one can modify Equation
%            3.3 in the same manner to obtain a solution under Itô
%            calculus, which is what is done here. This requires pBpy and
%            that d=m. Given diagonal noise, the result should be the same
%            as algorithm 3.
%          5 Use the strong order 1.5 Taylor method for autonomous scalar
%            problems from Equation 4.1 of Chapter 10.4 of [1]. This
%            requires pBpy,p2Bpypy,papy, p2apypy, that d=m=1, and that a
%            and B don't depend on t.
%          6 Use the strong order 1.5 Taylor method for additive noise from
%            Equation 4.10 of Chapter 10.4 of [1]. This requires papy,
%            p2apypy, papt, and pBpt and that B not depend on x.
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
% p, ItoAlg These values are only used if algorithm=2. These values
%          correspond to the inputs p and algorithm of the function
%          sampleItoIntegrals, which is called by this function for
%          necessary terms. Defaults if omitted or an empty matrix is
%          passed are 15 and 0.
%
%OUTPUTS: y The estimated value of the process after taking a step of
%           deltaT. This is a random value.
%   fixedTerm For all algorithms except 6, this is the value of the terms
%           that also appear in the implicit form of the algorithm (Chapter
%           12 of [1]) that do not change during iteration. For algorithm
%           6, this is a structure with members val, L0a and L0j, where val
%           is the value of common terms that do not change and L0a and L0b
%           are the values of the LoOperator and LjOpererator applied to
%           the drift function a. This output is used by the
%           semiImplicitStrongTaylorStep function to avoid recomputing
%           certain values when iterating the implicit equation. 
% randTerms This is a structure whose elements are the random terms used in
%           the approximation. These elements are typically deltaW, the
%           change in the Wiener process over time deltaT, and sometimes
%           deltaZ, which is correlated with deltaW and is an approximation
%           of I_{1,0} (a multiple Itô integral as used in [1]) and Ij1j2,
%           which is a matrix of approximations of double Itô integrals.
%           These values can be useful when this function is used to take
%           an explicit step for an initial estimate in an implicit
%           strong Taylor algorithm.
%
%EXAMPLE 1:
%This is an example of a nonlinear scalar problem with non-additive noise
%where an explicit solution is available as a basis of comparison. In
%Chapter 4.4 of [1], the stochastic differential equation and its implicit
%solution are from Equation 4.40.
% numMC=1e5;
% numSteps=1;
% algorithm=5;
% deltaT=1.1;
% y0=0.5;
% aDrift=@(x,t)((1/3)*x^(1/3));
% BDiff=@(x,t)(x^(2/3));
% explSim=@(W)(y0^(1/3)+(1/3)*W)^3;
% 
% valsSim=zeros(1,numMC);
% valsRK=zeros(1,numMC);
% for curMC=1:numMC
%     W=sqrt(deltaT)*randn(1);
%     valsSim(curMC)=explSim(W);
% 
%     y=y0;
%     for curStep=1:numSteps
%         aCur=aDrift(y);
%         BCur=BDiff(y);
%         
%         dBdx=2/(3*y^(1/3));
%         d2Bdxdx=-(2/(9*y^(4/3)));
%         dadx=1/(9*y^(2/3));
%         d2adxdx=-(2/(27*y^(5/3)));
% 
%         y=strongStochTaylorStep(y,aCur,BCur,deltaT/numSteps,algorithm,dBdx,d2Bdxdx,dadx,d2adxdx);
%     end
%     valsRK(curMC)=y;
% end
% [muSim,PSim]=calcMixtureMoments(valsSim);
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-muSim)./muSim)%Relative mean error.
% abs((PRK-PSim)./PSim)%Relative variance error.
%One will find that the relative mean error is often better than 1% and the
%relative  variance error is often better than 10%.
%
%EXAMPLE 2:
%Here, we consider a multivariate linear dynamic model with multivariate
%noise. The mean and covariance matrix of the sample paths are compared to
%the exact solution.
% algorithm=6;
% numMC=1e4;
% deltaT=1/3;
% numSteps=10;
% x0=[1/4;-12];
% A=[1.1,0.1;
%    -0.2,2.2];
% D=[1.5,-0.4;
%    0.1,1];
% [F,Q]=linDynMod2Disc(deltaT,A,D);
% mu=F*x0;
% P=Q;
% aDrift=@(x,t)(A*x);
% BDiff=@(x,t)(D);
% 
% dadx=A;
% 
% valsRK=zeros(2,numMC);
% for curMC=1:numMC
%     y=x0;
%     for curStep=1:numSteps
%         aCur=aDrift(y);
%         BCur=BDiff(y);
%         y=strongStochTaylorStep(y,aCur,BCur,deltaT/numSteps,algorithm,[],[],dadx);
%     end
%     valsRK(:,curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-mu)./mu)
% abs((PRK-P)./P)
%The error in the estimate of the mean and the variance should be less than
%10%.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Note that changing the number associated with an algorithm will require
%that one also edit the semiImplicitStrongTaylorStep function, which uses
%this function for the initial estimate.

if(nargin<13||isempty(ItoAlg))
    ItoAlg=0;
end

if(nargin<12||isempty(p))
    p=15;
end

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

d=size(aCur,1);%Dimensionality of the state
m=size(BCur,2);%Dimensionality of the noise.

if(nargin<6||isempty(pBpy))
    pBpy=zeros(d,m,d);
end

if(nargin<7||isempty(p2Bpypy))
    p2Bpypy=zeros(d,m,d,d);
end

if(nargin<8||isempty(papy))
    papy=zeros(d,d);
end

if(nargin<9||isempty(p2apypy))
    p2apypy=zeros(d,d,d);
end

if(nargin<10||isempty(papt))
    papt=zeros(d,1);
end

if(nargin<11||isempty(pBpt))
    pBpt=zeros(d,m);
end

switch(algorithm)
    case 0%The Euler-Maruyama method; Equation 2.5 of Chapter 10.2 of [1].
        deltaW=sqrt(deltaT)*randn(m,1);
        fixedTerm=y+BCur*deltaW;

        y=aCur*deltaT+fixedTerm;
        
        if(nargout>2)
            randTerms.deltaW=deltaW;
        end
    case 1%The Milstein scheme for scalar noise; Equation 3.2 of Chapter
          %10.3 of [1].
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
        deltaW=sqrt(deltaT)*randn(1);
        
        milsteinCoeff=zeros(d,1);
        for l=1:d
            milsteinCoeff=milsteinCoeff+(1/2)*BCur(l)*pBpy(:,1,l);
        end    
        
        fixedTerm=y+BCur*deltaW+milsteinCoeff.*(deltaW^2-deltaT);
        
        y=aCur*deltaT+fixedTerm;
        
        if(nargout>2)
            randTerms.deltaW=deltaW;
        end
    case 2%The Milstein scheme for general multivariate noise; Equation 3.3
          %of Chapter 10.3 of [1].

        %Stochastic integrals are needed.
        [I1Idx,I2Idx]=sampleItoIntegrals(m,deltaT,p,ItoAlg);
        
        deltaW=I1Idx.Ij;
        Ij1j2=I2Idx.Ij1j2;

        LjVal=LjOperator(BCur,pBpy);
        
        fixedTerm=y+BCur*deltaW+sum(sum(bsxfun(@times,LjVal,reshape(Ij1j2,[1,m,m])),2),3);
        y=aCur*deltaT+fixedTerm;
        
        if(nargout>1)
            randTerms.deltaW=deltaW;
            randTerms.Ij1j2=Ij1j2;
        end
    case 3%The Milstein scheme for diagonal noise; Equation 3.12 of Chapter
        %10.3 of [1].

        if(d~=m)
            error('This algorithm is only valid for d=m.') 
        end
        
        deltaW=sqrt(deltaT)*randn(m,1);
        fixedTerm=y+BCur*deltaW;
        for k=1:d
            fixedTerm(k)=fixedTerm(k)+(1/2)*BCur(k,k)*pBpy(k,k,k)*(deltaW(k).^2-deltaT);
        end
        
        y=aCur*deltaT+fixedTerm;
        
        if(nargout>2)
            randTerms.deltaW=deltaW;
        end
    case 4%The Milstein scheme for commutative noise; Equation 3.16 of
        %Chapter 10.3 of [1] provides the solution under Stratonovich
        %calculus. However, one can modify Equation 3.3 in the same manner
        %to obtain a solution under Itô calculus.
          
        if(d~=m)
            error('This algorithm is only valid for d=m.') 
        end
          
        deltaW=sqrt(deltaT)*randn(m,1);
        LjBar=LjOperator(BCur,pBpy);
        
        term3=zeros(d,1);
        for j1=1:m
            for j2=1:m
                if(j1==j2)
                    term3=term3+LjBar(:,j2,j1)*(deltaW(j1)^2-deltaT);
                else
                    term3=term3+LjBar(:,j2,j1)*deltaW(j1)*deltaW(j2);
                end
            end
        end
        fixedTerm=y+BCur*deltaW+term3/2;
        
        y=aCur*deltaT+fixedTerm;
        if(nargout>2)
            randTerms.deltaW=deltaW;
        end
    case 5%The strong order 1.5 Taylor method for autonomous scalar
          %problems; Equation 4.1 of Chapter 10.4 of [1].
        if(d~=1||m~=1)
            error('This algorithm is only for scalar problems.')
        end

        U1=randn(1);
        U2=randn(1);
        deltaW=U1*sqrt(deltaT);
        deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);
        
        fixedTerm=y+(1/2)*BCur*pBpy*(deltaW^2-deltaT)...
                +(aCur*pBpy+(1/2)*BCur^2*p2Bpypy)*(deltaW*deltaT-deltaZ)...
                +(1/2)*BCur*(BCur*p2Bpypy+pBpy^2)*((1/3)*deltaW^2-deltaT)*deltaW;
        
        y=aCur*deltaT+BCur*deltaW...
           +papy*BCur*deltaZ+(1/2)*(aCur*papy+(1/2)*BCur^2*p2apypy)*(deltaT^2)...
           +fixedTerm;

        if(nargout>2)
            randTerms.deltaW=deltaW;
            randTerms.deltaZ=deltaZ;
        end
    case 6%The strong order 1.5 Taylor method for additive noise; Equation
        %4.10 of Chapter 10.4 of [1].
        U1=randn(m,1);
        U2=randn(m,1);
        deltaW=U1*sqrt(deltaT);
        deltaZ=(1/2)*deltaT^(3/2)*(U1+(1/sqrt(3))*U2);

        L0=L0Operator(aCur,BCur,papt,papy,p2apypy);
        Lj=LjOperator(BCur,papy);
        
        val=y+BCur*deltaW+pBpt*(deltaW*deltaT-deltaZ);

        y=aCur*deltaT+(1/2)*L0*(deltaT^2)+Lj*deltaZ+val;
        
        if(nargout>1)
            fixedTerm.val=val;
            fixedTerm.L0a=L0;
            fixedTerm.Lja=Lj;
            
            if(nargout>2)
                randTerms.deltaW=deltaW;
                randTerms.deltaZ=deltaZ;
            end
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
