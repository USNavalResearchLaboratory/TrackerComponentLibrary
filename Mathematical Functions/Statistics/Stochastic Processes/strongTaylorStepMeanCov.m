function [yHat,Sigma]=strongTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt,getSqrtCov)
%%STRONGTAYLORSTEPMEANCOV Determine the mean and covariance matrix of a
%           single step of a strong Itô-Taylor expansion of a stochastic
%           differential equation under Itô calculus, such as one could
%           take using the strongStochTaylorStep function. The stochastic
%           differential equations under consideration take the form:
%           dy=a(y,t)*dt+B(y,t)*dW
%           where dW is the differential of an m-dimensional Wiener
%           process and t is time.
%
%INPUTS: y The dX1 initial value of the random process.
%     aCur The value of the drift function a at a(y,t), where t is the
%          current time.
%     BCur The value of the diffusion matrix function B at B(y,t).
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the type of explicit strong Itô-Taylor
%          step to take. Possible values are:
%          0 Use the Euler-Maruyama method from Equation 2.4 in Chapter
%            10.2 of [1]. This is an order 0.5 method. The other algorithms
%            will typically outperform the Euler-Maruyama method.
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
%            p2apypy, papt, and pBpt and that B not depend on y.
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
% getSqrtCov If this parameter is true, then the return parameter Sigma
%          is a square root of the covariance matrix such that Sigma*Sigma'
%          is the covariance matrix. The default if omitted or an empty
%          matrix is passed is false.
%
%OUTPUTS: yHat The dX1 mean of the strong stochastic Taylor step.
%        Sigma The dXd covariance matrix of the strong stochastic Taylor
%              step or, if getSqrtCov=true, the lower-triangular
%              Cholesky-style decompostion of the covariance matrix.
%
%The derivation of the moments comes directly from the fact that with the
%exception of algorithm=2, all of the noise terms in the expressions are
%Gaussian (sometimes correlated) and thus, with some algebra and the known
%central moments of the Gaussian distribution, one can compute the mean and
%covariance matrix of the stochastic Itô-Taylor expansions. Specific
%derivations are given in [2].
%
%For algorithm 2, there are double Itô stochastic integral terms. To
%expand the cross product of the Wiener process and the double stochastic
%integral, Equation 2.16 in Chapter 5.2 of [1] is used. To take the
%expected value of that expression (which is zero), Lemma 5.7.1 in Chapter
%5.7 of [1] is used. Lemma 5.7.2 in Chapter 5.7 of [1] also provides the
%second noncentral moment of the double Itô stochasic integrals.
%
%EXAMPLE 1:
%The values of a, B and their derivatives are just constants in this
%algorithm. Thus, in this example, we just generate random values for those
%rather than using an actual model. Here, we compare the moments returned
%by this function to those obtained by direct simulation of a single step
%for a model with additive noise with algorithm 6.
% numMC=1e4;
% algorithm=6;
% deltaT=1.1;
% d=4;
% m=3;
% y=randn(d,1);
% aCur=randn(d,1);
% BCur=randn(d,m);
% papy=randn(d,d);
% p2apypy=zeros(d,d,d);%Second derivaties are symmetric.
% for j1=1:d
%     for j2=1:j1
%         p2apypy(:,j1,j2)=randn(d,1);
%         p2apypy(:,j2,j1)=p2apypy(:,j1,j2);
%     end
% end
% papt=randn(d,1);
% pBpt=randn(d,m);
% [yHat,Sigma]=strongTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,[],[],papy,p2apypy,papt,pBpt);
% ySamp=zeros(d,numMC);
% for curRun=1:numMC
%     ySamp(:,curRun)=strongStochTaylorStep(y,aCur,BCur,deltaT,algorithm,[],[],papy,p2apypy,papt,pBpt);
% end
% [yHatMC,SigmaMC]=calcMixtureMoments(ySamp);
% norm(yHatMC-yHat)/norm(yHat)%Relative error of the mean.
% norm(SigmaMC-Sigma,'fro')/norm(Sigma,'fro')%Relative error of the covariance.
%The relative errors will tend to be around 1-3% for the mean and the
%covariance matrix norms. Increasing numMC will reduce the error indicating
%that the values are converging.
%
%EXAMPLE 2:
%This is the same as example 1, except for the general Milstein method,
%algorthm=2.
% numMC=1e4;
% algorithm=2;%Order 1.0, general noise
% deltaT=1.1;
% d=4;
% m=3;
% y=randn(d,1);
% aCur=randn(d,1);
% BCur=randn(d,m);
% pBpy=randn(d,m,d);
% [yHat,Sigma]=strongTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,pBpy);
% ySamp=zeros(d,numMC);
% for curRun=1:numMC
%     ySamp(:,curRun)=strongStochTaylorStep(y,aCur,BCur,deltaT,algorithm,pBpy);
% end
% [yHatMC,SigmaMC]=calcMixtureMoments(ySamp);
% norm(yHatMC-yHat)/norm(yHat)%Relative error of the mean.
% norm(SigmaMC-Sigma,'fro')/norm(Sigma,'fro')%Relative error of the covariance.
%The relative errors in the mean and covariance matrix will tend to be
%around 1-3%. Increasing numMC and increasing the number of the terms used
%in the approximation of the Itô integral in strongStochTaylorStep (the p
%input) can improve the estimates.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] D. F. Crouse, "Itô-Taylor expansion moments for continuous-time state
%    propagation," NRL Memo, 2019.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

d=size(BCur,1);
m=size(BCur,2);

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

if(nargin<12||isempty(getSqrtCov))
    getSqrtCov=false;
end

if(algorithm==2||algorithm==4)
    %Algorithm 2 is the Milstein Scheme (order 1.0) for general
    %multivariate noise from from Equation 3.3 of Chapter 10.3 of [1].
    %Algorithm 4 is the Milstein scheme for commutative noise; Equation
    %3.16 of Chapter 10.3 of [1] provides the solution under Stratonovich
    %calculus. However, one can modify Equation 3.3 in the same manner to
    %obtain a solution under Itô calculus.
    
    if(algorithm==4)
        if(d~=m)
            error('This algorithm is only valid for d=m.') 
        end 
    end

    yHat=y+aCur*deltaT;
    
    if(nargout>1)    
        LjVal=LjOperator(BCur,pBpy);
        D=reshape(permute(LjVal,[1,3,2]),[d,m^2]);

        if(getSqrtCov)
            Sigma=tria([sqrt(deltaT)*BCur,sqrt(deltaT^2/2)*D]);
        else
            Sigma=deltaT*(BCur*BCur')+(deltaT^2/2)*(D*D');
        end
    end
    return
end

switch(algorithm)
    case 0%The Euler-Maruyama method; Equation 2.5 of Chapter 10.2 of [1].
        yHat=y+aCur*deltaT;
        
        if(nargout>1)
            if(getSqrtCov)
                Sigma=tria(sqrt(deltaT)*BCur);
            else
                Sigma=deltaT*(BCur*BCur');
            end
        end
    case 1%The Milstein scheme for scalar noise; Equation 3.2 of Chapter
          %10.3 of [1].
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
        yHat=y+aCur*deltaT;

        if(nargout>1)
            D=zeros(d,1);
            for l=1:d
                D=D+BCur(l)*pBpy(:,1,l);
            end
            if(getSqrtCov)
                Sigma=tria([sqrt(deltaT)*BCur,sqrt(deltaT^2/4)*D]);
            else
                Sigma=deltaT*(BCur*BCur')+(deltaT^2/2)*(D*D');
            end
        end
    case 3%Use the Milstein scheme (order 1.0) for diagonal noise from
          %Equation 3.12 of Chapter 10.3 of [1]. This requires
          %pBpy and that d=m.
        if(d~=m)
            error('This algorithm is only valid for d=m.') 
        end
          
        yHat=y+aCur*deltaT;

        if(nargout>1)
            D=zeros(d,d);
            for k=1:d
                D(k,k)=BCur(k,k)*pBpy(k,k,k);
            end

            if(getSqrtCov)
                Sigma=tria([sqrt(deltaT)*BCur,sqrt(deltaT^2/2)*D]);
            else
                Sigma=deltaT*(BCur*BCur')+(deltaT^2/2)*(D*D');
            end
        end
    case 5%Use the strong order 1.5 Taylor method for autonomous scalar
          %problems from Equation 4.1 of Chapter 10.4 of [1]. This requires
          %pBpy,p2Bpypy,papy, p2apypy, that d=m=1, and that a and B don't
          %depend on t.
        if(d~=1||m~=1)
            error('This algorithm is only for scalar problems.')
        end
        
        term1=(aCur*pBpy+(1/2)*BCur^2*p2Bpypy);
        term2=BCur*(BCur*p2Bpypy+pBpy^2);

        c0=y+aCur*deltaT+(1/2)*(aCur*papy+(1/2)*BCur^2*p2apypy)*deltaT^2;

        yHat=c0;
        
        if(nargout>1)
            c1=BCur+term1*deltaT-(1/2)*term2*deltaT;
            c2=(1/2)*BCur*pBpy;
            c3=(1/6)*term2;
            c4=papy*BCur-term1;

            Sigma=(1/3)*deltaT*(3*c1^2+3*c1*(6*c3+c4)*deltaT+deltaT*(6*c2^2+(45*c3^2+9*c3*c4+c4^2)*deltaT));

            if(getSqrtCov)
                Sigma=sqrt(Sigma); 
            end
        end
    case 6%The strong order 1.5 Taylor method for additive noise; Equation
        %4.10 of Chapter 10.4 of [1].
        
        L0=L0Operator(aCur,BCur,papt,papy,p2apypy);

        c0=y+aCur*deltaT+(1/2)*L0*(deltaT^2);

        yHat=c0;
        
        if(nargout>1)
            Lj=LjOperator(BCur,papy);
            C1=BCur+pBpt*deltaT;
            C2=Lj-pBpt;
            
            if(getSqrtCov)
                Sigma=tria([sqrt(deltaT)*(C1+(deltaT/2)*C2),sqrt(deltaT^3/12)*C2]);
            else
                C1C2=C1*C2';
                Sigma=(C1*C1')*deltaT+(deltaT^3/3)*(C2*C2')+(deltaT^2/2)*(C1C2+C1C2');    
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
