function [yHat,Sigma]=weakTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,simplified,pBpy,p2Bpypy,papy,p2apypy,papt,pBpt,getSqrtCov)
%%WEAKTAYLORSTEPMEANCOV Determine the mean and covariance matrix of a
%           single step of a weak It�-Taylor expansion of a stochastic
%           differential equation under It� calculus, such as one could
%           take using the weakStochTaylorStep function. The stochastic
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
% algorithm A parameter specifying the type of explicit weak It�-Taylor
%          step to take. Possible values are:
%          0 Use the order 1.0 simplified weak Euler scheme from Equation
%            1.2 in Chapter 14.1 of [1] for the simplified scheme; the
%            unsimplified scheme is just the Euler-Maruyama method and both
%            have the same first and second moments.
%          1 Use the order 2.0 weak Taylor scheme for autonomous scalar
%            problems from Equation 2.1 in Chapter 14.2 of [1] or Equation
%            2.2 for the simplified scheme. This requires pBpy, p2Bpypy,
%            papy, p2apypy and that d=m=1.
%          2 Use the order 2.0 weak Taylor scheme for general problems from
%            Equation 2.6 in Chapter 14.2 of [1] if simplified=0 or
%            Equation 2.7 for the simplified scheme if simplified=1. This
%            requires pBpy, p2Bpypy, papy, p2apypy, papt, and pBpt.
%          3 Use the order 2.0 weak Taylor scheme for problems with scalar
%            noise from Equation 2.5 in Chapter 14.2 of [1]. This requires
%            pBpy, p2Bpypy, papy, p2apypy, papt, pBpt and that m=1.
% simplified A boolean parameter indicating whether the simplified version
%          of the algorithm should be used. Except for algorithm=0,3,
%          this changes the results.
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
%OUTPUTS: yHat The dX1 mean of the weak stochastic Taylor step.
%        Sigma The dXd covariance matrix of the weak stochastic Taylor
%              step, or, if getSqrtCov=true, the lower-triangular
%              Cholesky-style decompostion of the covariance matrix.
%
%The derivation of the moments proceeds using a number of identities
%involving the products of It� stochastic integrals and the products of the
%Wiener process and an It� stochastic integral. This includes Equation
%2.160 in Chapter 5.2 of [1] for the product of the Wiener process and the
%It� stochastic integral, Lemma 5.7.1 for the expected value of such
%integrals and Lemma 5.7.2 in Chapter 5.7 for the second noncentral moment
%of such integrals. Also, Corollary 5.2.4 in Chapter 5.2 is used to
%evaluate I_{0,0} when necessary. Specific derivations are given in [2].
%
%When getting the square root covariance matrix, for algorithm=1, the
%square root is directly taken. For All other algorithms except the
%unsimplified algorithm 2, a factorted formulation in the tria function is
%used. For the unsimplified algorithm 2, no clear factorization was
%possible, so Sigma is just passed to cholSemiDef.
%
%EXAMPLE:
%Here, we demonstrate algorithm 2 and show that the expected values equal
%those obtained through Monte Carlo simulations of weakStochTaylorStep and
%that the variance of the estimates from the simplified method are not the
%same as those from the unsimplified method.
% numMC=1e4;
% algorithm=2;
% deltaT=1.1;
% d=4;
% m=3;
% p=50;
% y=randn(d,1);
% aCur=randn(d,1);
% BCur=randn(d,m);
% papy=randn(d,d);
% pBpy=randn(d,m,d);
% p2apypy=zeros(d,d,d);%Second derivaties are symmetric.
% p2Bpypy=zeros(d,m,d,d);%Second derivaties are symmetric.
% for j1=1:d
%     for j2=1:j1
%         p2apypy(:,j1,j2)=randn(d,1);
%         p2apypy(:,j2,j1)=p2apypy(:,j1,j2);
%         
%         p2Bpypy(:,:,j1,j2)=randn(d,m);
%         p2Bpypy(:,:,j2,j1)=p2Bpypy(:,:,j1,j2);
%     end
% end
% papt=randn(d,1);
% pBpt=randn(d,m);
% [yHat,Sigma]=weakTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,0,pBpy,p2Bpypy, papy, p2apypy,papt,pBpt);
% [yHat1,Sigma1]=weakTaylorStepMeanCov(y,aCur,BCur,deltaT,algorithm,1,pBpy,p2Bpypy, papy, p2apypy,papt,pBpt);
% ySamp=zeros(d,numMC);
% ySamp1=zeros(d,numMC);
% for curRun=1:numMC
%     ySamp(:,curRun)=weakStochTaylorStep(y,aCur,BCur,deltaT,algorithm,0,pBpy, p2Bpypy, papy, p2apypy,papt,pBpt,p);
%     ySamp1(:,curRun)=weakStochTaylorStep(y,aCur,BCur,deltaT,algorithm,1,pBpy, p2Bpypy, papy, p2apypy,papt,pBpt,p);
% end
% [yHatMC,SigmaMC]=calcMixtureMoments(ySamp);
% [yHatMC1,SigmaMC1]=calcMixtureMoments(ySamp1);
% %Relative errors for the mean and covariance matrix of the unsimplified
% %method.
% norm(yHatMC-yHat)/norm(yHat)
% norm(SigmaMC-Sigma,'fro')/norm(Sigma,'fro')
% %Relative errors for the mean and covariance matrix of the unsimplified
% %method.
% norm(yHatMC1-yHat1)/norm(yHat1)
% norm(SigmaMC1-Sigma1,'fro')/norm(Sigma1,'fro')
% %Relative difference in the true covariance matrices between the simplified
% %and the unsimplified methods.
% norm(Sigma-Sigma1,'fro')/norm(Sigma,'fro')
%One will typically see that the estimation errors are less than 2% for any
%given set of Monte Carlo runs, but the difference in the covariance matrix
%between the simplified and the unsimplified methods is on the order of
%10-20%.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] D. F. Crouse, "It�-Taylor Expansion Moments for Continuous-Time State
%    Propagation," Naval Research Laboratory: Washington, DC., No.
%    NRL/MR/5344--19-9881, 20 Aug. 2019.
%
%December 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

d=size(BCur,1);
m=size(BCur,2);

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

if(nargin<13||isempty(getSqrtCov))
    getSqrtCov=false; 
end

switch(algorithm)
    case 0%The order 1.0 simplified weak Euler scheme; Equation 1.2 in
          %Chapter 14.1 of [1] for the simplified scheme; the unsimplified
          %scheme is just the Euler-Maruyama method.
        yHat=y+aCur*deltaT;
        if(getSqrtCov)
            Sigma=tria(sqrt(deltaT)*BCur);
        else
            Sigma=deltaT*(BCur*BCur');
        end
    case 1%The order 2.0 weak Taylor scheme for scalar problems; Equation
          %2.1 in Chapter 14.2 of [1] or Equation 2.2 for the simplified
          %scheme.
        if(m~=1||d~=1)
            error('This algorithm is only for scalar problems.')
        end
          
        if(simplified==0)%Equation 2.1
            yHat=y+aCur*deltaT+(1/2)*(aCur*papy+(1/2)*p2apypy*BCur^2)*deltaT^2;

            term1=(1/2)*BCur*pBpy;
            term2=aCur*pBpy+(1/2)*p2Bpypy*BCur^2;

            c0=-deltaT*term1;
            c1=BCur+deltaT*term2;
            c2=term1;
            c3=papy*BCur-term2;
            
            Sigma=c0^2+(c1^2+2*c0*c2)*deltaT+(3*c2^2+c1*c3)*deltaT^2+(c3^2/3)*deltaT^3;
        else
            yHat=y+aCur*deltaT+(1/2)*(aCur*papy+(1/2)*p2apypy*BCur^2)*deltaT^2;
            
            term0=(1/2)*BCur*pBpy;
            c0=-term0;
            c1=BCur+(deltaT/2)*(papy*BCur+aCur*pBpy+(1/2)*p2Bpypy*BCur^2);
            c2=term0;
            
            Sigma=c0^2+(c1^2+2*c0*c2)*deltaT+3*c2^2*deltaT^2;
        end
        if(getSqrtCov)
            Sigma=sqrt(Sigma);
        end
    case 2%The order 2.0 weak Taylor scheme for general problems; Equation
          %2.6 in Chapter 14.2 of [1] or Equation 2.7 for the simplified
          %scheme.
        L0a=L0Operator(aCur,BCur,papt,papy,p2apypy);
        L0B=L0Operator(aCur,BCur,pBpt,pBpy,p2Bpypy);
        Lja=LjOperator(BCur,papy);
        LjB=LjOperator(BCur,pBpy);
  
        yHat=y+aCur*deltaT+(1/2)*L0a*deltaT^2;
        if(simplified==0)
            C1=L0B;
            C2=Lja;
            C3=reshape(permute(LjB,[1,3,2]),[d,m^2]);
            
            Sigma=deltaT*(BCur*BCur')+(deltaT^3/3)*(C1*C1'+C2*C2')+(deltaT^3/6)*(C1*C2'+C2*C1')...
                  +(deltaT^2/2)*((C3*C3')+BCur*C1'+BCur*C2'+C1*BCur'+C2*BCur');
            if(getSqrtCov)
                %Unsure of a quadratic factorication of Sigma, so just use
                %cholSemiDef.
                Sigma=cholSemiDef(Sigma,'lower');
            end
        else
            C1=BCur+(deltaT/2)*(L0B+Lja);
            C2=(1/2)*reshape(permute(LjB,[1,3,2]),[d,m^2]);
            
            if(getSqrtCov)
                Sigma=tria([sqrt(deltaT)*C1,sqrt(2)*deltaT*C2]);
            else
                Sigma=deltaT*(C1*C1')+2*deltaT^2*(C2*C2');
            end
        end
    case 3%The order 2.0 weak Taylor scheme for scalar noise; Equation
          %2.5 in Chapter 14.2 of [1].
        if(m~=1)
            error('This algorithm is only for scalar noise.')
        end
          
        L0a=L0Operator(aCur,BCur,papt,papy,p2apypy);
        L0B=L0Operator(aCur,BCur,pBpt,pBpy,p2Bpypy);
        Lja=LjOperator(BCur,papy);
        LjB=LjOperator(BCur,pBpy);
          
        yHat=y+aCur*deltaT+(1/2)*L0a*deltaT^2;
        
        c0=-(deltaT/2)*LjB;
        c1=BCur+deltaT*L0B;
        c2=(1/2)*LjB;
        c3=Lja-L0B;
        
        if(getSqrtCov)
            Sigma=tria([(c0+deltaT*c2),sqrt(2)*deltaT*c2,sqrt(deltaT)*(c1+deltaT/2*c3),(1/2)*sqrt(deltaT^3/3)*c3]);
        else
            Sigma=c0*c0'+(c2*c2')*3*deltaT^2+(c3*c3')*deltaT^3/3+(c1*c1'+c0*c2'+c2*c0')*deltaT+(c1*c3'+c3*c1')*(deltaT^2/2);
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
