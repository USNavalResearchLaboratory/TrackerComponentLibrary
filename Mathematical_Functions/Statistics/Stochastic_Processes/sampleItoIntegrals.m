function [I1Idx,I2Idx]=sampleItoIntegrals(m,deltaT,p,algorithm,deltaW,deltaP,Wj)
%%SAMPLEITOINTEGRALS Return a random samples of Itô integrals with 1 and 2
%          subscripts. These are multiple Itô integrals of 1 (not of a
%          general function). The assumed Wiener process can be
%          multidimensional, in which case this function returns the
%          necessary cross terms too. All values are generated with the
%          same assumed Wiener process. These are useful in Itô-Taylor
%          expansions of stochastic differential equations. Integrals of a
%          jump-diffusion process with specified jumps can also be
%          approximated using algorithm 0.
%
%INPUTS: m The scalar dimensionality of the Wiener process, m>=1.
%   deltaT The time interval over which the integrals are taken.
%        p The number of terms to use in the approximation. The default if
%          omitted or an empty matrix is passed is 100 if algorithm=0,2 and
%          1000 if algorithm=1.
% algorithm An optional parameter specifying the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            expansion based on a component-wise Fourier expansion of the
%            process, which is given in Chapter 10.4 (pages 353-354) of
%            [1]. This algorithm is also provided in [3] in Chapter 5.3,
%            pages 263-264 along with modifications to account for jump 
%            process terms. If deltaW is given, the integral is computed 
%            using that Wiener process difference. If deltaP and W are 
%            given, the Itô integrals involving Poisson measures with mark 
%            independent jump coefficients are also calculated.
%          1 Simulate discrete steps in the underlying Wiener process
%            and use the defintion of the Itô integral to approximate the
%            values using multiple sums.
%          2 Use the algorithm described in [2]. This has better asymptotic
%            convergence as a function of p than 0. However, it is slower
%            than 0 due to operations with large matrices. Note that this
%            algorithm does not provide Ij0 or I0j.
%   deltaW An mX1 vector of the difference between the position of each 
%          Wiener process at the beginning of the current time step and the 
%          position at the end of the current time step. If omitted, a
%          scaled normal random variable is generated according to the time
%          step deltaT.
%   deltaP A scalar representing the number of Poissonian jumps which have
%          occurred during the time interval. If given, W is expected to be
%          defined as well. Otherwise, an error will occur. This input is
%          only used if algorithm=0.
%       Wj An mXdeltaP matrix where each row corresponds to an individual 
%          Wiener process and each column corresponds to the values of the 
%          processes at the time a Poissonian jump occurs. This input is
%          only used if algorithm=0. Note this must contain only the Wiener
%          process values at the jump times, not the process values at the
%          beginning and end of the time interval.
%
%OUTPUTS: I1Idx A structure whose members are all Itô integrals with 1
%               subscript. Members are:
%               I0 for I_0 (scalar) This is the deterministic integral of
%                                   the region 0->deltaT.
%               Ij for I_j (mX1)
%               Im1 for I_{-1} (scalar) This is the number of Poissonian
%                                       jumps which occurred over the time
%                                       interval.
%         I2Idx A structure whose members are all Itô integrals with 2
%               subscripts. Members are:
%               I00 for I_{0,0} (scalar) This is the deterministic second
%                                        integral of the region 0->deltaT.
%               I0j for I_{0,j} (mX1)
%               Ij0 for I_{j,0} (mX1)
%               Ij1j2 for I_{j1,j2} (mXm)
%               Ij1m1 for I_{j1,-1} (mX1)
%               Im1j1 for I_{-1,j1} (mX1)
%               Im1m1 for I_{-1,-1} (scalar)
%
%For algorithm 1, the Itô integral is defined in a manner similar to a
%Riemann-Stieljes integral. That is, the integral of some function f(x),
%when taking N discrete steps, is
%sum_{j=0}^{N-1} f(t_j)*(W_{t_{j+1}}-W_{t_j})
%where W_{t_j} is the value of a Wiener process at discrete time t_j. The
%integral comes as the number of steps goes to Inf. Note that the Itô
%integral evaluates the function at the beginning of the interval. f could
%also be another integral.
%
%A single Itô integral of a the vector Wiener process with m dimensions is:
%I_{j}=sum(deltaW(j,:))
%where j is the selected dimension, deltaW is a vector of the discrete
%steps (W_{t_{j+1}}-W_{t_j}). Next, when doing a double integral where the
%inner integral is deterministic (integral of the constant 1, we have
%I_{j,0}=sum_{i=1}^N t(i-1)*deltaW(j,i) where t(0)=0;
%and when the outer integral is deterministic but the innter one is
%stochastic:
%I_{0,j}=sum_{i1=1}^N sum_{i_2=1}^{i_1}deltaW(j,i1)*dt;
%When dealing with with double stochastic integrals, we have: 
%I_{j1,j2}=sum_{i1=1}^N sum_{i_2=1}^{i_1-1}deltaW(j1,i1)*deltaW(j2,i2)
%More information on stochastic integrals is in [1].
%
%When using Algorithm 2, Equation 4.9 in [2] notes that p can be chosen as
%p>=sqrt(5*m^2*(m-1)/(24*pi^2))/(sqrt(C*deltaT))
%where C is some constant related to the desired accuracy. This compares to
%Equation 4.9 in Chapter 10.4 of [1], for Algorithm 0, where p>=C/deltaT^2.
%
%EXAMPLE:
%Here, we demonstrate that some of the second moments match the value
%returned by ItoIntegral2ndMoment. This might take a few seconds to run.
% m=2;
% deltaT=2;
% p=50;
% numMonteCarlo=1e5;
% Ij1j2Moment=zeros(m,m,numMonteCarlo);
% for curRun=1:numMonteCarlo
%     [~,I2Idx]=sampleItoIntegrals(m,deltaT,p);
%     
%     Ij1j2Moment(:,:,curRun)=I2Idx.Ij1j2.^2;
% end
% meanSim=mean(Ij1j2Moment,3)
% exactSol=ItoIntegral2ndMoment([1,1],[1,1],deltaT)
%The exact solution is 2 (It is the same for all elements) and the meanSim
%will have values that are usually within 0.03 of 2.
%
%REFERENCES:
%[1] P. E. Kloeden and E. Platen, Numerical Solution of Stochastic
%    Differential Equations. Berlin: Springer, 1999.
%[2] M. Wiktorsson, "Joint characteristic function and simultaneous
%    simulation of iterated Itô integrals for multiple independent Brownian
%    motions," The Annals of Applied Probability, vol. 11, no. 2, pp.
%    470-487, May 2001.
%[3] Platen, Eckhard, and Nicola Bruti-Liberati. Numerical solution of 
%    stochastic differential equations with jumps in finance. Vol. 64. 
%    Springer Science & Business Media, 2010.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%Support for jump processes added by Codie T. Lewis July 2019 , Naval
%Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(algorithm))
    algorithm=0; 
end

if(nargin<3||isempty(p))
    if(algorithm==0||algorithm==2)
        p=100;
    else
        p=1000;
    end
end

switch(algorithm)
    case 0%From [1] and [3]
        sqrtDeltaT=sqrt(deltaT);
        
        %Check if the Wiener processes were generated externally.
        if(~exist('deltaW','var')||isempty(deltaW))
            xi=randn(m,1);
            deltaW=sqrtDeltaT*xi;
        else
            xi=deltaW/sqrtDeltaT;
        end
        
        %Equations 6.2.5 in Chapter 6.2 of [3].
        if(exist('deltaP','var')&&~isempty(deltaP))
            I1Idx.Im1=deltaP;
            if(nargout>1)
                I2Idx.Im1m1=0.5*(deltaP^2-deltaP);
                if(exist('Wj','var')&&~isempty(Wj))
                    I2Idx.Ij1m1=sum(Wj,2)-deltaP*deltaW;
                else
                    I2Idx.Ij1m1=-deltaP*deltaW; % A vector of zeros.
                end
                I2Idx.Im1j1=deltaP*deltaW-I2Idx.Ij1m1;
            end
        end

        I0=deltaT;
        %Equation 4.7 in Chapter 10.4 of [1].
        Ij=deltaW;

        I1Idx.I0=I0;
        I1Idx.Ij=Ij;
        
        if(nargout>1)
            I00=deltaT^2/2;
            I2Idx.I00=I00;

            if(m==1)
                %Use the explicit solutions for scalar problems. The
                %solution for Ij0 for scalar problems is from equation 4.3
                %of Chapter 10.4 of [1].
                I2Idx.Ij0=(1/2)*deltaT^(3/2)*(xi+(1/sqrt(3))*randn(1));
                I2Idx.I0j=deltaW*deltaT-I2Idx.Ij0;
                I2Idx.Ij1j2=(1/2)*(deltaW^2-deltaT);
                return;
            end

            zeta=randn(m,p);
            eta=randn(m,p);
            mu=randn(m,1);

            rInv=1./(1:p);
            rhop=1/12-(1/(2*pi^2))*sum(rInv.*rInv);
            a0=-(sqrt(2*deltaT)/pi)*sum(bsxfun(@times,rInv,zeta),2)-2*sqrt(deltaT*rhop)*mu;

            Ij0=(1/2)*deltaT*(deltaW+a0);
            I0j=deltaW*deltaT-Ij0;
            
            Ij1j2=zeros(m,m);
            for j1=1:m
                for j2=1:m
                    if(j1==j2)
                       Ij1j2(j1,j1)=(1/2)*(deltaW(j1)^2-deltaT);
                       continue;
                    end

                    Aj1j2=sum(bsxfun(@times,rInv,(zeta(j1,:).*eta(j2,:)-eta(j1,:).*zeta(j2,:))))/(2*pi);

                    Ij1j2(j1,j2)=(1/2)*deltaT*xi(j1)*xi(j2)-(1/2)*sqrtDeltaT*(xi(j1)*a0(j2)-xi(j2)*a0(j1))+Aj1j2*deltaT;
                end
            end
            I2Idx.I0j=I0j;
            I2Idx.Ij0=Ij0;
            I2Idx.Ij1j2=Ij1j2;
        end
    case 1
        %Time of the intervals.
        intT=(deltaT/p);

        deltaW=sqrt(intT)*randn(p,m);
        
        I0=deltaT;
        I1=sum(deltaW,1).';
        
        I1Idx.I0=I0;
        I1Idx.I1=I1;
        
        if(nargout>1)
            Ij1j2=zeros(m,m);
            
            %The times of the beginnings of the intervals.
            t=linspaceNoEnd(0,deltaT,p).';
            
            I00=deltaT^2/2;
            Ij0=zeros(m,1);
            I0j=zeros(m,1);
            
            for j1=1:m
                Ij0(j1)=sum(t.*deltaW(:,j1));
                wCumSum=cumsum(deltaW(:,j1));
                %Note that intT=t(2)-t(1)
                I0j(j1)=sum(wCumSum)*intT;
            end

            for j1=1:m
                w1=deltaW(:,j1);
                for j2=1:m
                    w2CumSum=[0;cumsum(deltaW(1:(p-1),j2))];
                    Ij1j2(j1,j2)=sum(w1.*w2CumSum);
                end
            end

            I2Idx.I00=I00;
            I2Idx.I0j=I0j;
            I2Idx.Ij0=Ij0;
            I2Idx.Ij1j2=Ij1j2;
        end
    case 2%The algorithm of [2].
        n=p;
        h=deltaT;
        
        I0=deltaT;
        
        %Step 1.
        deltaW=sqrt(h)*randn(m,1);
        I1=deltaW;

        M=m*(m-1)/2;
        K=zeros(M,m^2);
        %Fill in K from the bottom up.
        startRow=1;
        for k=1:(m-1)
            numRows=m-k;
            numColsOver=(k-1)*m+k;

            rowSpan=startRow:(startRow+numRows-1);
            colSpan=(numColsOver+1):(numColsOver+numRows);

            K(rowSpan,colSpan)=eye(numRows,numRows);
            startRow=startRow+numRows;
        end

        %The Kronecker product reversing permutation matrix.
        P=commutationMatrix(m,m);
        Im2=eye(m^2,m^2);
        %The coefficient after (1/k) in the sum is step 2.
        KPi=K*(P-Im2);%Note that P==P' so (Im2-P)*K'=-KPi'.

        %Note that as mentioned in [1]:
        %1) K*Im2*K'==eye(M,M)
        %2) K*P*K'=zeros(M,M);
        %3) K*K'==eye(M,M)

        Im=eye(m,m);
        IM=eye(M,M);

        %Step 2
        ATilde=zeros(M,1);
        for k=n:-1:1
            Xk=randn(m,1);
            Yk=randn(m,1);
            
            ATilde=ATilde+(1/k)*KPi*kron((Yk+sqrt(2/h)*deltaW),Xk);
        end
        ATilde=(h/(2*pi))*ATilde;

        alphaVal=polygamma(1,n+1);
        deltaWOuterProd=deltaW*deltaW';

        SigmaInf=2*IM+(2/h)*KPi*kron(Im,deltaWOuterProd)*KPi';

        deltaWMag2=dot(deltaW,deltaW);
        temp=sqrt(1+deltaWMag2/h);
        %Equation 4.7
        SigmaInfRoot=(SigmaInf+2*temp*IM)/(sqrt(2)*(1+temp));

        G=randn(M,1);
        %Step 3
        ATilde=ATilde+(h/(2*pi))*sqrt(alphaVal)*SigmaInfRoot*G;
        
        %Step 4.
        %Im2 in [2], for the equation below, should have been Im. This is
        %consistent with Equation 2.2 in [2].
        Ij1j2=(deltaWOuterProd-h*Im)/2-reshape(KPi'*ATilde,[m,m]);
        
        I00=h^2/2;

        I1Idx.I0=I0;
        I1Idx.I1=I1;
        
        I2Idx.I00=I00;
        I2Idx.I0j=[];
        I2Idx.Ij0=[];
        I2Idx.Ij1j2=Ij1j2;
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
