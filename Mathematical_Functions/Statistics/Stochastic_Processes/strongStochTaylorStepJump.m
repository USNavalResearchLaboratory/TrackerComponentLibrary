function [approx,drivers]=strongStochTaylorStepJump(y,aCur,BCur,cCur,lambda,deltaT,algorithm,pBpy,pcpy,deltaP,W,Ts)
%%STRONGSTOCHTAYLORSTEPJUMP Perform a step of an explicit strong 
%           Ito-Taylor expansion for a Levy process. This integrates a 
%           d-dimensional stochastic differential equation of the form:
%           dy=a(y,t)*dt+B(y,t)*dW+c(y)*dP
%           where dW is the differential of an m-dimensional Wiener
%           process, dP is the differential of a Poisson jump process and 
%           t is time. Strong methods converge to an optimal path as the 
%           stepsize decreases, not considering finite precision errors.
%           The Poisson process intensities for all jumps are assumed to be
%           time-invariant and mark independent during the given time 
%           interval with the exception of algorithm 1 (see details below). 
%           This means that the value of the jump term coefficient does not
%           change based on how many jumps occur in a given time step 
%           (essentially imposing a piecewise constant assumption on the 
%           c(y) term in the generated solution).
%
%INPUTS: y The dX1 initial value of the random process.
%     aCur The drift coefficient function which takes as input the vector y
%          and outputs a dX1 vector.
%     BCur The diffusion coefficient function which takes as input the
%          vector y and outputs a dXm matrix where m is the number of
%          independent Wiener processes driving Y.
%     cCur The jump coefficient function which takes as input the vector y
%          and outputs a dX1 vector. For algorithm 1, three inputs will be
%          given including the position y, jump time, and jump index with 1
%          being the first jump for the current time step.
%   lambda The frequency parameter for the Poisson distribution for a time
%          increment of 1.
%   deltaT The time increment over which the step is taken.
% algorithm A parameter specifying the algorithm to use. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Euler-Maruyama method from equation 6.1.27, page 277 in 
%            Chapter 6.1 of [1]. This is an order 0.5 method.
%          1 Use the Euler-Maruyama method from equation 6.1.22, page  277 
%            in Chapter 6.1 of [1]. This is an order 0.5 method. It does
%            not assume mark independence (as the other methods do). This
%            means the jump term coefficient can vary according to the
%            index of a jump during the step. As noted above, the cCur
%            function is required to accept 3 inputs, the process value,
%            jump time, and the jump index in that order. Otherwise an
%            error referencing insufficient inputs will occur.
%          2 Use the Milstein scheme (order 1) for general multivariate
%            noise with jumps from equation 6.2.9 of Chapter 6.2 of [1]. 
%            This requires pBpy and pcpy. If W if also given, be sure to
%            include Ts if the process values begin at a time other than 0.
%          3 Use the order 1 Taylor scheme for mark indpenedent jump 
%            coefficients given in equation 6.3.15 of Chapter 6.3, page 
%            289 of [1], assuming the diffusion commutivity condition given
%			 in equation 5.3.25 of Chapter 5, page 258 of [1], and the jump
%            commutivity condition, given in equation 6.3.7 of Chapter 6.3,
%            page 287 of [1].
%     pBpy A dXmXd matrix of the partial derivatives of B with respect to
%          the d elements of y. pBpy(:,:,i) is the derivative with respect
%          to the ith element of y. If this input is omitted or an empty
%          matrix is passed, then a matrix of zeros will be used.
%     pcpy A dXd matrix of the partial derivatives of c with respect to
%          the d elements of y. pcpy(:,i) is the derivative with respect to
%          the ith element of y. If this input is omitted or an empty
%          matrix is passed, then a matrix of zeros will be used.
%   deltaP A scalar representing the number of Poissonian jumps which have
%          occurred during the time interval. If omitted, a value is
%          generated internally.
%        W An mXdeltaP matrix where each row corresponds to an individual 
%          Wiener process and each column corresponds to the values of the 
%          processes at the time a Poissonian jump occurs. If omitted, the
%          matrix is generated internally. It is assumed that the given
%          matrix includes the initial value at the beginning of the time
%          step in column 1, then the ordered values for each jump
%          occurring in the time interval, and finally the values at the
%          end of the time step in the last column.
%       Ts The time at the beginning of the step. If not given, defaults to
%          0, indicating that the Wiener process jump times will only be
%          scaled to the interval (0,deltaT].
%
%OUTPUTS: approx The estimated value of the process after taking a step of
%                 deltaT. This is a random value.
%        drivers A structure with members containing the driving processes
%                 used in computations. Note that not all memebers are
%                 needed for some methods and therefore are not created.
%                 Refer to the last lines of the method you are using to see
%                 which members are created. Members are:
%                 dW (mX1) The Wiener process differences over the given
%                    time interval.
%                 dP (scalar) The number of Poisson jumps occuring during
%                    the given time interval.
%                  W (mXdeltaP+2) The Wiener process values used in the
%                     calculation.
%                 jT (1XdeltaP+2) The jump times used in the calculation.
%
%EXAMPLE 1: A simple plot of several realized paths dominated by the jump
%         process. The expected value of a compound Poisson process is the
%         product of the expected number of jumps times the expected value
%         of the jump size (since these are both independent variables).
%         Here we fixed the jump size to 10 and the number of expected
%         jumps is 10*0.25 = 2.5. So the average should be around 25. 100
%         paths, with the given presets, are usually good enough to get 
%         within the interval [20,30] and takes about a second to run with
%         plotting. The variance should be the expected number of jumps
%         times the expectation of the jump size squared. Therefore, a
%         variance of 250 is expected. With the presets, the simulated
%         variance will be within [190,310] most of the time.
% 
% numpaths = 100; % Number of paths to simulate
% Tf = 10; % End time (assume start time is 0
% dt = 0.1; % Time step
% steps = ceil(Tf/dt)+1; % Number of steps to perform assuming time [0,100]
% 
% hold on
% t = 0:dt:Tf;
% X = zeros(steps,numpaths);
% for n = 1:numpaths
%     X(1) = 0;
%     for i = 2:steps
%         X(i,n) = strongStochTaylorStepJump(X(i-1,n),@(x)0,@(x)1,@(x)10,0.25,dt,0);
%     end
%     %plot(t,X(:,n))
% end
% hold off
% 
% [avg,P] = calcMixtureMoments(X(end,:));
% fprintf(strcat('The average of all %d processes is %.5f \n',...
%         'with a variance of %.5f.\n'),numpaths,avg,P);
%
%EXAMPLE 2: We test that the extrapolated paths from one of the algorithms
%           achieves similar mixture statistics to those generated by the 
%           exact Merton model solution. This can take a couple seconds to 
%           run.
% 
%For this particular example, the order 1 Taylor method and the order 1
% derivative-free method in strongRungeKStep.m should return similar values
% values since the order 1 expansion of the Merton model is identical for
% both methods as [1] notes in Chapter 7.1, page 313 comparing equation
% 7.1.9 on page 312 to equation 6.3.6 in Chapter 6.3, page 287. This is due
% to the particular choices of diffusion and jump coefficients.
% 
%The Merton model is discussed in detail in Chapter 1.7, pages 50-52 and 
% Chapter 9.6, page 414 of [1]. See equations 9.6.3 and 9.6.4 and the
% discussion on page 414. Note we simplified the jump process term's 
% coefficient to c(y)=psi*y as suggested in Chapter 7.1, page 313.  
% 
% numsamples = 10000; %Number of samples to collect for moments
%                     %calculation
% 
% %Model parameters
% X0=1; % Initial path value
% mue=0.05; % Drift mean
% sig=0.2; % Wiener process standard deviation
% psi=-0.25; % Jump intensity coefficient
% T=1; % End time starting from 0
% lambda=0.3; % Poisson jump frequency for unit time interval
% 
% %Define Discretization Coefficients and Derivatives
% % Assume a time invariant SDE of the form: 
% % dY_{n+1} = a(Y_n) dt + b(Y_n) dW + c(Y_n) dP
% a = @(x) mue*x;
% B = @(x) sig*x;
% c = @(x) psi*x;
% pBpy = sig;
% pcpy = psi;
% 
% disp(' ')
% disp('Computing samples...')
% 
% %Create sample vectors
% Euler = zeros([1,numsamples]);
% Taylor1 = zeros([1,numsamples]);
% X = zeros([1,numsamples]);
% for i = 1:numsamples
%     %Get approximation samples
%     Euler(i) = strongStochTaylorStepJump(X0,a,B,c,lambda,T,0);
%     Taylor1(i) = strongStochTaylorStepJump(X0,a,B,c,lambda,T,2,pBpy,pcpy);
%     
%     %Compute Merton Process Values using Exact Solution 
%     % Equation 1.7.46 from [1]
%     X(i) = X0.*exp((mue-0.5*sig.^2)*T+sig.*randn()*sqrt(T)+psi*PoissonD.rand(1,lambda*sqrt(T)));
% end
% 
% %Compute moments and print to screen
% [avg,P] = calcMixtureMoments(Euler);
% fprintf(strcat('The Euler samples have a mixture mean of %.10f and\n',...
%      'a mixture variance of %.10f.\n'),avg,P);
% [avg,P] = calcMixtureMoments(Taylor1);
% fprintf(strcat('The Taylor1 samples have a mixture mean of %.10f and\n',...
%      'a mixture variance of %.10f.\n'),avg,P);
% [avg,P] = calcMixtureMoments(X);
% fprintf(strcat('\nThe true samples have a mixture mean of %.10f and\n',...
%      'a mixture variance of %.10f.\n'),avg,P);
%
%EXAMPLE 3: Here, we consider a multivariate linear dynamic model with 
%           multivariate noise. The mean and covariance matrix of the 
%           sample paths are compared to the exact solution.
% algorithm=0;
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
% aDrift=@(x)(A*x);
% BDiff=@(x)(D);
% 
% valsRK=zeros(2,numMC);
% for curMC=1:numMC
%     y=x0;
%     for curStep=1:numSteps
%         y=strongStochTaylorStepJump(y,aDrift,BDiff,@(x)0,0,deltaT/numSteps,algorithm);
%     end
%     valsRK(:,curMC)=y;
% end
% [muRK,PRK]=calcMixtureMoments(valsRK);
% abs((muRK-mu)./mu)%Relative mean error.
% abs((PRK-P)./P)%Relative variance error.
% 
% algorithm=1;
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
%REFERENCES:
%[1] Platen, Eckhard, and Nicola Bruti-Liberati. Numerical solution of 
%    stochastic differential equations with jumps in finance. Vol. 64. 
%    Springer Science & Business Media, 2010.
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7 || isempty(algorithm))
    algorithm=0;
end

ay = aCur(y);
By = BCur(y);
[d,m] = size(By);
if(algorithm~=1)
    cy = cCur(y);
    cyc = cCur(y+cy);
    Byc = BCur(y+cy);
end

if(~exist('pBpy','var') || isempty(pBpy))
    pBpy = zeros([d,m,d]);
end

if(~exist('pcpy','var') || isempty(pcpy))
    pcpy = zeros([d,1]);
end

if(~exist('Ts','var')||isempty(Ts))
    Ts = 0;
end

switch(algorithm)
    case 0 %The Euler-Maruyama method; Equation 6.1.27, page 277 of [1].
        if(~exist('W','var')||isempty(W))
            deltaW = sqrt(deltaT)*randn([m,1]);
        else
            deltaW = W(:,end)-W(:,1);
        end
        
        if(~exist('deltaP','var')||isempty(deltaP))
            deltaP = PoissonD.rand(1,lambda*deltaT);
        end
        
        wiener = y+ay*deltaT+By*deltaW;
        
        poiTerm = cy*deltaP;

        approx = wiener+poiTerm;
        
        %Drivers members
        drivers.dW = deltaW;
        drivers.dP = deltaP;
    case 1 %The mark dependent Euler-Maruyama method; Equation 6.1.22, 
           % page 277 of [1]. Note that cCur must be a function of the
           % process value (first input), the jump time (second input),
           % and the index of the jump for the current time step (third
           % input).
        if(~exist('deltaP','var')||isempty(deltaP))
            deltaP = PoissonD.rand(1,lambda*deltaT);
        end
        
        jumpTimes = sort([0,rand([1,deltaP]),1]*(deltaT+Ts));
        
        if(~exist('W','var')||isempty(W))
            W = zeros([m,deltaP+2]);
            for i = 1:m
                W(i,:) = [0,sqrt(diff(jumpTimes)).*randn([1,deltaP+1])];
            end
            W = cumsum(W,2);
        end
        
        deltaW = W(:,end)-W(:,1);
        
        wiener = y+ay*deltaT+By*deltaW;
        
        poiTerm = zeros(size(y));
        for i = 2:length(W)-1
            poiTerm = poiTerm + cCur(y,jumpTimes(i-1),i-1);
        end

        approx = wiener+poiTerm;
        
        %Drivers members
        drivers.dW = deltaW;
        drivers.dP = deltaP;
        drivers.W = W;
        drivers.jT = jumpTimes;     
    case 2 %The Milstein scheme for mark independent jump coefficients; 
           % Equation 6.2.9 of Chapter 6.2, page 283 of [1]
        if(~exist('deltaP','var')||isempty(deltaP))
            deltaP = PoissonD.rand(1,lambda*deltaT);
        end
        
        jumpTimes = sort([0,rand([1,deltaP]),1]*(deltaT+Ts));
        
        if(~exist('W','var')||isempty(W))
            W = zeros([m,deltaP+2]);
            for i = 1:m
                W(i,:) = [0,sqrt(diff(jumpTimes)).*randn([1,deltaP+1])];
            end
            W = cumsum(W,2);
        end
        
        deltaW = W(:,end)-W(:,1);
        [I1Idx,I2Idx] = sampleItoIntegrals(m,deltaT,[],0,deltaW,deltaP,W(:,2:end-1));

        noisenoise = zeros(size(y));
        for j2 = 1:m
            noisenoise = noisenoise+squeeze(pBpy(:,j2,:))*By*I2Idx.Ij1j2(:,j2);
        end
        noisejump = pcpy.'*By*I2Idx.Ij1m1;
        jumpnoise = (Byc-By)*I2Idx.Im1j1;

        jumpjump = (cyc-cy)*I2Idx.Im1m1;

        approx = y+ay*I1Idx.I0+sum(By*I1Idx.Ij,2)+...
                       cy*I1Idx.Im1+noisenoise+noisejump+...
                       jumpnoise+jumpjump;
        
        %Drivers members
        drivers.dW = deltaW;
        drivers.dP = deltaP;
        drivers.W = W;
        drivers.jT = jumpTimes;
    case 3 %The order 1 Taylor scheme for mark indpenedent jump 
           % coefficients assuming the diffusion commutivity condition;
           % from equation 6.3.15 of Chapter 6.3, page 289 of [1]
        if(~exist('deltaP','var')||isempty(deltaP))
            deltaP = PoissonD.rand(1,lambda*deltaT);
        end
        
        if(~exist('W','var')||isempty(W))
            deltaW = sqrt(deltaT)*randn([m,1]);
        else
            deltaW = W(:,end)-W(:,1);
        end

        noisenoise = zeros(size(y));
        for j2 = 1:m
            noisenoise = noisenoise+squeeze(pBpy(:,j2,:))*By*(deltaW*deltaW(j2)-deltaT*ones([m,1]));
        end
        noisejump = sum((Byc-By)*deltaP*deltaW,2);

        jumpjump = (cyc-cy)*(deltaP^2-deltaP);
        
        approx = y+ay*deltaT+By*deltaW+...
                       cy*deltaP+0.5*noisenoise+noisejump+0.5*jumpjump;
        
        %Drivers members
        drivers.dW = deltaW;
        drivers.dP = deltaP;
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
