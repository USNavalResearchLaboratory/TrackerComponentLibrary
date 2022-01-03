function demoMertonModel()
%%DEMOMERTONMODEL Demonstrates scenario simulation of two independent
%             Merton processes using the strongStochTaylorStepJump
%             function. A plot comparing the approximated paths to the true
%             path is generated. In order to demonstrate ideal results, the
%             same collection of driving processes are used for both the
%             true path and the simulations. The first line may be
%             uncommented and modified if reproducible plots are desired.
%
%The Merton model is discussed in detail in Chapter 1.7, pages 50-52 and 
%Chapter 9.6, page 414 of [1]. See equations 9.6.3 and 9.6.4 and the
%discussion on page 414. Note we simplified the jump process term's 
%coefficient to c(y)=psi*y as suggested in Chapter 7.1, page 313.  
%
%REFERENCES:
%[1] Platen, Eckhard, and Nicola Bruti-Liberati. Numerical solution of 
%    stochastic differential equations with jumps in finance. Vol. 64. 
%    Springer Science & Business Media, 2010.
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %rng(16,'twister'); % Make the example reproducible
    close all;

    %Define Parameters
    X0=[1;1]; % Initial path value
    mue=[0.05;0.05]; % Drift mean
    sig=[0.2;0.2]; % Wiener process standard deviations
    psi=[-0.25;-0.25]; % Jump intensity coefficients
    T=10; % End time starting from 0
    lambda=0.3; % Poisson jump frequency for unit time interval

    %Define Discretization Coefficients and Derivatives
    % Assume a time invariant SDE of the form: 
    % dY_{n+1} = a(Y_n) dt + b(Y_n) dW + c(Y_n) dP
    a = @(x) diag(mue)*x;
    B = @(x) diag(sig)*diag(x);
    c = @(x) diag(psi)*x;
    pBpy = cat(3,sig(1)*[1 0; 0 0],sig(2)*[0 0; 0 1]);
    pcpy = diag(psi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%% Exact Solution %%%%%%%%%%%%%%%%%%%%%
    %The Merton model is discussed in detail on p.50-52 and p.414 of [1].
    % We simplified the compound jump process to have constant jump size.
    % See equations 9.6.3 and 9.6.4.

    dt = 0.001; % True path time step
    t = 0:dt:T;
    X = zeros([2,length(t)]); % A vector to contain the true path values
    X(:,1) = X0; % Initial condition

    %Generate Poisson samples as 2 independent compound Poisson processes
    CP = [0,genMarkedPointProc(lambda,dt,0+dt,T,0)];

    %Simulate 2 independent Wiener Processes
    W = [[0;0],sqrt(dt).*randn([2,length(t)-1])];
    W = cumsum(W,2);

    %Compute Merton Process Values using Exact Solution from [1]
    for i = 2:length(X)
        X(:,i) = X0.*exp((mue-0.5*sig.^2)*t(i)+sig.*W(:,i)+psi*CP(i));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%% Approximations %%%%%%%%%%%%%%%%%%%%%

    %Define Approximation Parameters
    dtapprox = 0.5; % Approximation time step.
    tapprox = 0:dtapprox:T;
    numruns = 50; % Number of times to perform the random simulation steps.
    
    XEuler = zeros([2,length(tapprox)]);
    XEuler(:,1) = X0;
    XTaylor = zeros([2,length(tapprox)]);
    XTaylor(:,1) = X0;
    XDerFree = zeros([2,length(tapprox)]);
    XDerFree(:,1) = X0;
    
    for i = 2:length(tapprox)
        [deltaP,Wj] = getDrivers(W,CP,dtapprox,dt,i);
        
        %Get Euler Approximation
        for run = 1:numruns
            XEuler(:,i) = XEuler(:,i) + strongStochTaylorStepJump(...
                XEuler(:,i-1),a,B,@(x,k,idx)c(x),lambda,dtapprox,1,[],[],deltaP,Wj);
        end
        XEuler(:,i) = XEuler(:,i)/numruns;
            
        %Get Taylor Approximation
        for run = 1:numruns
            XTaylor(:,i) = XTaylor(:,i) + strongStochTaylorStepJump(...
                XTaylor(:,i-1),a,B,c,lambda,dtapprox,3,pBpy,pcpy,deltaP,Wj,tapprox(i));
        end
        XTaylor(:,i) = XTaylor(:,i)/numruns;
        
        %Get Derivative Free Taylor Approximation
        for run = 1:numruns
            XDerFree(:,i) = XDerFree(:,i) + strongRungeKStepJump(...
                XDerFree(:,i-1),a,B,c,lambda,dtapprox,3,deltaP,Wj,tapprox(i));
        end
        XDerFree(:,i) = XDerFree(:,i)/numruns;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%

    %Plot the Paths
    figure('Position',[10 10 900 600])
    plot3(t,X(1,:),X(2,:)); label11='Exact';
    hold on
    plot3(tapprox,XEuler(1,:),XEuler(2,:)); label12='Euler';
    plot3(tapprox,XTaylor(1,:),XTaylor(2,:)); label13='Taylor';
    plot3(tapprox,XDerFree(1,:),XDerFree(2,:)); label14='DF';
    legend(label11,label12,label13,label14)
    xlabel('Time (t)','FontSize',16)
    ylabel('Position ($X_t^{(1)}$ or $\tilde{X}_t^{(1)})$','Interpreter','latex','FontSize',16)
    zlabel('Position ($X_t^{(2)}$ or $\tilde{X}_t^{(2)})$','Interpreter','latex','FontSize',16)
    title('A Simulated Merton Model','FontSize',16)
    hold off

    %Plot the Errors
    figure('Position',[10 10 900 600])
    errorEuler = sum((X(:,1:floor(dtapprox/dt):end)-XEuler).^2)/2;
    errorTaylor = sum((X(:,1:floor(dtapprox/dt):end)-XTaylor).^2)/2;
    errorDerFree = sum((X(:,1:floor(dtapprox/dt):end)-XDerFree).^2)/2;
    plot(tapprox,errorEuler); label21='Euler';
    hold on
    plot(tapprox,errorTaylor); label22='Taylor';
    plot(tapprox,errorDerFree); label23='DF';
    legend(label21,label22,label23)
    xlabel('Time (t)','FontSize',16)
    ylabel('$\frac{1}{2}\left[(X_t^{(1)}-\tilde{X}_t^{(1)})^2+(X_t^{(2)}-\tilde{X}_t^{(2)})^2\right]$','Interpreter','latex','FontSize',16)
    title('Mean Squared Error for Each Time Step','FontSize',16)
    hold off
end

function [deltaP,Wj]=getDrivers(W,CP,dtapprox,dt,iter)
%%GETDRIVERS A helper function to automate retrieval of the correct
%            interval of values from the driving processes used for the
%            generation of the true path.
%
%INPUTS: W An mXlength(t) vector containing realized Wiener processes in
%          each row.
%       CP A 1Xlength(t) row vector containing a realized Poisson process.
% dtapprox The coarse approximation time step (assumed to be constant).
%       dt The fine true path time step (assumed to be constant).
%     iter The current iteration index (assumed indices are positive and
%          begin at 1).
%
%OUTPUTS: deltaP The number of jumps which occurred in CP during the time
%                interval.
%             Wj An mXdeltaP+2 vector which contains the values from W
%                which most closely occured at the same time as each jump
%                counted in deltaP. This assumes the distribution of jumps
%                over the time interval are uniformly random.
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.

    if(iter>1)
        currentIdx=floor(dtapprox/dt)*(iter-1)+1;
        prevIdx=floor(dtapprox/dt)*(iter-2)+1;
        deltaP=CP((iter-1)*floor(dtapprox/dt)+1)-CP((iter-2)*floor(dtapprox/dt)+1);
        if(deltaP~=0)
            jumpIdx=floor((prevIdx+rand([deltaP,1])*dtapprox/dt));
            Wj=W(:,[prevIdx;jumpIdx;currentIdx]);
        else
            Wj=W(:,[prevIdx;currentIdx]);
        end
    else
        error('Iter must be >1.')
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
