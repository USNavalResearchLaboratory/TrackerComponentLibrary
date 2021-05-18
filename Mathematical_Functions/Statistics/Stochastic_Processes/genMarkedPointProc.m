function [X,proc]=genMarkedPointProc(lambda,dt,ts,tf,X0,mark,eventDist,timeDist,staycell)
%%GENMARKEDPOINTPROC A function to simulate a marked point process with
%                    rate parameter lambda, initial value X0, and given 
%                    time parameters. Using the defaults for the last 3 
%                    parameters will generate a compound Poisson process
%                    provided the mark distribution is properly defined
%                    (see note below). The null mark is assumed to be zero.
%                    See [1] for more details on general marked point
%                    processes.
%
%INPUTS: lambda The unit rate parameter for the default Poisson
%               distribution. This may be arbitrarily defined if eventDist
%               is given.
%            dt Time step for event draws.
%            ts Start time for process.
%            tf End time for process.
%            X0 The initial value for the process. Defaults to 0 if
%               omitted.
%          mark An optional function handle f(x,t,j,proc) to be used for
%               sampling mark values. This must accept three scalar inputs
%               and a structure input:
%               x The process value at the beginning of a dt time step. 
%               t The elapsed time at the beginning of the dt time step,
%                 starting from ts. 
%               j The index corresponding to the ordinal number of the 
%                 jump since the beginning of the current time step.
%               proc A structure containing information about the process as
%                 currently constructed. See output information below for
%                 more details.
%               If mark is not given or is empty, it defaults to the unity
%               function for all potential values (jumps will have mark 1 
%               regardless of the variables).
%     eventDist An optional function handle for the distribution of events.
%               This must output a scalar value in the interval [0,inf). If
%               not given, it defaults to a Poisson distribution with rate
%               parameter lambda*dt. This function handle must take as
%               input:
%               x The process value at the beginning of a dt time step. 
%               t The elapsed time at the beginning of the dt time step,
%                 starting from ts.
%               i The index corresponding to the ordinal number of the 
%                 time step since the beginning of the process at time 
%                 ts.
%               proc A structure containing information about the process as
%                 currently constructed. See output information below for
%                 more details.
%     timeDist An optional function handle for the distribution of
%              times during a step. This must output a vector of values in
%              the interval (0,1). If not given, it defaults to a uniform
%              random distribution on (0,1). This function handle must take
%              as input:
%              x The process value at the beginning of a dt time step. 
%              t The elapsed time at the beginning of the dt time step,
%                 starting from ts.
%              i The index corresponding to the ordinal number of the 
%                 time step since the beginning of the process at time 
%                 ts.
%              proc A structure containing information about the process as
%                 currently constructed. See output information below for
%                 more details.
%     staycell A logical value. If this is equivalent to true, then X is 
%              returned as a cell array. Otherwise, X is returned as a 
%              vector. Defaults to false.
%
%OUTPUTS: X Defaults to a 1Xnsteps+1 vector of the process values after
%           each step. X(1) = X0. If the staycell variable evaluates to
%           true, then X is returned as a 1Xnsteps+1 cell array. 
%      proc A structure containing information which can be used to
%           reconstruct the exact process X. Members are:
%           nEvents A 1Xnsteps+1 vector of the number of jumps which
%                   occured during each step prior to the index. So 
%                   nEvents(1) = 0 and nEvents(i) will be the number of 
%                   jumps which occured during the (i-1) time step 
%                   corresponding to the change X{i}-X{i-1}.
%               mks A 1Xnsteps cell array where each cell contains a row
%                   vector with the mark values generated for a step.
%              times A 1X((tf-ts)/dt) cell array with each cell 
%                   containing a row vector containing the jump times for
%                   the process at the corresponding time step.
%
%Note: If the mark distribution is defined such that all marks are
%independently and identically distributed and the events are Poisson
%distributed, then the generated process will be a compound Poisson
%process. For more on compound Poisson processes, see pages 10-12 of [2].
%
%WARNING: No check is performed to ensure that eventDist and timeDist
%output reasonable values. This function could technically work with marks
%which are not numbers, provided some addition operation is defined for
%them. However, read  the documentation above to ensure your inputs will
%generate a proper process as intended.
%
%EXAMPLE 1: Generating several realizations of a time dependent marked
%           Poisson process and comparing it with the expected path.
% %Set path parameters
% ts = 40; % Start time
% tf = 50; % End time
% dt = 1; % Time step
% numpaths = 100; % Number of simulated paths
% t = ts:dt:tf;
% lambda1 = 5; % Poisson rate parameter for unit time
% lambda2 = 6; % Exponential rate parameter
% X0 = 100; % Initial value for the process
% 
% close all
% figure
% hold on
% 
% pathends = zeros([1,numpaths]);
% %Generate paths
% for i = 1:numpaths
%     [cp,proc] = genMarkedPointProc(lambda1,dt,ts,tf,X0,@(x,t,j,proc)(tf-t)*ExponentialD.rand(1,lambda2));
%     plot(t,cp)
%     pathends(i) = cp(end);
% end
% avgend = mean(pathends);
% stdend = std(pathends);
% 
% %Compute the expected path
% expected = X0+[0,cumsum(lambda1*dt*(tf-(t(1:end-1)+t(2:end))/2)/lambda2)];
% expectedstd = sqrt(cumsum(lambda1*dt*((tf-(t(1:end-1)+t(2:end))/2)).^2*2/(lambda2)^2));
% plot(t,expected,'--k','LineWidth',5)
% xlabel('time (t)')
% 
% %Compare simulated path statistics to the expected value
% fprintf("The average path value at t = %0.0f is : %0.5f\n",tf,avgend)
% fprintf("The path standard deviation at t = %0.0f is : %0.5f\n",tf,stdend)
% fprintf("The expected path value at t = %0.0f is : %0.5f\n",tf,expected(end))
% fprintf("The expected standard deviation at t = %0.0f is : %0.5f\n",tf,expectedstd(end))
%
%REFERENCES:
%[1] Stover, Christopher. "Marked Point Process." From MathWorld--A 
%    Wolfram Web Resource, created by Eric W. Weisstein. 
%    http://mathworld.wolfram.com/MarkedPointProcess.html
%[2] Platen, Eckhard, and Nicola Bruti-Liberati. Numerical solution of 
%    stochastic differential equations with jumps in finance. Vol. 64. 
%    Springer Science & Business Media, 2010. 
%
%July 2019 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(~exist('X0','var')||isempty(X0))
    X0 = 0;
end
if(~exist('mark','var')||isempty(mark))
    mark = @(x,t,j,proc) 1;
end
if(~exist('eventDist','var')||isempty(eventDist))
    eventDist = @(x,t,i,proc) PoissonD.rand(1,lambda*dt);
end
if(~exist('timeDist','var')||isempty(timeDist))
    timeDist = @(x,t,j,proc) rand([1,proc.nEvents(j)]);
end
if(~exist('staycell','var'))
    staycell = false;
elseif(~islogical(staycell))
    error('The value in staycell is not a logical value.')
end

nsteps = floor((tf-ts)/dt);
X = cell(1,nsteps+1);
X{1} = X0;

proc.nEvents = zeros([1,length(X)]);
proc.mks = cell([1,length(X)]);
proc.times = cell([1,length(X)]);

for i = 1:nsteps
    t = (i-1)*dt+ts;
    proc.nEvents(i) = eventDist(X{i},t,i,proc);
    stepmks = cell([1,proc.nEvents(i)]);
    steptimes = sort(timeDist(X{i},t,i,proc)*dt+t);
    proc.times{i} = steptimes;
    
    stepIntensity = 0;
    for j = 1:proc.nEvents(i)
        stepmks{j} = mark(X{i},steptimes(j),j,proc);
        proc.mks{i}{j} = stepmks{j};
        stepIntensity = stepIntensity+stepmks{j};
    end
    X{i+1} = X{i}+stepIntensity;
end

if(~staycell)
    X = cell2mat(X);
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
