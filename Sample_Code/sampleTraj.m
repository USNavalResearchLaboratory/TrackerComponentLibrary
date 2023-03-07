function [traj,times,mT,mTi] = sampleTraj(dt,choice)
%%SAMPLETRAJ A function for generating several example trajectories noted
%            in the tracking literature. Trajectories are returned
%            as a table with columns corresponding to each state variable
%            and rows corresponding to each time instance.
%
%INPUT:
% dt: The scalar time step to use for sampling the chosen trajectory.
% choice: A nonnegative integer indicating which trajectory to generate.
%         There are five test trajectories, all for 2D motion. POssible
%         values are:
%         0 Scenario 1 in Section 10.8.3 of [1]. This is non-maneuvering
%           flight, a 2G turn, and then non-maneuvering flight.
%         1 Scenario 2 in Section 10.8.3 of [1]. This is the same as 0 but
%           with different initial conditions and a 2.5G turn.
%         2 This is Scenario 3 in Section 10.8.3 of [1]. This is a sequence
%           of four 180 degree turns with magnitudes 1G, 1.5G, 3G, and
%           2.5G.
%         3 This is similar to Scenario 3 in Section 10.8.3 of [1], except
%           the initial state differs and the duration of the turns differs
%           so that the target performs a full loop. This is a sequence
%           of two turns with magnitudes 2.5G and 1.5G.
%         4 This is Scenario 3 in Section 10.8.3 of [1]. This is a sequence
%           of two turns with magnitudes 2.5G and 1.5G.
%
%OUTPUT:
% traj: A table with columns corresponding to each state variable and
%       rows corresponding to each time instance. The table2array function
%       can be used to change this output into a matrix.
% mT: A vector of the times at which changes in motion model occur.
% mTi: A vector of the time indicies at which changes in motion model
%      occur.
%
%EXAMPLE: Plots every trajectory on one plot.
% figure(1)
% clf;hold on;
% t = cell(1,4);
% for i = 0:4
%     t{i+1} = sampleTraj(1,i);
%     tarr = table2array(t{i+1});
%     plot(tarr(:,1),tarr(:,2),'Linewidth',2)
% end
% hold off
%
%REFERENCES:
%[1] X.-R. Li, Y. Bar-Shalom, and W. D. Blair, "Engineer's guide to
%    variable-structure multiple-model estimation for tracking,"
%    Multitarget-multisensor tracking: Applications and advances.,
%    vol. 3, pp. 499-567, 2000.
%
%June 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if ~isscalar(dt) || dt<=0
    error('Invalid time step was given.')
end

switch choice
    case 0
        % Times
        T = 150;
        TM1 = 61;
        TM2 = 106;
        mT = [TM1,TM2];

        % Time indices
        mTi = floor(TM1/dt);
        iTM2 = floor(TM2/dt);
        mTi = [mTi,iTM2];

        % Motion models
        FCV = FPolyKal(dt,4,1);
        initTR = deg2rad(3.74);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR]);

        % Generating trajectory
        times = 0:dt:T;
        initPos = [30e3;30e3];
        initVel = [-172;-246];
        tar = zeros(5,length(times));
        tar(:,1) = [initPos;initVel;0];
        idx = 1;
        for t = dt:dt:TM1
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR;
        for t = TM1+dt:dt:TM2
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = TM2+dt:dt:T
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        featureNames = {'xpos','ypos','xvel','yvel','turnrate'};
        traj = table(tar(1,:)',tar(2,:)',tar(3,:)',tar(4,:)',tar(5,:)',...
                     'VariableNames',featureNames,'RowNames',string(times),...
                     'DimensionNames',{'times','state variables'});
    case 1
        % Times
        T = 150;
        TM1 = 61;
        TM2 = 96;
        mT = [TM1,TM2];

        % Time indices
        mTi = floor(TM1/dt);
        iTM2 = floor(TM2/dt);
        mTi = [mTi,iTM2];

        % Motion models
        FCV = FPolyKal(dt,4,1);
        initTR = deg2rad(4.68);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR]);

        % Generating trajectory
        times = 0:dt:T;
        initPos = [20e3;50e3];
        initVel = [-172;-246];
        tar = zeros(5,length(times));
        tar(:,1) = [initPos;initVel;0];
        idx = 1;
        for t = dt:dt:TM1
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR;
        for t = TM1+dt:dt:TM2
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = TM2+dt:dt:T
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        featureNames = {'xpos','ypos','xvel','yvel','turnrate'};
        traj = table(tar(1,:)',tar(2,:)',tar(3,:)',tar(4,:)',tar(5,:)',...
                     'VariableNames',featureNames,'RowNames',string(times),...
                     'DimensionNames',{'times','state variables'});
    case 2
        % Times
        mT = cumsum([50,96,30,64,30,32,30,38,30]);

        % Time indices
        for i = 1:length(mT)
            mTi = floor(mT(i)/dt);
        end

        % Motion models
        FCV = FPolyKal(dt,4,1);
        initTR = [deg2rad(1.87),...
                  deg2rad(-2.8),...
                  deg2rad(5.6),...
                  deg2rad(-4.68)];

        % Generating trajectory
        times = 0:dt:mT(end);
        initPos = [60e3;40e3];
        speed = 300;
        initVel = [-speed/sqrt(2);speed/sqrt(2)];
        tar = zeros(5,length(times));
        tar(:,1) = [initPos;initVel;0];
        idx = 1;
        for t = dt:dt:mT(1)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(1);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(1)]);
        for t = mT(1)+dt:dt:mT(2)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(2)+dt:dt:mT(3)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(2);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(2)]);
        for t = mT(3)+dt:dt:mT(4)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(4)+dt:dt:mT(5)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(3);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(3)]);
        for t = mT(5)+dt:dt:mT(6)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(6)+dt:dt:mT(7)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(4);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(4)]);
        for t = mT(7)+dt:dt:mT(8)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(8)+dt:dt:mT(end)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        
        featureNames = {'xpos','ypos','xvel','yvel','turnrate'};
        traj = table(tar(1,:)',tar(2,:)',tar(3,:)',tar(4,:)',tar(5,:)',...
                     'VariableNames',featureNames,'RowNames',string(times),...
                     'DimensionNames',{'times','state variables'});
    case 3
        % Times
        mT = cumsum([50,96,30,64,30]);

        % Time indices
        for i = 1:length(mT)
            mTi = floor(mT(i)/dt);
        end

        % Motion models
        FCV = FPolyKal(dt,4,1);
        initTR = [deg2rad(4.68),...
                  deg2rad(-2.8)];

        % Generating trajectory
        times = 0:dt:mT(end);
        initPos = [60e3;40e3];
        speed = 300;
        initVel = [-speed/sqrt(2);speed/sqrt(2)];
        tar = zeros(5,length(times));
        tar(:,1) = [initPos;initVel;0];
        idx = 1;
        for t = dt:dt:mT(1)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(1);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(1)]);
        for t = mT(1)+dt:dt:mT(2)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(2)+dt:dt:mT(3)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(2);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(2)]);
        for t = mT(3)+dt:dt:mT(4)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(4)+dt:dt:mT(5)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        
        featureNames = {'xpos','ypos','xvel','yvel','turnrate'};
        traj = table(tar(1,:)',tar(2,:)',tar(3,:)',tar(4,:)',tar(5,:)',...
                     'VariableNames',featureNames,'RowNames',string(times),...
                     'DimensionNames',{'times','state variables'});
    case 4
        % Times
        mT = cumsum([50,28,82,50]);

        % Time indices
        for i = 1:length(mT)
            mTi = floor(mT(i)/dt);
        end

        % Motion models
        FCV = FPolyKal(dt,4,1);
        initTR = [deg2rad(4.68),...
                  deg2rad(-2.8)];

        % Generating trajectory
        times = 0:dt:mT(end);
        initPos = [10e3;50e3];
        speed = 300;
        initVel = [speed/sqrt(2);-speed/sqrt(2)];
        tar = zeros(5,length(times));
        tar(:,1) = [initPos;initVel;0];
        idx = 1;
        for t = dt:dt:mT(1)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        tar(5,idx) = initTR(1);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(1)]);
        for t = mT(1)+dt:dt:mT(2)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = initTR(2);
        FCT = FCoordTurn2D(dt,[zeros(4,1);initTR(2)]);
        for t = mT(2)+dt:dt:mT(3)
            idx = idx+1;
            tar(:,idx) = FCT*tar(:,idx-1);
        end
        tar(5,idx) = 0;
        for t = mT(3)+dt:dt:mT(4)
            idx = idx+1;
            tar(1:4,idx) = FCV*tar(1:4,idx-1);
        end
        
        featureNames = {'xpos','ypos','xvel','yvel','turnrate'};
        traj = table(tar(1,:)',tar(2,:)',tar(3,:)',tar(4,:)',tar(5,:)',...
                     'VariableNames',featureNames,'RowNames',string(times),...
                     'DimensionNames',{'times','state variables'});
    otherwise
        error('Invalid choice value was given.')
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