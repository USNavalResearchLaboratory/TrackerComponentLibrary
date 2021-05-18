function demo2DIntegratedDataAssociation()
%%DEMO2DINTEGRATEDDATAASSOCIATION Demonstrate how the Global Nearest
%                       Neighbor (GNN) Joint Integrated Probabilistic Data
%                       Association Filter (JIPDAF) can be used to track
%                       targets performing automatic track initiation and
%                       termination. This is a simple scenario where no
%                       gating/ clustering is performed, so approximate
%                       data association probabilities must be used for the
%                       problem to be computationally practicable.
%
%This is a simple two-dimensional (x,y) simulation scenario involving two
%maneuvering aircraft that cross, come within close range of each other and
%then separate. The scenario involves the trajectory used in the air
%traffic control (ATC) scenario discussed in Chapter 11.7.4 of [1] and a
%shifted, mirrored version of it. As measurement filtering is not the focus
%of this problem, a converted-measurement Kalman filter with a first-order
%Gauss-Markov (integrated Ornstein-Uhlenbeck) dynamic model is used with
%polar-Cartesian converted measurements.
%
%Such a simulation scenario might be representative of what one gets from a
%rotating radar when all contacts from a single rotation are collected into
%a scan rather than performing updates for every single dwell. In such a
%scenario, not all of the measurements will be taken at the same time,
%meaning that different targets will be predicted to different times based
%on where they are. To simplify the presentation here, all measurements are
%taken to be at the same time.
%
%As approximate target-measurement association probabilities are used here,
%it was observed that the algorithm is sufficiently fast without performing
%gating and then jointly processing groups of targets that gate with common
%measurements together. That is, all targets are processed jointly in this
%implementation, which also simplifies it. If gating and grouping targets
%are neded to reduce computational complexity in other applications (such
%as if one wishes to use exact target-measurement association
%probabilities), then one can use the DisjointSetM class to cluster the
%targets and measurements into groups. Such clustering, without discussion
%of data structures to perform it, is mentioned in the original paper on
%the JIPDA in [2].
%
%The scenario is run with false alarms and missed detections. To illustrate
%timely termination of tracks, an extra 50 time-steps of simulation are
%performed after the trajectories end to demonstrate that they are
%terminated in a timely manner.
%
%This function plots all of the tracks found by the tracking algorithm in
%green. Tracks are only initially displayed once their existence
%probability exceeds 95% (the PDisp parameter). After that, they continue
%to be displayed until they are terminated. Tracks are terminated if their
%probability of existence is less than 0.01%.
%
%Note that performance could be improved if range-rate (Doppler)
%information were simulated and used. Specifically, even if not used in the
%measurement update step, range-rate information could improve the
%likelihoods used for target-measurement association, reducing the
%occurrence of false tracks.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%[2] D. Musicki and R. Evans, "Joint integrated probabilistic data
%    association: JIPDA," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 3, pp. 1093-1099, Jul. 2004.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('An example of a simple Kalman filter with Gauss-Markov model, data association, and track initiation/ termination in 2D') 

%%%SETUP THE MODEL
disp('Setting parameters and creating the simulated trajectories.') 
PD=0.8;%Detection probability --same for all targets.

%Assumed clutter parameter in terms of false alarms per meter-radian
%(polar coordinates). False alarms are generated in the measurement
%domain of the radar. The tracker can still work if the value for lambda
%does not perfectly match the value simulated/ the real value. However, one
%should note that lambda in the tracker cannot be arbitrarily large lest
%the tracker eventually ignore all measurements (false alarms become always
%more likely than true tracks).
lambda=2e-5;

%The AlgSel parameters are inputs to the singleScanUpdate function.
%Use the GNN-JIPDAF with approximate association probabilities.
algSel1=3;
algSel2=0;
param3=[];

%Assumed standard deviations of the measurement noise components.
sigmaR=10;
sigmaAz=0.1*(pi/180);
%Square root measurement covariance matrix; assume no correlation.
SR=diag([sigmaR,sigmaAz]);

%The viewing region range. This is important for dealing with the clutter
%model.
minRange=1e3;
maxRange=40e3;
%The angular region is two pi (all around).

%lambda times the "volume" in the receiver's polar coordinate system needed
%for the Poisson clutter model
lambdaV=lambda*(maxRange-minRange)*2*pi;

%The transition probability when predicting forward the target states that
%the track did not end during the transition. This is a design parameter. 
PStayAlive=0.98;

%This is the initial probability assigned that the detection might be a
%real, new track. This value is a design parameter.
PIsRealInit=0.1;

%The probability of existence below which a track is terminated.
PTerminate=1e-4;

xDim=4;
zDim=2;

%Potential new tracks are started using one-point differencing. Usually,
%this means that the velocity value is set to zero and one uses a
%standard deviation value (in each dimension) based loosely on a standard
%maximum velocity.
velInitStdDev=500;%(m/s)

%Cubature points for measurement conversion.
[xi,w]=fifthOrderCubPoints(zDim);
H=[1,0,0,0;
   0,1,0,0];%Cartesian measurement matrix

%Generate the two true targets. The time between samples is T; the total
%number of samples is in numSamples.
[xTrue1,T,numSamplesTrack]=genManeuveringTrack();
%We will go for 50 extra samples (with no tracks) to make sure that the
%trackers terminate the tracks in a timely manner.
numSamples=numSamplesTrack+50;

xTrue2=xTrue1;
xTrue2(1,:)=-xTrue2(1,:)+6.5e3;%x position component.
xTrue2(3,:)=-xTrue2(3,:);%Velocity component.
xTrue=cat(3,xTrue1,xTrue2);
numTargets=2;

%Parameters for the dynamic model. We are using a first-order Gauss-Markov
%model.
tau=20;%20 seconds; the assumed maneuver decorrelation time.
maxAccel=9.8*3;%Assume max 3G turn
%Rule-of-thumb process noise suggestion.
q=processNoiseSuggest('PolyKal-ROT',maxAccel,T);
Q=QGaussMarkov(T,xDim,q,tau,1);%Process noise covariance matrix.
SQ=chol(Q,'lower');

F=FGaussMarkov(T,xDim,tau,1);%State transition matrix

disp('Generating measurements.') 

%Generate measurements and false alarms for each scan.
zMeasCart=cell(numSamples,1);
SRMeasCart=cell(numSamples,1);
zMeasPolar=cell(numSamples,1);
zMeasJacobDet=cell(numSamples,1);
for curScan=1:numSamples
    %Determine the number of false alarms to generate.
    numFalse=PoissonD.rand(1,lambdaV);

    if(curScan<=numSamplesTrack)
        %Determine which, if any, targets should be detected.
        isDet=rand(numTargets,1)<PD;
    else
        %We are after the end; the tracks are now gone.
        isDet=[0;0];
    end

    %Allocate space for the detections.
    numMeas=numFalse+sum(isDet);
    zCur=zeros(2,numMeas);
    curDet=1;

    if(numMeas>0)
        %Generate the detection from the targets, if any.
        for curTar=1:numTargets
            if(isDet(curTar))
                zCur(:,curDet)=Cart2Pol(H*xTrue(:,curScan,curTar))+SR*randn(2,1);
                curDet=curDet+1;
            end
        end

        %Generate the false alarm detections, if any. 
        rClutBounds=[minRange;maxRange];
        azBounds=[-pi;pi];
        for curFalse=1:numFalse
            r=UniformD.rand(1,rClutBounds);
            az=UniformD.rand(1,azBounds);

            zCur(:,curDet)=[r;az];
            curDet=curDet+1;
        end
    
        zMeasPolar{curScan}=zCur;
        
        %We will now convert the measurements into Cartesian coordinates as
        %we are using a converted-measurement filter.
        [zMeasCart{curScan},RMeasCart]=pol2CartCubature(zCur,SR,0,true,[],[],[],xi,w);
        
        %Take the lower-triangular square root of the covariance matrices.
        measDetCur=zeros(numMeas,1);
        for curMeas=1:numMeas
            RMeasCart(:,:,curMeas)=chol(RMeasCart(:,:,curMeas),'lower');
            measDetCur(curMeas)=det(calcPolarConvJacob(zCur(:,curMeas),0,true));
        end
        SRMeasCart{curScan}=RMeasCart;
        zMeasJacobDet{curScan}=measDetCur;
    end
end

%Now for the tracker with integrated track initiation/ termination.
disp('Running the integrated tracking algorithm.') 

%We will save the value of each track at each time. Thus, a cell array
%holds the track states at each time. The first column in xStates is the
%actual values, the second column in xStates is an ID so that we can
%associate the tracks across time.
xStates=cell(numSamples,2);
%Lower-triangular square root covariance matrices.
SStates=cell(numSamples,1);
%Existence probabilities
rStates=cell(numSamples,1);
for curScan=1:numSamples
    zCur=zMeasCart{curScan};
    SRCur=SRMeasCart{curScan};
    numMeas=size(zCur,2);
    measJacobDets=zMeasJacobDet{curScan};
    
    %First, we create potential tracks for all of the observations.
    xNew=zeros(xDim,numMeas);
    SNew=zeros(xDim,xDim,numMeas);
    xNew(1:zDim,:)=zCur;
    SNew(1:zDim,1:zDim,:)=SRCur;
    %The uncertainty for the unknown velocity
    for curDim=(zDim+1):xDim
        SNew(curDim,curDim,:)=velInitStdDev;
    end
    %Initialization existence probabilities
    rNew=PIsRealInit*ones(numMeas,1);
    
    %Generate a UUID for each potential track so that we can associate
    %tracks over time to draw lines for display. The UUIDs are 36-character
    %strings. For simplicity, we are just using the hash values of the
    %UUIDs so that they can be easily compared with >,=,< for sorting.
    xIDNew=zeros(numMeas,1);
    for curNewTrack=1:numMeas
        [~,xIDNew(curNewTrack)]=genUUID();
    end
    
    %Next, if there are any existing states, we want to predict them to the
    %current time step and update them with the measurements
    if(curScan>1&&~isempty(xStates{curScan-1,1}))
        x=xStates{curScan-1,1};
        S=SStates{curScan-1};
        r=rStates{curScan-1};
        xID=xStates{curScan-1,2};
        numTargetsCur=size(x,2);
        
        %Predict the tracks to the current time.
        for curTar=1:numTargetsCur
            [x(:,curTar),S(:,:,curTar)]=sqrtDiscKalPred(x(:,curTar),S(:,:,curTar),F,SQ);
        end
        %Update the target existence probabilities with the Markov 
        %switching model.
        r=PStayAlive*r;

        %The inclusion of r takes into account the track existence
        %probabilities.
        [A,xHyp,PHyp]=makeStandardCartOnlyLRMatHyps(x,S,zCur,SRCur,[],PD,lambda,r,[],measJacobDets);
        
        [xUpdate,PUpdate,rUpdate,probNonTargetMeas]=singleScanUpdateWithExistence(xHyp,PHyp,PD,r,A,algSel1,algSel2,param3);
        
        %Determine which tracks to drop because their existence
        %probabilities are below the termination probability.
        sel=rUpdate>PTerminate;
        numTargetsCur=sum(sel);
        xUpdate=xUpdate(:,sel);
        PUpdate=PUpdate(:,:,sel);
        rUpdate=rUpdate(sel);
        xID=xID(sel);

        SUpdate=zeros(size(PUpdate));
        for curTar=1:numTargetsCur
            SUpdate(:,:,curTar)=chol(PUpdate(:,:,curTar),'lower');
        end
        
        rNew=probNonTargetMeas.*rNew;
        
        sel=rNew>PTerminate;
        xNew=xNew(:,sel);
        SNew=SNew(:,:,sel);
        rNew=rNew(sel);
        xIDNew=xIDNew(sel);
        
        xStates{curScan,1}=[xNew,xUpdate];
        xStates{curScan,2}=[xIDNew;xID];
        SStates{curScan}=cat(3,SNew,SUpdate);
        rStates{curScan}=[rNew;rUpdate(:)];
    else
        xStates{curScan,1}=xNew;
        xStates{curScan,2}=xIDNew;
        SStates{curScan}=SNew;
        rStates{curScan}=rNew;
    end
end

disp('Displaying contacts, true tracks, and tracks found.') 

figure(1)
clf
hold on
plot(xTrue1(1,:),xTrue1(2,:),'-b','linewidth',4)
plot(xTrue2(1,:),xTrue2(2,:),'-r','linewidth',4)

PDisp=0.98;

trackList=AVLTree();
for curScan=1:numSamples
    zCur=zMeasCart{curScan};
    if(~isempty(zCur))
        scatter(zCur(1,:),zCur(2,:),'ok');
    end
    
    stateCur=xStates{curScan,1};
    IDCur=xStates{curScan,2};
    rCur=rStates{curScan};
    
    numHyp=length(rCur);
    
    trackListNew=AVLTree();
    
    for curHyp=1:numHyp
        curID=IDCur(curHyp);
        prevTrackLoc=trackList.find(curID);
        
        %If it is an existing track, then connect the points
        if(~isempty(prevTrackLoc))
            prevTrackLoc=prevTrackLoc.value;
            
            curTrackLoc=stateCur(1:2,curHyp);
            plot([prevTrackLoc(1);curTrackLoc(1)],[prevTrackLoc(2);curTrackLoc(2)],'-g','linewidth',2);
            trackListNew.insert(KeyVal(curID,curTrackLoc));
        elseif(rCur(curHyp)>PDisp)
            %We only start the tentative track when it is above the
            %detection threshold for track existence.
            curTrackLoc=stateCur(1:2,curHyp);
            
            if(~isempty(prevTrackLoc))
                prevTrackLoc=prevTrackLoc.value;
                plot([prevTrackLoc(1);curTrackLoc(1)],[prevTrackLoc(2);curTrackLoc(2)],'-g','linewidth',2);
            else
                scatter(curTrackLoc(1),curTrackLoc(2),'.g')
            end
            
            trackListNew.insert(KeyVal(curID,curTrackLoc));
        end
    end
    
    trackList=trackListNew;
end

h1=xlabel('x meters');
h2=ylabel('y meters');
h3=title('True and Estimated Trajectories');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
axis square

end

function [xTrue,T,numSamples]=genManeuveringTrack()
T=5;
TTotal=125+90+125+30+125;
%The extra sample is the very first step in the trajectory.
numSamples=TTotal/T+1;
xTrue=zeros(4,numSamples);

%An indicator of the maneuvering mode 1 for not maneuvering, and 2 for
%maneuvering.
mode=zeros(numSamples,1);

%The plane first flies westward at 120 m/s for 125 seconds.
xInit=[25e3;10e3; -120; 0];
xTrue(:,1)=xInit;
mode(1)=1;
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

%Propagate the constant velocity model for 125 seconds.
baseStep=1;
for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
end
baseStep=baseStep+numStepSeg;

%Next, the plane performs a 1 degree per second coordinated turn for 90
%seconds.
turnRate=1*(pi/180);%Convert from degrees per second to radians per second.

numStepSeg=90/T;
for curStep=(baseStep+1):(baseStep+numStepSeg)
    %Get the discrete-time propagation matrix for the state. The
    %FCoordTurn2D function assumes that the turn rate is part of the state.
    %In this instance, we just add it to the state, and then remove the
    %extra elements from the returned matrix.
    F=FCoordTurn2D(T,[xTrue(:,curStep-1);turnRate]);
    F=F(1:4,1:4);
    
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=2;
end
baseStep=baseStep+numStepSeg;

%Next, go straight (South after the previous turn) for another 125 seconds.
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
end
baseStep=baseStep+numStepSeg;

%Make another turn. This time, it is at a rate of -three degrees per second
%for 30 seconds.
turnRate=-3*(pi/180);%Degrees per second to radians per second.

numStepSeg=30/T;%The number of steps in this segment.
for curStep=(baseStep+1):(baseStep+numStepSeg)
    F=FCoordTurn2D(T,[xTrue(:,curStep-1);turnRate]);
    F=F(1:4,1:4);
    
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=2;
end
baseStep=baseStep+numStepSeg;

%Finally, go straight again for another 125 seconds.
numStepSeg=125/T;%The number of steps in this segment.

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

for curStep=(baseStep+1):(baseStep+numStepSeg)
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
    mode(curStep)=1;
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
