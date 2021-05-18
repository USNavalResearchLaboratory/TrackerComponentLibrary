function demo2DDataAssociation(algSel1,algSel2,param3)
%%DEMO2DDATAASSOCIATION Demonstrate how variants of the JPDAF and nearest
%              neighbor filters can be used to track three maneuvering
%              targets in the presence of false alarms and possible missed
%              detections. The JPDAF* is used by default. Clustering is
%              performed so as to make the exact computation of the
%              target-measurement association probabilities computationally
%              practicable.
%
%INPUTS: algSel1,algSel2,param3 These are the same inputs as in the
%                 function singleScanUpdate and they can be used to select
%                 the specific method for handling target-measurement
%                 association uncertainty. If these parameters are omitted,
%                 then the JPDAF* is used (algSel1=5 and the others inputs
%                 do not matter).
%
%This is a simple two-dimensional (x,y) simulation scenario involving three
%maneuvering aircraft that cross, come within close range of each other and
%then separate. The scenario involves the trajectory used in the air
%traffic control (ATC) scenario discussed in Chapter 11.7.4 of [1] and a
%shifted, mirrored version of it as well as a straight trajectory that goes
%between the two. As measurement filtering is not the focus of this
%problem, a converted-measurement Kalman filter with a second-order
%Gauss-Markov (Singer) dynamic model is used with polar-Cartesian converted
%measurements. The dynamic model consists of elements of position, velocity
%and acceleration.
%
%Such a simulation scenario might be representative of what one gets from a
%rotating radar when all contacts from a single rotation are collected into
%a scan rather than performing updates for every single dwell. In such a
%scenario, not all of the measurements will be taken at the same time,
%meaning that different targets will be predicted to different times based
%on where they are. To simplify the presentation here, all measurements are
%taken to be at the same time.
%
%This function plots the true tracks in black, the estimated tracks in red,
%and the measurents are plotted as circles. A green X marks the location of
%the sensor.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%%%%%%
%%%SET UP THE SCENARIO AND SIMULATION PARAMETERS
%%%%%%

%This selects which type of JPDAF algorithm to run. Here, we select the
%JPDAF* by default.
if(nargin<1||isempty(algSel1))
    algSel1=5;
end
if(nargin<2||isempty(algSel2))
    algSel2=0;
end
if(nargin<3||isempty(param3))
    param3=[];
end

numTargets=3;
xDimEst=6;
xDimTrue=4;
zDim=2;
[xTrue1,T,numSamples]=genManeuveringTrack();

xTrue2=xTrue1;
xTrue2(1,:)=-xTrue2(1,:)+6.5e3;
xInit3=[3124.50645843013-250;
        28924.50645843013;
        0;
        -120];
xTrue3=genStraightTrack(xInit3,T,numSamples);
xTrue=cat(2,reshape(xTrue1,xDimTrue,1,numSamples),reshape(xTrue2,xDimTrue,1,numSamples),reshape(xTrue3,xDimTrue,1,numSamples));

%Parameters for the dynamic model. We are using a second-order Gauss-Markov
%model (Singer dynamic model).
tau=20;%20 seconds; the assumed maneuver decorrelation time.
%Assumed maximum acceleration change in one second.
maxJerk=9.8/3;
%Rule-of-thumb process noise suggestion.
q0=processNoiseSuggest('PolyKal-ROT',maxJerk,T);
%Process noise covariance matrix.
Q=QGaussMarkov(T,xDimEst,q0,tau,2);
SQ=chol(Q,'lower');
F=FGaussMarkov(T,xDimEst,tau,2);%State transition matrix

%The viewing region range. This is important for dealing with the clutter
%model. The clutter is generated in the coordinate system of the polar
%measurements.
minRange=1e3;
maxRange=40e3;

%Assumed clutter parameter in terms of false alarms per meter-radian
%(polar coordinates). False alarms are generated in the measurement
%domain of the radar. The tracker can still work if the value for lambda
%does not perfectly match the value simulated/ the real value. However, one
%should note that lambda in the tracker cannot be arbitrarily large lest
%the tracker eventually ignore all measurements (false alarms become always
%more likely than true tracks).
lambda=1e-4;

%The bounds of the viewing region in range and angle.
bounds=[minRange,-pi;
        maxRange,pi];

%Gate probability
PG=0.9997;
%Threshhold for gating
gammaVal=ChiSquareD.invCDF(PG,2);

%Noise variances and square root covariance matrix.
%Assumed standard deviations of the measurement noise components.
sigmaR=10;
sigmaAz=0.2*(pi/180);
%Square root measurement covariance matrix; assume no correlation.
SR=diag([sigmaR,sigmaAz]);

sensorLoc=[-7e3;0];
PD=0.8;%Detection probability.

%Cubature points for measurement conversion.
[xi,w]=fifthOrderCubPoints(2);

%Allocate space to hold the track estimates.
xEst=zeros(xDimEst,numTargets,numSamples);
%Allocate space to hold the lower-triangular square roots of the
%covariances of the track estimates.
SEst=zeros(xDimEst,xDimEst,numTargets,numSamples);

%%%%%%
%%%INITALIZE THE TRACKS USING TWO POINTS
%%%%%%

%First, we initialize the tracks using hand-picked measurements. All of the
%non-integrated JPDAF variants do not perform initiation and termination.
%We get an initial position and velocity. The acceleration is just set to
%zero with a relatively large covariance to indicate that it is mostly
%unknown.
%The standard deviation of acceleration for initialization is 4G in each
%dimension.
initAccelStdDev=9.8;
for curTar=1:numTargets
    %Generate measurements for the first two steps.
    z=zeros(zDim,2);
    R=zeros(zDim,zDim,2);
    
    %Get the polar location and add noise.
    z1=Cart2Pol(xTrue(1:2,curTar,1)-sensorLoc)+SR*randn(2,1);
    %We are using a Cartesian-converted filter.
    [z(:,1),R(:,:,1)]=pol2CartCubature(z1,SR,0,true,sensorLoc,sensorLoc,[],xi,w);

    z2=Cart2Pol(xTrue(1:2,curTar,2)-sensorLoc)+SR*randn(2,1);
    [z(:,2),R(:,:,2)]=pol2CartCubature(z2,SR,0,true,sensorLoc,sensorLoc,[],xi,w);
    
    [xInit,PInit]=twoPointDiffInit(T,z,R,q0);
    xEst(1:4,curTar,2)=xInit;
    SEst(1:4,1:4,curTar,2)=chol(PInit,'lower');
    SEst(5,5,curTar,2)=initAccelStdDev;
    SEst(6,6,curTar,2)=initAccelStdDev;
end

%%%%%%
%%%GENERATE MEASUREMENT AND PERFORM FILTERING
%%%%%%

figure(3)
clf
hold on
for curStep=3:numSamples
    %GENERATE MEASUREMENTS
    zTrueCart=xTrue(1:2,:,curStep);
    zPol=genObs(PD,lambda,bounds,zTrueCart,sensorLoc,SR);

    %We are using a Cartesian-converted filter, so the measurements must be
    %converted.
    [zCart,RCart]=pol2CartCubature(zPol,SR,0,true,sensorLoc,sensorLoc,[],xi,w);
    numMeas=size(zCart,2);
    
    %We need determinants of the measuremnt Jacobian in order to transform
    %the clutter density at the measurements for assigning likelihoods.
    measJacobDet=zeros(numMeas,1);
    for k=1:numMeas
        measJacobDet(k)=det(calcPolarConvJacob(zPol(:,k),0,true,sensorLoc));
    end

    %Transform the covariance matrices to square-root covariance matrices
    SRCart=zeros(zDim,zDim,numMeas);
    for curMeas=1:numMeas
        SRCart(:,:,curMeas)=chol(RCart(:,:,curMeas),'lower');
    end
    
    %Display the measurements
    scatter(zCart(1,:),zCart(2,:),'ob')
    
    %%PREDICT TARGETS FORWARD
    for curTar=1:numTargets
        [xEst(:,curTar,curStep),SEst(:,:,curTar,curStep)]=sqrtDiscKalPred(xEst(:,curTar,curStep-1),SEst(:,:,curTar,curStep-1),F,SQ);
    end

    %%GET HEURISTIC MAXIMUM GATE SIZE
    %The function makeStandardCartOnlyLRMatHyps can adjust for the gate
    %size being finite. The maximum gate size to be considered for each
    %target could be approximated using a measurement covariance at the
    %predicted target location. However, here, we simply use a maximum
    %value over all of the measurements in each dimension as a heuristic.
    if(numMeas>0)
        SRMax=zeros(2,2);
        SRMax(1,1)=max(SRCart(1,1,:));
        SRMax(2,2)=max(SRCart(2,2,:));
    else%If there are no measurements, we assume that nothing useful could
        %have been outside the gates, so no covariance inflation is
        %performed. In practical scenarios, this may or may not be a good
        %approximation. However, with PG generally near 1, the covariance
        %inflation often does not matter.
        SRMax=[];
    end
    %In practice, however, this won't really make a difference if the gate
    %probability PG is very close to 1 and PD is not close to 1 (one could
    %have set SRMax=[], not done the correction, and gotten usually the
    %same result) 
    
    %%COMPUTE LIKELIHOODS AND GATE
    [A,xHyp,PHyp,GateMat]=makeStandardCartOnlyLRMatHyps(xEst(:,:,curStep),SEst(:,:,:,curStep),zCart,SRCart,SRMax,PD,lambda,[],gammaVal,measJacobDet);

    %PERFORM CLUSTERING AND UPDATE EACH CLUSTER INDEPENDENTLY
    %Next, we can to break the targets and measurements up into clusters
    %for the updating. We do this, because the computation of the
    %target-measurement association probabilities in the JPDAF is very slow
    %if one does not perform clustering or use an approximation. In
    %practice, if cluster sizes become too large, one would have to fall
    %back to an approximation. A DisjointSetM data structure is used to
    %perform the clustering. The result is then extracted into a ClusterSet
    %class, which allows for indexation by cluster number and by the index
    %of the thing in the cluster. The DisjointSetM class specifically is
    %used to obtain measurements in the clusters as well.
    theSet=DisjointSetM(numTargets,numMeas);
    theSet.unionFromBinMat(GateMat);
    %This method could also return a list of indices of measurements that
    %do not gate with anything.
    [targetCluster,measCluster]=theSet.createClusterSet();

    numClusters=targetCluster.numClusters;
    
    %We run the filter separately for each cluster. For clusters with no
    %measurements (individual tracks that did not gate with anything), the
    %prediction is the same as the update.
    for curClust=1:numClusters
        tarIdx=targetCluster(curClust,:);
        measIdx=measCluster(curClust,:);
 
        %If something gated, then update the estimates. Otherwise, the
        %prediction remains the estimate.
        if(~isempty(measIdx))
            %Select the measurement update hypotheses for the given
            %measurements and targets as well as the missed detection
            %hypothesis for each target.
            xHypClust=cat(3,xHyp(:,tarIdx,measIdx),xHyp(:,tarIdx,end));
            PHypClust=cat(4,PHyp(:,:,tarIdx,measIdx),PHyp(:,:,tarIdx,end));
            numTarCur=length(tarIdx);
            AClust=[A(tarIdx,measIdx),A(tarIdx,numMeas+tarIdx)];
            
            [xPost,PPost]=singleScanUpdate(xHypClust,PHypClust,AClust,algSel1,algSel2,param3);
            
            xEst(:,tarIdx,curStep)=xPost;
            for curTar=1:numTarCur
                SEst(:,:,tarIdx(curTar),curStep)=chol(PPost(:,:,curTar),'lower');
            end
        end
    end

    for curTar=1:3
        x=[xEst(1,curTar,curStep-1),xEst(1,curTar,curStep)];
        y=[xEst(2,curTar,curStep-1),xEst(2,curTar,curStep)];
        plot(x(:),y(:),'-r','linewidth',2)
    end
end

%Plot the true tracks and the estimated tracks. 
%Plot the estimated tracks
for curTrack=1:3
    x=xEst(1,curTrack,2:end);
    y=xEst(2,curTrack,2:end);
    plot(x(:),y(:),'-r','linewidth',4)
end

%Plot the true tracks
for curTrack=1:3
    x=xTrue(1,curTrack,:);
    y=xTrue(2,curTrack,:);
    plot(x(:),y(:),'--k','linewidth',2)
end

%Mark the location of the sensor
scatter(sensorLoc(1),sensorLoc(2),250,'xg','linewidth',4)

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

function xTrue=genStraightTrack(xInit,T,numSamples)

F=FPolyKal(T,4,1);

xTrue=zeros(4,numSamples);
xTrue(:,1)=xInit;

for curStep=2:numSamples
    xTrue(:,curStep)=F*xTrue(:,curStep-1);
end

end

function z=genObs(PD,lambda,bounds,zTrueCart,sensorLoc,SR)
%Generate observations in polar coordinates with false alarms.

numTar=size(zTrueCart,2);
tarsDetected=rand(numTar,1)<PD;
lambdaV=lambda*prod(bounds(2,:)-bounds(1,:));

numClutter=PoissonD.rand(1,lambdaV);

numObs=sum(tarsDetected)+numClutter;
z=zeros(2,numObs);

%First we will generate the target observations.
curObs=1;
for curTar=1:numTar
    %Everything is observed at the first step so that we can easily
    %do two point differencing later.
    if(tarsDetected(curTar))
        curLoc=zTrueCart(:,curTar);%Get the Cartesian Location

        %Get the polar location and add noise.
        z(:,curObs)=Cart2Pol(curLoc-sensorLoc)+SR*randn(2,1);
        curObs=curObs+1;
    end
end
        
%Now we will generate the false alarms. False alarms are generated in polar
%coordinates.
rClutBounds=bounds(:,1);
azBounds=bounds(:,2);
for curFalse=1:numClutter
    r=UniformD.rand(1,rClutBounds);
    az=UniformD.rand(1,azBounds);

    z(:,curObs)=[r;az];
    curObs=curObs+1;
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
