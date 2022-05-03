function demoMultipleModelFiltering()
%%DEMOMULTIPLEMODELFILTERING This demonstrates how multiple model filters
%           can be used. The use of multiple dynamic models can improve the
%           accuracy of tracking targets than can maneuver. The scenario
%           considered is the simple 2D tracking of an aircraft. The
%           separate problems of data association, scheduling, track
%           initiation and termination are not addressed here.
%
%The example scenario is the air traffic control (ATC) scenario discussed
%in Chapter 11.7.4 of [1]. However, rather than using a fixed mode
%transition probability matrix, the use of the function
%getMarkovPTransProbMat to get the mode transition probability matrix given
%mean sojorn times and the time between samples is demonstrated. This
%allows the algorithms to be used with variable sampling rates, though in
%this instance, the sampling rate is held constant.
%
%Six scenarios are considered. One is a baseline Kalman filter with a
%polynomial dynamic model and a large process noise. The second is a Kalman
%filter with a first-order Gauss-Markov model (the integrated Ornstein-
%Uhlenbeck model), which, unlike the standard polynomial model, includes a
%correlation constant to try to better cover turns. The third is the
%reduced state filter. The fourth is the separated covariance filter. The
%other two scenarios utilize the interacting multiple model (IMM) filter.
%The first IMM scenario is as in [1], where two nearly constant velocity
%dynamic models with differing process noises are used in a standard Kalman
%filter (with the IMM). The second IMM scenario is where a nearly constant
%velocity model in a Kalman filter is used as well as a coordinated turn
%model in an extended Kalman filter (EKF). The functions multipleModelPred
%and multipleModelUpdate are used to handle the prediction and mixing in
%the IMM, though a change of the algorithm selection would allow one to
%just as easily use a generalized pseudo-Bayesian 2 filter. The IMM
%estimators can be more difficult to use as there are more parameters to
%tune.
%
%The measurements are generated in Cartesian space as measurement
%conversion and tracking using non-Cartesian measurements is not the focus
%of this example.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%
%June 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%This can be changed to try other multiple model algorithms. If GPB1 is
%selected, the multiple model filter using the coordinate dturn model is
%not used, because the GPB1 requires that all models have the same
%dimensionality.
AlgSel='GPB2';

%Interval between samples in seconds.
T=5;

%%STEP 1: Generate the simulation scenario.
disp('Generating the "true" track for the simulation scenario.')
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

%%STEP 2: Generate measurements.
disp('Generating the measurements for the simulation; A simple Cartesian model is used.')
%The measurements are just Cartesian with 100 meters of noise added per
%dimension, with no correlation between them.
H=[1,0,0,0;
   0,1,0,0];%Measurement matrix

SR=[100, 0;
    0,  100];%Square root covariance matrix.
R=SR*SR';%The covariance matrix.

z=zeros(2,numSamples);
for curSamp=1:numSamples
    z(:,curSamp)=H*xTrue(:,curSamp)+SR*randn(2,1);
end

%Plot the trajectory and the measurements
disp('Plotting the trajectory the measurements, and the true maneuver mode.')
figure(1)
clf
hold on
plot(xTrue(1,:),xTrue(2,:),'-b','linewidth',2)
scatter(z(1,:),z(2,:),'.k')

h1=xlabel('Meters West->East');
h2=ylabel('Meters South->North');
h3=title('The True Trajectory and the Detections');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')

%Plot when it is in Mode 2, 1 for being in it and 2 for not.
figure(2)
clf
hold on
plot(mode-1,'-k')
h1=xlabel('Discrete Step');
h2=ylabel('Maneuver mode');
h3=title('Indicating when the target is performing maneuvers.');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')

%STEP 3: Get the baseline scenario: Run a Kalman filter with a standard
%first-order white-noise process model, a Kalman filter with a first-order
%Gauss-Markov model, a reduced state estimator, and the
%separated covariance filter.
disp('Computing the baseline scenarios with just a Kalman filter and with the')
disp('reduced state estimator.')

%A direct-discrete model with a 1m/2s^2 process noise standard deviation is
%used (normally, a discretized model is preferred, because it is
%more consistent when using a variable sampling rate, and the covariance
%matrix is not singular).
Q=QPolyKalDirectDisc(T,4,1,1^2);%1m/s^2 process noise

%The state transition matrix for constant velocity motion.
F=FPolyKal(T,4,1);

%Parameters for the Gauss-Markov dynamic model. We are using a first-order
%Gauss-Markov model.
tau=20;%20 seconds; the assumed maneuver decorrelation time.
maxAccel=9.8*3;%Assume max 3G turn
%Rule of thumb- process noise suggestion.
q=processNoiseSuggest('PolyKal-ROT',maxAccel,T);
QGM=QGaussMarkov(T,4,q,tau,1);%Process noise matrix.
FGM=FGaussMarkov(T,4,tau,1);%State Transition matrix

%Allocate space for the polynomial Kalman filter estimates
xKalman=zeros(4,numSamples);%The states
PKalman=zeros(4,4,numSamples);%The covariance matrices

%Allocate space for the Gauss-Markov Kalman filter estimates
xKalmanGM=zeros(4,numSamples);%The states
PKalmanGM=zeros(4,4,numSamples);%The covariance matrices

%Allocate space for the reduced state estimator.
xRedState=zeros(4,numSamples);%The states

%Allocate space for the separated covariance filter.
xSCFState=zeros(4,numSamples);%The states

%We are using the version of the reduced state estimator that needs a
%maximum linear acceleration (directed along the direction of motion of the
%target) for the target along with a maximum turn rate acceleration.
%A 0.5m/s^2 maximum linear acceleration with 1 deg/s^2 maximum turn rate
%acceleration is used.
ARedState=0.5;
OmegaRedState=1*(pi/180);

%The separated covariance filter needs a matrix introducing the effects of
%the maximum acceleration onto the state (for covariance computation
%results. For a state consisting of position and velocity in 2D, this can
%be
aMax=9.8/2;
Ba=[T^2/2,   0;
    0,       T^2/2;
	T,       0;
	0,       T]*aMax;
%where we chose the maximum acceleration in each dimension to be 0.5G
%(1G=9.8m/s^2).

%Track initiation is by one-point differencing using just assumed maximum
%values as the standard deviations for the velocity covariance matrix. The
%assumed maximum value is 300m/s.
[xKalman(:,1),PKalman(:,:,1)]=onePointCartInit(z(:,1),SR,300);
xKalmanGM(:,1)=xKalman(:,1);
PKalmanGM(:,1)=PKalman(:,1);

%The reduced state estimator keeps the covariance matrix broken into parts
%M and D, which will be propagated and updated.
xRedState(1:2)=z(:,1);
M=PKalman(:,:,1);
D=zeros(4,2);

%The separated covariance filter. It has been chosen to start the filter
%with zero lag and all uncertainty in the covariance matrix.
xSCFState(1:2)=z(:,1);
PSCF=PKalman(:,:,1);
LSCF=zeros(4,2);

absErrKalman=zeros(numSamples,1);%Allocate space
absErrKalman(1)=norm(xKalman(1:2,1)-xTrue(1:2,1));

absErrKalmanGM=zeros(numSamples,1);%Allocate space
absErrKalmanGM(1)=norm(xKalman(1:2,1)-xTrue(1:2,1));

absErrRedState=zeros(numSamples,1);%Allocate space
absErrRedState(1)=norm(xRedState(1:2,1)-xTrue(1:2,1));

absErrSCF=zeros(numSamples,1);%Allocate space
absErrSCF(1)=norm(xSCFState(1:2,1)-xTrue(1:2,1));

%Run the filter
for curSamp=2:numSamples
    %The Kalman filter with the polynomial model
    %Predict the state forward
    [xPred, PPred]=discKalPred(xKalman(:,curSamp-1),PKalman(:,:,curSamp-1),F,Q);
    %Update the state with a measurement
    [xKalman(:,curSamp),PKalman(:,:,curSamp)]=KalmanUpdate(xPred,PPred,z(:,curSamp),R,H);
    
    absErrKalman(curSamp)=norm(xKalman(1:2,curSamp)-xTrue(1:2,curSamp));
    
    %The Kalman filter with the Gauss-Markov model
    [xPred, PPred]=discKalPred(xKalmanGM(:,curSamp-1),PKalmanGM(:,:,curSamp-1),FGM,QGM);
    %Update the state with a measurement
    [xKalmanGM(:,curSamp),PKalmanGM(:,:,curSamp)]=KalmanUpdate(xPred,PPred,z(:,curSamp),R,H);
    
    absErrKalmanGM(curSamp)=norm(xKalmanGM(1:2,curSamp)-xTrue(1:2,curSamp));
    
    %The reduced state estimator.
    %Predict the state forward
    modParams=[];
    modParams.T=T;
    modParams.A=ARedState;
    modParams.Omega=OmegaRedState;
    [xPred,M,D,PPred]=reducedStateDiscPred(xRedState(:,curSamp-1),M,D,[],[],'GenTurn',modParams);
    %Update the state with a measurement
    [xRedState(:,curSamp),M,D]=reducedStateUpdate(xPred,PPred,M,D,z(:,curSamp),R);
    
    absErrRedState(curSamp)=norm(xRedState(1:2,curSamp)-xTrue(1:2,curSamp));

    %The separated covariance filter.
    %Predict the state forward
    [xPred,LSCF,TSCF]=separatedCovDiscPred(xSCFState(:,curSamp-1),PSCF,LSCF,F,Ba);
    %Update the state with a measurement
    [xSCFState(:,curSamp),LSCF,PSCF]=separatedCovUpdate(xPred,LSCF,TSCF,z(:,curSamp),R,H);
    absErrSCF(curSamp)=norm(xSCFState(1:2,curSamp)-xTrue(1:2,curSamp));
end

figure(1)
hold on 
plot(xKalman(1,:),xKalman(2,:),'-g','linewidth',2)
plot(xKalmanGM(1,:),xKalmanGM(2,:),'-c','linewidth',2)
plot(xRedState(1,:),xRedState(2,:),'--b','linewidth',2)
plot(xSCFState(1,:),xSCFState(2,:),'-.m','linewidth',2)

figure(3)
clf
hold on
plot(absErrKalman,'-g','linewidth',2)
plot(absErrKalmanGM,'-c','linewidth',2)
plot(absErrRedState,'--b','linewidth',2)
plot(absErrSCF,'-.m','linewidth',2)

h1=xlabel('Discrete Step');
h2=ylabel('Absolute distance error');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')

%%STEP 4: Run the algorithms with the two linear models.
%The linear model pair consists of trackers using direct discrete linear
%dynamic models.
disp('Computing the IMM scenario with two linear Kalman filters with different process noises.')

%The covariance matrix for the low-noise model.
QLow=QPolyKalDirectDisc(T,4,1,0.1^2);%0.1m/s^2 process noise
%The covariance matrix for the high-noice model.
QHigh=QPolyKalDirectDisc(T,4,1,2^2);%2m/s^2 process noise.

%The propagation routines for each model
transFuns{1}=@(xPrev,PPrev)discKalPred(xPrev,PPrev,F,QLow);
transFuns{2}=@(xPrev,PPrev)discKalPred(xPrev,PPrev,F,QHigh);

%The measurement update function is the same for both models.
measUpdateFuns=@(x,P,z,R)KalmanUpdate(x,P,z,R,H);

%Allocate space for the IMM model states. This is the size of the
%largest state by the number of models, except when using the GPB1
%estimator, in which case it is just the size of the state (all states have
%to be the same size).
if(~strcmp(AlgSel,'GPB1'))
    xIMML=zeros(4,2);
    PIMML=zeros(4,4,2);

    %Initialize the states using one-point differencing. The covariance matrix
    %is set using the assumed velocity for each of the models, which is
    %300m/s^2.
    %The first model.
    [xIMML(:,1),PIMML(:,:,1)]=onePointCartInit(z(:,1),SR,300);
    %The second model has the same initialization.
    xIMML(:,2)=xIMML(:,1);
    PIMML(:,:,2)=PIMML(:,:,1);
else
    [xIMML,PIMML]=onePointCartInit(z(:,1),SR,300);
end

%Allocate space for the states to display. These will be the merged values
%from the IMM model states.
xIMMLDisp=zeros(4,numSamples);
PIMMLDisp=zeros(4,4,numSamples);

%Initially, both models have the same initialization, which will thus be
%the one to display at the first step.
xIMMLDisp(:,1)=xIMML(:,1);
PIMMLDisp(:,:,1)=PIMML(:,:,1);

%Allocate space for the mode probabilities. These are saved over all scans
%so that they can be plotted.
muMode=zeros(2,numSamples);%There is one for each model at each step.
%Initially, the mode probabilities are uniform: we do not know which mode
%it is in.
muMode(:,1)=1/2;

%The mode transition probability matrix for the two linear models must be
%determined. This is based on a matrix relating the mean sojourn times of
%the models and the time between samples.
%The mean sojourn time matrix is a design parameter. The negative inverse
%of the diagonal terms is the mean time spent in that state. Each row must
%sum to one. If the matrix were larger, the relation between the
%non-diagonal values in the rows would affect to which state it
%transitions. This is a design parameter.
A=[-0.01/2, 0.01/2;
    0.02*2, -0.02*2];
 
%Compute the mode transition probability matrix for the given time between
%samples. Because the sample interval is constant, this does not change 
%over time.
LambdaL=getMarkovPTransProbMat(A,T);

absErrIMML=zeros(numSamples,1);%Allocate space
absErrIMML(1)=norm(xIMMLDisp(1:2,1)-xTrue(1:2,1));
%Run the filter
for curSamp=2:numSamples
   %Predict the states forward 
   [xIMML,PIMML]=multipleModelPred(AlgSel,xIMML,PIMML,transFuns);
   
   %Update the state with a measurement
   [xIMML,PIMML,muMode(:,curSamp),xIMMLDisp(:,curSamp),PIMMLDisp(:,:,curSamp)]=multipleModelUpdate(AlgSel,xIMML,PIMML,z(:,curSamp),R,measUpdateFuns,muMode(:,curSamp-1),LambdaL);
   absErrIMML(curSamp)=norm(xIMMLDisp(1:2,curSamp)-xTrue(1:2,curSamp));
end

figure(1)
plot(xIMMLDisp(1,:),xIMMLDisp(2,:),'-r','linewidth',2)

figure(2)
plot(muMode(2,:),'-r','linewidth',2)

figure(3)
plot(absErrIMML,'-r','linewidth',2)


if(strcmp(AlgSel,'GPB1'))
    return
end


%STEP 5: Run the algorithms with the linear and the coordinated turn model.
disp('Computing the IMM scenario with a linear Kalman filter and an EKF.')

%The covariance matrix for the linear model.
QLinear=QPolyKalDirectDisc(T,4,1,0.01^2);%0.01m/s^2 process noise.

%0.025m/s^2 linear acceleration (in 3D, not just directed along the
%direction of motion of the target) with 0.5 deg/s^2 turning acceleration
%for the coordinated turn model as the standard deviations going into the
%covariance matrix. Note that these assumptions differ from those used in
%the reduced state filter. The filters should be tuned separately for a
%fair comparison. Often, improving the IMM's ability to recognize a turn
%worsens its overall track performance. Thus, a tradeoff must be present.
QCT=QCoordTurn(T,zeros(5,1),0.025^2,(0.5*(pi/180))^2);

%The propagation routines for each model
transFuns=[];
transFuns{1}=@(xPrev,PPrev)discKalPred(xPrev,PPrev,F,QLinear);
f=@(x)(FCoordTurn2D(T,x)*x);
transFuns{2}=@(xPrev,PPrev)discEKFPred(xPrev,PPrev,f,@(x)JacobCoordTurn2D(T,x),QCT);

%The measurement update function differs for the different models, because
%H is different for the coordinated turn model, due to the extra element
%in the state.
HCT=[1,0,0,0,0;
     0,1,0,0,0];
measUpdateFuns=[];
measUpdateFuns{1}=@(x,P,z,R)KalmanUpdate(x,P,z,R,H);
measUpdateFuns{2}=@(x,P,z,R)KalmanUpdate(x,P,z,R,HCT);

%The number of elements in each state and the number of dimensions of each
%state that play a role in mixing the states.
numStateDims=[4;5];
numMixDims=[4;4];

%Allocate space for the IMM model states. The number of rows is the size of
%the longest state.
xLCT=zeros(5,2);
PLCT=zeros(5,2);

%Initialize the states using one-point differencing. The covariance matrix
%is set using the assumed velocity for each of the models, which is
%300m/s^2 and the maximum turn rate, which is assumed to be 20 degrees per
%second.
%The first model.
xLCT(1:2,1)=z(:,1);
PLCT(1:2,1:2,1)=R;
PLCT(3,3,1)=300^2;
PLCT(4,4,1)=300^2;
%The second model has the same initialization.
xLCT(:,2)=xLCT(:,1);
PLCT(:,:,2)=PLCT(:,:,1);
PLCT(5,5,2)=(20*(pi/180))^2;

%Allocate space for the states to display. These will be the merged values
%from the IMM model states. The size is the size of the state with the
%least number of mixing components.
xLCTDisp=zeros(4,numSamples);
PLCTDisp=zeros(4,4,numSamples);

%Initially, both models have the same initialization, which will thus be
%the one to display at the first step.
xLCTDisp(:,1)=xLCT(1:4,1);
PLCTDisp(:,:,1)=PLCT(1:4,1:4,1);

%Allocate space for the mode probabilities. These are saved over all scans
%so that they can be plotted.
muModeLCT=zeros(2,numSamples);%There is one for each model at each step.
%Initially, the mode probabilities are uniform: we do not know which mode
%it is in.
muModeLCT(:,1)=1/2;

%Compute the mode transition probability matrix.
A=[-0.01/2, 0.01/2;
    0.02*2, -0.02*2];
LambdaCT=getMarkovPTransProbMat(A,T);

absErrLCT=zeros(numSamples,1);%Allocate space
absErrLCT(1)=norm(xLCTDisp(1:2,1)-xTrue(1:2,1));
%Run the filter
for curSamp=2:numSamples
   %Predict the states forward.
   [xLCT,PLCT]=multipleModelPred(AlgSel,xLCT,PLCT,transFuns,numStateDims,numMixDims);

   %Update the state with a measurement.
   [xLCT,PLCT,muModeLCT(:,curSamp),xLCTDisp(:,curSamp),PLCTDisp(:,:,curSamp)]=multipleModelUpdate(AlgSel,xLCT,PLCT,z(:,curSamp),R,measUpdateFuns,muModeLCT(:,curSamp-1),LambdaCT,numStateDims,numMixDims);  
   absErrLCT(curSamp)=norm(xLCTDisp(1:2,curSamp)-xTrue(1:2,curSamp));
end

figure(1)
plot(xLCTDisp(1,:),xLCTDisp(2,:),'--c','linewidth',2)
legend('True Trajectory','Raw Detections','Basic Kalman Filter', 'Gauss-Markov Kalman Filter','Reduced State Filter', 'Separated Covariance Filter', 'Two Linear Models','Linear and Maneuvering Models','Location','SouthEast') 

figure(2)
plot(muModeLCT(2,:),'--c','linewidth',2)
legend('True Mode', 'Two Linear Models','Linear and Maneuvering Models','Location','NorthWest') 

figure(3)
plot(absErrLCT,'--c','linewidth',2)
legend('Polynomial Kalman Filter','Gauss-Markov Kalman Filter','Reduced State Filter', 'Separated Covariance Filter', 'Two Linear Models','Linear and Maneuvering Models','Location','NorthWest') 
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
