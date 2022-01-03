function demoNonlinearMeasFilters()
%%DEMONONLINEARMEASFILTERS This demonstrates the use of basic target
%               tracking filters on a simple scenario involving nonlinear
%               measurements and a linear dynamic model. The instance of a
%               moving observer tracking an object using in 2D angle-only
%               measurements is considered. The object is observable if the
%               observer is more maneuverable than the object.
%
%The problem of tracking using 2D angular measurements is given as an
%example in Chapter 10.3.4 of [1]. However, the example here is closer to
%that used in [2].
%
%A platform moves around, in this instances making instantaneous direction
%changes, and takes angle-only measurements of a sensor. To demonstrate the
%algorithms, no actual track initiation routine is used. Rather, a random
%initialization (given an accurate prior) is used. 
%
%The use of angular measurements also demonstrates the need for the
%innovation and measurement averaging transformations as described in [3],
%which apply to the pure propagation filter, the square root cubature
%Kalman filter, and the quasi-Monte Carlo Kalman Filter. The same
%transformation also applies to variants of the EKF, though it might take a
%few Monte Carlo runs to see how very large errors can arise when not
%properly wrapping innovation values.
%
%Running this displays the tracks from twelve filters for the scenario. In
%this instance, all of the algorithms have comparable accuracy.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%[2] O. Straka, J. Duník, and M. Simandl, "Design of pure propagation
%    unscented Kalman filter," in Proceedings of the 19th World Congress of
%    The International Federation of Automatic Control, Cape Town,
%    South Africa, 24-29 Aug. 2014, pp. 5399-5938.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

T=60;%60 second sampling interval

F=FPolyKal(T,4,1);%State transition matrix.
Gamma=[0.5*T^2, 0;
         0,     0.5*T^2;
         T,     0;
         0,     T];
%Singular process noise covariance with units of m^2/s^4;
sigma=sqrt(1.6e-6);
QReduced=sigma^2*eye(2,2);
SQReduced=chol(QReduced,'lower');
Q=Gamma*QReduced*Gamma';
SQ=cholSemiDef(Q,'lower');

%Measurement noise covariance.
R=(1*(pi/180))^2;%1 degree standard deviation (converted to radians).
RInv=inv(R);

%Create the true course of the target and of the platform.
numSteps=50;

%Initial target heading -140 degrees (clockwise from y axis).
theta=-140*(pi/180);
v=convertSpeedUnits(4,'nml','h','m','s');%knots to m/2
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xTargetInit=[6e3;0*2e3;xVel;yVel];

%The target just goes in a straight line.
xTarget=zeros(4,numSteps);
xTarget(:,1)=xTargetInit;
for k=2:numSteps
    xTarget(:,k)=F*xTarget(:,k-1)+SQ*randn(4,1);
end

%Initial heading +140 degrees (clockwise from y axis).
theta=140*(pi/180);
v=convertSpeedUnits(5,'nml','h','m','s');%knots to m/2
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatformInit=[0;0;xVel;yVel;0.002];

xPlatform=zeros(5,numSteps);
xPlatform(:,1)=xPlatformInit;
for k=2:numSteps
    xPlatform(:,k)=FPolarCoordTurn2D(T,xPlatform(:,k-1))*xPlatform(:,k-1);
end

%Create the trajectory of the platform. It makes instantaneous turns.
xPlatform=zeros(4,numSteps);
xPlatform(:,1)=xPlatformInit(1:4);
for k=2:13
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

theta=18*(pi/180);
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatform(3:4,13)=[xVel;yVel];
for k=14:24
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

theta=90*(pi/180);
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatform(3:4,24)=[xVel;yVel];
for k=25:30
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

theta=180*(pi/180);
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatform(3:4,30)=[xVel;yVel];
for k=31:35
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

theta=-18*(pi/180);
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatform(3:4,35)=[xVel;yVel];
for k=36:40
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

theta=-60*(pi/180);
xVel=v*cos(theta-pi/2);
yVel=-v*sin(theta-pi/2);
xPlatform(3:4,40)=[xVel;yVel];
for k=41:numSteps
    xPlatform(:,k)=F*xPlatform(:,k-1);
end

%Create the measurements
z=zeros(1,numSteps);
for k=1:numSteps
    z(:,k)=measFunc(xTarget(:,k),xPlatform(:,k))+sqrt(R)*randn(1,1);
end

%Just use a simple random initialization.
PInit=[1e5,    0,    0,   0;
       0,      1e5   0,   0;
       0,      0,    2    0;
       0,      0,    0    2];
SInit=chol(PInit,'lower');
xInit=xTarget(:,1)+SInit*randn(4,1);%Initial predicted value.

%The cubature points for most filters will be fifth-order.
[xi,w]=fifthOrderCubPoints(4);

%The cubature points for the ensemble Kalman filter will be a set of the
%same number of points as used for the other filters, but with uniform
%weights.
xiLCD=GaussianLCDSamples(4,length(w));

f=@(x)(F*x);%The dynamic model.
%The wrapping function for the innovation.
innovTrans=@(a,b)wrapRange(bsxfun(@minus,a,b),-pi,pi);
%The function for taking an average of weighted angular measurements.
measAvgFun=@(z,w)meanAng(z,w');
%The function for taking an average of unweighted angular measurements.
measAvgFunUnweighted=@(z)meanAng(z);

%The pure propagation filter
xEstPP=zeros(4,numSteps);
xEstPP(:,1)=xInit;
%The initial set of points to propagate for the pure propagation filter.
xiPP=transformCubPoints(xi,xInit,SInit);

%The square root CKF
xEstCKF=zeros(4,numSteps);
xEstCKF(:,1)=xInit;
SCKF=chol(PInit,'lower');
SR=chol(R,'lower');

%The extended Kalman filter
xEstEKF=zeros(4,numSteps);
xEstEKF(:,1)=xInit;
PEKF=PInit;

%The square root extended Kalman filter
xEstSEKF=zeros(4,numSteps);
xEstSEKF(:,1)=xInit;
SEKF=chol(PInit,'lower');

%The extended square root information filter
xEstESRIF=zeros(4,numSteps);
PInvSqrtESRIF=inv(chol(PInit,'lower'));
yEstESRIF=PInvSqrtESRIF*xInit;
xEstESRIF(:,1)=xInit;

%The Ensemble Kalman Filter. It needs a set of points to start. These have
%to all have the same weighting, unlike in the pure propagation filter.
%Thus, we will generate from cubature points from the GaussianLCDSamples
%function. We will use the same number of points as the pure propagation
%filter uses.
xEstEnsemb=zeros(4,numSteps);
xEstEnsemb(:,1)=xInit;
xEnsemb=bsxfun(@plus,SInit*xiLCD,xInit);

%The progressive Gaussian filter (PGF). This filter is more appropriate for
%non-Gaussian measurements, but it can be used here as well.
xEstPGF=zeros(4,numSteps);
xEstPGF(:,1)=xInit;
SPGF=chol(PInit,'lower');

%The Gaussian particle filter (GPF). This filter is also more appropriate
%for non-Gaussian measurements. 1e4 particles are used.
numPartGPF=1e3;
xEstGPF=zeros(4,numSteps);
xEstGPF(:,1)=xInit;
PGPF=PInit;

%The cubature information filter (CIF)
xEstCIF=zeros(4,numSteps);
xEstCIF(:,1)=xInit;
PInvCIF=inv(PInit);
yEstCIF=PInvCIF*xEstCIF(:,1);

%The extended information filter (EIF)
xEstEIF=zeros(4,numSteps);
xEstEIF(:,1)=xInit;
PInvEIF=inv(PInit);
yEstEIF=PInvEIF*xEstEIF(:,1);

%The first-order square-root divided difference filter
xEstDDF1=zeros(4,numSteps);
xEstDDF1(:,1)=xInit;
SDDF1=chol(PInit,'lower');

%The second-order square-root divided difference filter
xEstDDF2=zeros(4,numSteps);
xEstDDF2(:,1)=xInit;
SDDF2=chol(PInit,'lower');

%The quasi-Monte Carlo Kalman Filter 
xEstQMC=zeros(4,numSteps);
xEstQMC(:,1)=xInit;
PQMC=PInit;
numQMCSamples=500;

for curStep=2:numSteps
    %The measurement function for all filters.
    h=@(xState)(measFunc(xState,xPlatform(:,curStep)));
    %The Jacobian for the EKF.
    HJacob=@(xState)(measFuncJacob(xState,xPlatform(:,curStep)));
    
    %%Predict and update the pure propagation filter.
    xiPP=purePropDiscPred(xiPP,w,f,Q);
    [xiPP,xEstPP(:,curStep)]=purePropUpdate(xiPP,w,z(curStep),R,h,innovTrans,measAvgFun);
    
    %%Predict and update the square root Cubature Kalman filter.
    [xPred, SCKF]=sqrtDiscKalPred(xEstCKF(:,curStep-1),SCKF,F,SQ);
    [xEstCKF(:,curStep), SCKF]=sqrtCubKalUpdate(xPred,SCKF,z(curStep),SR,h,xi,w,innovTrans,measAvgFun);
    
    %Predict and update the EKF
    [xPred, PEKF]=discKalPred(xEstEKF(:,curStep-1),PEKF,F,Q);
    [xEstEKF(:,curStep),PEKF]=EKFUpdate(xPred,PEKF,z(curStep),R,h,HJacob,[],[],innovTrans);
    
    %Predict and update the square root EKF
    [xPred, SEKF]=sqrtDiscKalPred(xEstSEKF(:,curStep-1),SEKF,F,SQ);
    [xEstSEKF(:,curStep),SEKF]=sqrtEKFUpdate(xPred,SEKF,z(curStep),SR,h,HJacob,innovTrans);
    
    %Predict and update the extended square root information filter.
    [ySqrtPred, PInvSqrtESRIF]=sqrtInfoFilterDiscPred(yEstESRIF,PInvSqrtESRIF,F,SQReduced,[],Gamma);
    [yEstESRIF,PInvSqrtESRIF]=ESRIFUpdate(ySqrtPred,PInvSqrtESRIF,z(curStep),SR,h,HJacob,innovTrans);
    %Get the state estimate from the square root information estimate.
    xEstESRIF(:,curStep)=PInvSqrtESRIF\yEstESRIF;
    
    %Predict and update the ensemble Kalman filter.
    xEnsemb=EnKFDiscPred(xEnsemb,f,SQ);
    [xEnsemb,xEstEnsemb(:,curStep)]=EnKFUpdate(xEnsemb,z(curStep),SR,h,0,innovTrans,measAvgFunUnweighted);
    
    %Predict and update the progressive Gaussian filter. Fifth-order
    %cubature points are used.
    [xPred, SPGF]=sqrtDiscKalPred(xEstPGF(:,curStep-1),SPGF,F,SQ);
    %The PGF needs the conditional likelihood of the state given the
    %measurement. However, the conditional likelihood does not need to be
    %normalized. Thus, we will use an unnormalized approximation to the
    %wrapped normal distribution.
    zLike=@(x)measLike(z(curStep),x,RInv,xPlatform(:,curStep));
    [xEstPGF(:,curStep),SPGF]=progressivGaussUpdate(xPred,SPGF,zLike,xi,w);
    
    %Predict and update the Gaussian particle filter.
    [xPred, PGPF]=discKalPred(xEstGPF(:,curStep-1),PGPF,F,Q);
    %The GPF uses the same likelihood state ad the PGF.
    zLike=@(x)measLike(z(curStep),x,RInv,xPlatform(:,curStep));
    %The predicted PDF is used as the importance sampling density.
    [xEstGPF(:,curStep),PGPF]=GaussPartFilterUpdate(xPred,PGPF,zLike,numPartGPF);
    
    %Predict and update the cubature information filter.
    [yPred, PInvPred]=infoFilterDiscPred(yEstCIF,PInvCIF,F,Q);
    [yEstCIF,PInvCIF]=cubInfoUpdate(yPred,PInvPred,z(curStep),RInv,h,xi,w,innovTrans,measAvgFun);
    %Get the state estimate from the information estimate.
    xEstCIF(:,curStep)=PInvCIF\yEstCIF;

    %Predict and update the extended information filter.
    [yPred, PInvPred]=infoFilterDiscPred(yEstEIF,PInvEIF,F,Q);
    [yEstEIF,PInvEIF]=EIFUpdate(yPred,PInvPred,z(curStep),RInv,h,HJacob,innovTrans);
    %Get the state estimate from the information estimate.
    xEstEIF(:,curStep)=PInvEIF\yEstEIF;
    
    %Predict and update the first-order divided difference filter.
    [xPred,SPred]=sqrtDiscKalPred(xEstDDF1(:,curStep-1),SDDF1,F,SQ);
    [xEstDDF1(:,curStep),SDDF1]=sqrtDDFUpdate(xPred,SPred,z(curStep),SR,h,1,innovTrans);
    
    %Predict and update the second-order divided difference filter.
    [xPred,SPred]=sqrtDiscKalPred(xEstDDF2(:,curStep-1),SDDF2,F,SQ);
    [xEstDDF2(:,curStep),SDDF2]=sqrtDDFUpdate(xPred,SPred,z(curStep),SR,h,2,innovTrans);
    
    %Predict and update the quasi-Monte Carlo filter
    [xPred, PQMC]=discKalPred(xEstQMC(:,curStep-1),PQMC,F,Q);
    [xEstQMC(:,curStep),PQMC]=QMCKalUpdate(xPred,PQMC,z(curStep),R,h,numQMCSamples,innovTrans,measAvgFun);
end

figure(1)
clf
hold on
plot(xPlatform(1,:)/1e3,xPlatform(2,:)/1e3,'-r','linewidth',4)%The sensor
plot(xTarget(1,:)/1e3,xTarget(2,:)/1e3,'--b','linewidth',4)%The target

plot(xEstPP(1,:)/1e3,xEstPP(2,:)/1e3,'-c','linewidth',2)%The estimate from pure propagation
plot(xEstCKF(1,:)/1e3,xEstCKF(2,:)/1e3,'-m','linewidth',2)%The estimate from the CKF
plot(xEstEKF(1,:)/1e3,xEstEKF(2,:)/1e3,'-y','linewidth',2)%The estimate from the EKF
plot(xEstSEKF(1,:)/1e3,xEstSEKF(2,:)/1e3,'-k','linewidth',2)%The estimate from the square root EKF
plot(xEstESRIF(1,:)/1e3,xEstESRIF(2,:)/1e3,'-r','linewidth',2)%The estimate from the extended square root information filter
plot(xEstEnsemb(1,:)/1e3,xEstEnsemb(2,:)/1e3,'-g','linewidth',2)%The estimate from the ensemble information filter
plot(xEstPGF(1,:)/1e3,xEstPGF(2,:)/1e3,'-b','linewidth',2)%The estimate from the progressive Gaussian filter
plot(xEstGPF(1,:)/1e3,xEstGPF(2,:)/1e3,'-.c','linewidth',2)%The estimate from the Gaussian particle filter
plot(xEstCIF(1,:)/1e3,xEstCIF(2,:)/1e3,'-.m','linewidth',2)%The estimate from the cubature information filter
plot(xEstEIF(1,:)/1e3,xEstEIF(2,:)/1e3,'-.y','linewidth',2)%The estimate from the extended information filter
plot(xEstDDF1(1,:)/1e3,xEstDDF1(2,:)/1e3,'-.k','linewidth',2)%The estimate from the first order divided difference filter.
plot(xEstDDF2(1,:)/1e3,xEstDDF2(2,:)/1e3,'-.r','linewidth',2)%The estimate from the second order divided difference filter.
plot(xEstQMC(1,:)/1e3,xEstQMC(2,:)/1e3,'--g','linewidth',2)%The estimate from the quasi-Monte Carlo Kalman filter.

h1=xlabel('x (km)');
h2=ylabel('y (km)');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')

legend('Sensor Trajectory','Target Trajectory','Pure Propagation Filter', 'Square Root CKF','EKF','Square Root EKF', 'ESRIF','Ensemble KF','PGF','GPF','CIF','EIF','DDF1','DDF2','QMC','Location','SouthWest') 

end

function z=measFunc(xTarget,xPlatform)
%The measurement function; it returns an angle-only measurement between the
%platform and the target.
    x=xTarget(1);
    y=xTarget(2);
    xp=xPlatform(1);
    yp=xPlatform(2);

    z=atan2((x-xp),(y-yp));
end

function J=measFuncJacob(xTarget,xPlatform)
    x=xTarget(1);
    y=xTarget(2);
    xp=xPlatform(1);
    yp=xPlatform(2);

    J=zeros(1,4);
    J(1)=(y-yp)/((x-xp)^2+(y-yp)^2);
    J(2)=(-x+xp)/((x-xp)^2+(y-yp)^2);
end

function likeVals=measLike(z,x,RInv,xPlatform)
%The non-normalized likelihood function of the measurement given the state,
%evaluated for a matrix of states. Since a Gaussian random variable is
%added to the angular measuremnet as noise, the result is a wrapped normal
%distribution. We will approximate it by just wrapping the difference in
%the normal PDF; this avoids an infinite sum and works well unless to
%variance is large.
numX=size(x,2);

likeVals=zeros(numX,1);
for curX=1:numX
    diff=wrapRange(z-measFunc(x(:,curX),xPlatform),-pi,pi);

    likeVals(curX)=exp(-(1/2)*diff'*RInv*diff);
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
