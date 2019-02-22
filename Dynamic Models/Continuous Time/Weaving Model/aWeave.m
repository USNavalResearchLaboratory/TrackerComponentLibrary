function [aVal,aJacob,aHess,papt]=aWeave(x,t,A,alpha,theta0,uTurn)
%%AWEAVE The drift function for a somewhat sinusoidal weaving motion model
%        under a 3D flat-Earth approximation. The third example below shows
%        how this flat-Earth model can be used on a curved Earth.
%
%INPUTS: x The 6X1 target state at time t consisting of 3 position and 3
%          velocity components.
%        t The time at which the drift function is to be evaluated.
%        A A scalar term determining the magnitude of the deflection during
%          turns. It should be less than alpha*pi/2 to keep the target from
%          doubling back on its previous path.
%    alpha The positive scalar weave period in radians per second.
%   theta0 The scalar phase at time 0 of the weave. A cosine function is
%          used for the weave. The default if this parameter is omitted or
%          an empty matrix is passed is 0. To start a track going in the
%          overall direction of progression of the target, set theta=0 and
%          the initial velocity direction to the overall direction that you
%          want the target to go (weaves are orthogonal to that direction).
%    uTurn A 3X1 unit vector pointing in the direction of the turn axis
%          used for the weaves. It is orthogonal to the turn plane. For
%          example, if a target is going to weave in the x-y plane, then
%          this should be [0;0;1]. The default if this parameter is omitted
%          or an empty matrix is passed is [0;0;1].
%
%OUTPUTS: aVal The 6X1 flat-Earth time-derivative of the state.
%       aJacob The 6X6 matrix of partial derivatives of aVals such that
%              aJacob(:,k) is the partial derivative of aVals(:,k) with
%              respect to x(k).
%        aHess The 6X6X6 matrix of second partial derivatives of aVals such
%              that aHess(:,k1,k2) is the second partial derivative of
%              aVals with respect to x(k1) and x(k2).
%         papt The 6X1 partial derivative with resect to time of aVals.
%
%This model is derived in [1]. It is parameterized in a manner that allows
%one to create a target trajectory that performed weaves of a desired
%amplitude in a particular direction, progressing at a desired speed
%overall.
%
%Equation 78 and 79 in [1] are for a general turning model. Equation 88
%and 89 specifically show how the velocity is modulated, which in turn can
%lead to using a time-varying angular momentum vector modified from 79 to
%Omega=A*cos(alpha*t+theta0)*uTurn. Unlike Equation 79, this is generalized
%to allow uTurn to be any rotation axis, not just the z-axis. This is the
%model as implemented here.
%
%EXAMPLE 1:
%Here, we implement code that draws Figure 11 from [1]. We want the overall
%motion direction of the target to be along the x-axis. We also want the
%target to fly travel at a specified speed and go a specified distance
%along the x axis over the observation period. The one parameter that is
%left free is the magnitude of the weaving. Here, we plot two trajectories
%that have the same starting and ending points and make the same number of
%weaves, but have different weaving magnitudes. Since both targets are
%going the same speed, we need to determine the time duration of the
%trajectories so that we can show them advancing the same amount along the
%x-axis.
% Nw=6;%The number of weave cycles.
% speed=10;%Meters per second
% s=100;%The desired distance we want to cover along the x-axis (meters).
% uTurn=[0;0;1];%Turns are in the x-y plane.
% %The initial velocity of the target is along the x axis. By setting
% %theta0-0 and the velocity of xInit to point in the desired direction of
% %progression, we can easily set the direction of progression.
% theta0=0;
% xInit=[0;0;0;speed*[1;0;0]];
% 
% %Parameters for the first target
% %beta determines the weave magnitude, as in Equation 101 of [1]. This value
% %can be from 0 to 1.
% betaVal1=1;
% %Solve Equation 102 in [1] to determine how long one must travel to go a
% %distance of s meters. Also determine the weave period in meters per second
% %as well as the amplitude parameter A of the weaves.
% [tEnd1,alphaVal1,A1]=determineWeaveTime(s,speed,Nw,betaVal1);
% 
% %Parameters for the second target. This target does not weave as strongly.
% betaVal2=1/2;
% [tEnd2,alphaVal2,A2]=determineWeaveTime(s,speed,Nw,betaVal2);
% 
% %Now, create the trajectories. We use a Runge-Kutta algorithm with a fixed
% %step to get the values. The step has been chosen to be sufficiently small.
% %However, if we only wanted a few samples we could use RungeKAtTimes for a
% %fixed-step size algorithm, or we could use RKAdaptiveAtTimes.
% numSteps=3000;
% deltaT1=tEnd1/numSteps;
% deltaT2=tEnd2/numSteps;
% aDyn1=@(x,t)aWeave(x,t,A1,alphaVal1,theta0,uTurn);
% xList1=RungeKSteps(xInit,0,aDyn1,deltaT1,numSteps);
% aDyn2=@(x,t)aWeave(x,t,A2,alphaVal2,theta0,uTurn);
% xList2=RungeKSteps(xInit,0,aDyn2,deltaT2,numSteps);
% 
% figure(1)
% clf
% hold on
% plot(xList1(1,:),xList1(2,:),'-r','linewidth',6)
% plot(xList2(1,:),xList2(2,:),'--g','linewidth',6)
% axis([0 100 0 14])
% h1=xlabel('x');
% h2=ylabel('y');
% 
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%In this example, we want a target to travel a distance of s meters in tMax
%seconds while weaving Nw times. As in the previous example, the target is
%moving along the x axis. We ned to use the function determineWeaveSpeed to
%determine the speed of the target that is required to make this happen. 
% Nw=4;
% betaVal=0.5;
% tEnd=180;%The time over which the distance was traveled.
% theta0=0;
% uTurn=[0;0;1];%Turns are in the x-y plane.
% s=1.085e6;%The distance traveled.
% 
% [speed,alphaVal,A]=determineWeaveSpeed(s,tEnd,Nw,betaVal);
% xInit=[0;0;0;speed*[1;0;0]];
% aDyn=@(x,t)aWeave(x,t,A,alphaVal,theta0,uTurn);
% 
% numSteps=1000;
% deltaT=tEnd/numSteps;
% xList=RungeKSteps(xInit,0,aDyn,deltaT,numSteps);
% figure(2)
% clf
% hold on
% plot(xList(1,:),xList(2,:),'-r','linewidth',6)
% axis([0 s 0 8e4])
% h1=xlabel('x');
% h2=ylabel('y');
% 
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 3:
%Here, we demonstrate how to use this flat-Earth model on a curved Earth.
%Some place in the ocean off of Hawaii.
% latLonStart=[20.3924;-155.4976]*(pi/180);
% %The location of the Ahihi-Kinau Natural Area Reserve
% latLonEnd=[20.6050;-156.4353]*(pi/180);
% %The altitude of the target-- it will remain constant.
% altitudeInit=20e3;
% 
% tMax=180;
% times=0:tMax;
% 
% %Get an initial heading and distance. We are ignoring the altitude of the
% %target, though it would be taken into account using the (slower) function
% %indirectGeodeticProbGen
% [azi1,s]=indirectGeodeticProb(latLonStart,latLonEnd);
% 
% pStart=[phi0;lambda0;altitudeInit];
% uInit=getENUAxes(pStart);
% 
% xyzStart=ellips2Cart(pStart);
% %The initial heading in the local tangent plane coordinates.
% uh=[sin(azi1);cos(azi1);0];
% 
% %These two parameters make the sinusoidal deviations approximately
% %centered about the nominal trajectory 
% theta0=0;
% Nw=4;%The number of weave cycles.
% beta=0.5;
% %We want the plane to travel a distance of about s in tMax seconds
% %while weaving Nw times. This means that the speed along the
% %trajectory is
% [speed,alphaVal,A]=determineWeaveSpeed(s,tMax,Nw,beta);
% xInit=[xyzStart;speed*uh];
% 
% %To determine the realism of this being a manned maneuver, since
% %this is level flight, the maximum G-force felt can be easily found
% %to be
% %MaxGForce=sqrt((A)^2*norm(xInit(4:6,end))^2+9.81^2)/9.81
% %which is 7.4G and is within a reasonable tolerance for a pilot
% %wearing an "anti-G" suit.
% 
% uTurn=[0;0;1];%The local vertical is the turn axis.
% aDyn=@(x,t)aWeave(x,t,A,alphaVal,theta0,uTurn);
% 
% deltaTMax=0.05;
% [xList,uList]=RungeKCurvedAtTimes(xInit,uInit,times,aDyn,[],deltaTMax);
% xList=xList(1:6,:);
% %If one wanted to get the velocity information the global coordinates,
% %then the following line is necessary:
% xList(4:6,:)=getGlobalVectors(xList(4:6,:),uList);
% 
% %Bounds for the axes on the map in case this trajectory is plotted.
% axisMapBounds=[-156-0.6,-156+0.6,20.5-0.6,20.5+0.6];
% axisHeightBounds=[0 times(end) 0 25e3];
% 
% %Convert the Cartesian trajectory into latitude, longitude and
% %altitude.
% points=Cart2Ellipse(xList(1:3,:));
% 
% figure(1)
% clf
% hold on
% axis(axisMapBounds)
% axis square
% 
% %Plot the trajectory of the aircraft with the curved-Earth model.
% plot(points(2,:)*180/pi,points(1,:)*180/pi,'-r','linewidth',6);
% 
% h1=xlabel('Longitude');
% h2=ylabel('Latitude');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
% 
% %Plot the altitude as a function of time.
% figure(2)
% clf
% plot(times,points(3,:),'-r','linewidth',6);
% axis(axisHeightBounds)
% 
% h1=xlabel('Time (Seconds)');
% h2=ylabel('Altitude (Meters)');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<6||isempty(uTurn))
        uTurn=[0;0;1]; 
    end
    
    if(nargin<5||isempty(theta0))
       theta0=0; 
    end

    Omega=A*cos(alpha*t+theta0)*uTurn;
    aVal=[x(4:6);cross(Omega,x(4:6))];
    
    if(nargout>1)
        OmegaX=Omega(1);
        OmegaY=Omega(2);
        OmegaZ=Omega(3);
        
        aJacob=[0,0,0,1,        0,          0;
                0,0,0,0,        1,          0;
                0,0,0,0,        0,          1;
                0,0,0,0,        -OmegaZ,    OmegaY;
                0,0,0,OmegaZ,   0,          -OmegaX;
                0,0,0,-OmegaY,  OmegaX,     0];

        if(nargout>2)
            aHess=zeros(6,6,6);
            if(nargout>3)
                sinTerm=sin(alpha*t+theta0);
                
                uTurnX=uTurn(1);
                uTurnY=uTurn(2);
                uTurnZ=uTurn(3);
                
                vx=x(4);
                vy=x(5);
                vz=x(6);
                
                papt=[0;
                      0;
                      0;
                      A*alpha*(uTurnZ*vy-uTurnY*vz)*sinTerm;
                      A*alpha*(-uTurnZ*vx+uTurnX*vz)*sinTerm;
                      A*alpha*(uTurnY*vx-uTurnX*vy)*sinTerm];
            end
        end
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
