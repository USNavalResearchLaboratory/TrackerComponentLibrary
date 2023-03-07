function [a,c,xRange]=findHermiteInterpPolySet(x,y,numDims,groupSize)
%%FINDHERMITEINTERPPOLYSET Given the output of a multidimensional function
%           y(x) that takes a scalar parameter x at a number of points,
%           generate a set of interpolation coefficients so that one may
%           perform Hermite interpolation over the region covered by the
%           samples (and ot a limited extent, outside the region). The
%           function HermiteSetInterpVal is used to perform interpolation
%           using the output of this function. This function is useful when
%           one wishes to interpolate a multidimensional target state, for
%           example, consisting of position, velocity and acceleration in
%           [x;y;z] when given samples of the state. The difference between
%           this function and HermiteInterpMultiPoly is that this function
%           produces multiple low order Hermite interpolating polynomials,
%           whereas HermiteInterpMultiPoly tries to produce a single high
%           order one fitting all points. However, trying to fit too many
%           points will eventually lead to finite precision errors becoming
%           very significant.
%
%INPUTS: x An NpX1 or 1XNp vector of real, scalar values at which the
%          function y and its derivative are given. The values are assumed
%          to be given in increasing order. Np>=2.
%        y A (numMoments*numDims)XNp matrix of values of the
%          multidimensional function y(x) and its derivatives evaluated at
%          the points in x. All values for a particular derivative are
%          given before any of the next derivative. For example, to
%          utilize 3D position and velocity information, the ordering would
%          be [x;y;z;xDot;yDot;zDot]. Note that even if one is only
%          interested in interpoalting position values, passing known
%          moments (velocity, acceleration, jerk, etc.) will improve the
%          interpolation of the position values.
%  numDims The scalar number of dimensions present. This will typically be
%          2 or 3 when dealing with target states. If this parameter is
%          omitted or an empty matrix is passed, then numDims=1 is used.
% groupSize To perform interpolation between points x(i) and x(i+1), the
%          values of y at those points must be taken. However, extra points
%          can be used to improve the estimates. This is the total number
%          of points to use. floor((groupSize-1)/2) are used before x(i)
%          and ceil((groupSize-1)/2) are used after x(i). The minimum value
%          of groupSize is 2. If this parameter is omitted or an empty
%          matrix is passed, then 2 is used. Using a large group allows for
%          the interpolation of higher-order moments.
%
%OUTPUTS: a The numDimsXnumCoeffXnumInterpRegions hypermatrix of the
%           interpolating polynomial coefficients for each dimension and
%           region. This and the following two other outputs can be used in
%           the function HermiteSetInterpVal to perform interpolation.
%         c The numDimsX(numCoeff-1)XnumInterpRegions hypermatrix of
%           control points for the interpolating polynomials.
%    xRange The length (numInterpRegions+1) vector of boundary points of
%           the interpolating regions when using the HermiteSetInterpVal
%           function. The vector has the same orientation as x.
%
%This function uses a sliding window of points for generating multivariate
%Hermite interpolating polynomials using the HermiteInterpMultiPoly
%function. The interpolating polynomials will match all given points. All
%of the coefficients for the polynomials are stored in a and c. The values
%in xRange indicate where the interpolation regions for the windows are
%centered, so that the function HermiteSetInterpVal can choose the correct
%values in a and c to perform interpolation.
%
%EXAMPLE:
%Here, we consider the case where we are given samples of the position,
%velocity and acceleration of a target in a flat-Earth weaving trajectory.
%We would like to interpolate state vectors consisting of position,
%velocity and acceleration. To demonstrate that one can interpolate moments
%higher than those provided to this function, we only provide position and
%velocity, even though we interpolate acceleration. Here, we generate a
%trajectory 2000 samples long. We then only take every 25th sample.
%This function is used to generate the interpolating polynomials. The true
%and interpolated trajectories along with the samples values are plotted.
%Also, the interpolation error in position, velocity, and acceleration are
%plotted. One will see that the interpolation error is quite low
%considering the scale of the values.
%
% %The following are all parameters for the dynamic model.
% Nw=4;
% betaVal=0.5;
% t0=0;
% tEnd=180;%The time over which the distance was traveled.
% theta0=0;
% uTurn=[0;0;1];%Turns are in the x-y plane.
% s=1.085e6;%The distance traveled.
% [speed,alphaVal,A]=determineWeaveSpeed(s,tEnd,Nw,betaVal);
% xInit=[0;0;0;speed*[1;0;0]];
% aDyn=@(x,t)aWeave(x,t,A,alphaVal,theta0,uTurn);
% numDims=3;%It is in 3 dimensions.
% %The "true" trajectory will be generated using 2000 steps.
% numSteps=2000;
% deltaT=tEnd/numSteps;
% xTrue=RungeKSteps(xInit,t0,aDyn,deltaT,numSteps);
% t=t0:deltaT:tEnd;%Times of the refined estimates
% %The true state will have the acceleration from the model in it (which we
% %are appending here). However, we will not give this information to the
% %interpolation algorithm. This is just recorded to demonstrate how well
% %the interpolation algorithm can interpolate moments higher than the
% %given data.
% xTrueAug=zeros(9,numSteps+1);
% xTrueAug(1:6,:)=xTrue;
% for curStep=1:(numSteps+1)
%     derivVal=aDyn(xTrue(:,curStep),t(curStep));
%     xTrueAug(7:9,curStep)=derivVal(4:6);
% end
% xTrue=xTrueAug;
% 
% %The samples will be taken every 25 steps.
% sampIdx=1:25:(numSteps+1);
% xSamp=xTrue(1:6,sampIdx);
% tSamp=t(sampIdx);
% 
% %Generate the interpolating values from the sparsely sampled trajectory.
% groupSize=4;
% [a,c,xRange]=findHermiteInterpPolySet(tSamp,xSamp,numDims,groupSize);
% 
% %Interpolate the state up to the scale factor.
% numDerivs2Interp=2;%Get coefficients for position and velocity.
% xInterp=HermiteSetInterpVal(t,a,c,xRange,numDerivs2Interp);
% 
% %Display the position estimates
% figure(1)
% clf
% hold on
% plot(xTrue(1,:),xTrue(2,:),'-r','linewidth',6)
% plot(xInterp(1,:),xInterp(2,:),'--g','linewidth',3)
% scatter(xSamp(1,:),xSamp(2,:),100,'xb','linewidth',2)
% axis([0 s 0 8e4])
% title('Position Estimates')
% h1=xlabel('x');
% h2=ylabel('y');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
% legend('True Track','Interpolated Track','Samples for Interpolation')
% 
% interpErr=sqrt(sum((xTrue(1:3,:)-xInterp(1:3,:)).^2,1));
% figure(2)
% clf
% hold on
% plot(t,interpErr,'-b','linewidth',3)
% axis([0 tEnd 0 1e-2])
% title('Position Estimate Error')
% h1=xlabel('t');
% h2=ylabel('Error Magnitude');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
% 
% interpErr=sqrt(sum((xTrue(4:6,:)-xInterp(4:6,:)).^2,1));
% figure(3)
% clf
% hold on
% plot(t,interpErr,'-b','linewidth',3)
% axis([0 tEnd 0 1e-2])
% title('Velocity Estimate Error')
% h1=xlabel('t');
% h2=ylabel('Error Magnitude');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
% 
% interpErr=sqrt(sum((xTrue(7:9,:)-xInterp(7:9,:)).^2,1));
% figure(4)
% clf
% hold on
% plot(t,interpErr,'-b','linewidth',3)
% axis([0 tEnd 0 0.08])
% title('Purely Interpolated Acceleration Estimate Error')
% h1=xlabel('t');
% h2=ylabel('Error Magnitude');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%May 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(numDims))
    numDims=1; 
end

if(nargin<4||isempty(groupSize))
    groupSize=2; 
end

Np=length(x);
stateSize=size(y,1);

numMatch=stateSize/numDims;

deltaStart=floor((groupSize-1)/2);
deltaEnd=max(1,ceil((groupSize-1)/2));

xRangeMin=deltaStart+1;
xRangeMax=Np-deltaEnd+1;

xIdxRange=xRangeMin:xRangeMax;

numCoeff=numMatch*(groupSize);

xRange=x(xIdxRange);

numInterpRegions=xRangeMax-xRangeMin;

a=zeros(numDims,numCoeff,numInterpRegions);
c=zeros(numDims,numCoeff-1,numInterpRegions);

for curRegion=1:numInterpRegions
    %This is for interpolation centered about the interval
    %xRange(curRegion) to xRange(curRegion+1).    
    idx=(xIdxRange(curRegion)-deltaStart):(xIdxRange(curRegion)+deltaEnd);
    xCur=x(idx);
    yCur=y(:,idx);
    
    [a(:,:,curRegion),c(:,:,curRegion)]=HermiteInterpMultiPoly(xCur,yCur,numDims);
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
