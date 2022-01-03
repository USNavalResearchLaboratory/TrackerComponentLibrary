function [zCart,exitCode]=rangeRateRatio2StaticPos2D(fRat,sRRef,sRx,c,AbsTol,RelTol,params)
%%RANGERATERATIO2STATICPOS2D It is possible to localize a stationary
%               emitter using only Doppler measurements from moving sensors
%               without knowing the transmission frequency of the emitter.
%               In such an instance, one would measure the (baseband or
%               passband) frequency of the emitter at multiple location,
%               and then from the ratio of the frequencies, one can
%               determine the location of the emitter. This function
%               localizes an emitter in 2D using only the ratios of
%               frequency measurements.
%
%INPUTS:fRat A 2X1 vector of ratios of frequency measurements from a
%            statonary emitter. In both instances, the numerator is the
%            frequency measurement from the reference sensor and the
%            denominator is the frequency measurement of a different
%            sensor. This means, there must be a total of three sensors.
%      lRRef The 4X1 state (position and velocity) of the reference
%            sensor. This consists of components [x,y,xDot,yDot].
%        lRx The 4X2 states (2D position and velocity, stacked) of the
%            two other sensors. lRx(:,i) is the state corresponding to
%            the reference sensors in the measurement rrRat(i).
%          c The speed of propagation of whatever signal is received. If
%            this parameter is omitted or an empty matrix is passed, a
%            default value of Constants.speedOfLight is used.
% AbsTol,RelTol In some geometries of targets and recievers and when noise
%           is added to the measurements, meaningless complex solutions
%           can be produced. AbsTol and RelTol are tolerances for
%           determining whether is solution is numerically real and
%           should be kept or numerically complex and should be
%           discarded. A solution x is kept if abs(imag(x))<AbsTol
%           or if abs(imag(x))<RelTol*abs(real(x)). If omitted or empty
%           matrices are passed, the default values of AbsTol=1e-9 and
%           RelTol=1e-7 are used. Additionally, AbsTol and RelTol are
%           used to eliminate false solutions that arise due to the range
%           rates being squared when formulating the polynomial problem.
%           That is, some solutions can arise having the wrong sign. These
%           solutions are eliminated if ~(absDiff<AbsTol || 
%           absDiff<RelTol*abs(fRat(curEq))) is true. This criterion is
%           used rather than just comparing the signs of range rate
%           computations to improve performance when range rate values are
%           small/ zero.
%    params An optional structure of parameters for the polyRootsMultiDim
%           function, which is used by this function. possible members are
%           maxDegIncreases and useMotzkinNull, both of which are described
%           in the function polyRootsMultiDim. Default values for the
%           function polyRootsMultiDim will be used if this parameter is
%           omitted or an empty matrix is passed.
%
%OUTPUTS: zCart The 2XnumSol set of 2D real solutions for the location of
%               the target.
%      exitCode The exit code returned by the function polyRootsMultiDim,
%               which is used to solve the problem.
%
%This function implements some of the concepts described in [1], but
%without making use of external simultaneous multivariate polynomial
%solvers. 
%
%The problem is formulated as simultaneous multivariate polynomials and is
%solved using the polyRootsMultiDim function.
%
%%EXAMPLE:
%A 2D example
% uTrue=[1e3;5e3];
% 
% lRx1=[1000;3000];
% lRx1Dot=[150;-150];
% lRx=[500, 1100;
%     2500, 2500]; 
% lRxDot=[300, 300;
%       0, 0];
% sRRef=[lRx1;lRx1Dot];
% sRx=[lRx;lRxDot];
% 
% rr1=-lRx1Dot'*(uTrue-lRx1)/norm(uTrue-lRx1);
% rr=zeros(2,1);
% for curMeas=1:2
%    rr(curMeas)=-lRxDot(:,curMeas)'*(uTrue-lRx(:,curMeas))/norm(uTrue-lRx(:,curMeas));
% end
% c=Constants.speedOfLight;
% fRat=(1-rr1/c)./(1-rr/c);
% [zCart,exitCode]=rateRateRatio2StaticPos2D(fRat,sRRef,sRx)
% %One will get the correct solution plus three other ones due to
% %ambiguity.
%
%REFERENCES:
%[1] D. F. Crouse, "General Multivariate Polynomial Target Localization and
%    Initial Estimation," Journal of Advances in Information Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(c))
   c=Constants.speedOfLight;
end

if(nargin<5||isempty(AbsTol))
    AbsTol=1e-9;
end

if(nargin<6||isempty(RelTol))
   RelTol=1e-7; 
end

if(nargin<7)
   params=[]; 
end

maxDegIncreases=[];
useMotzkinNull=[];
if(~isempty(params))
    if(isfield(params,'maxDegIncreases'))
        maxDegIncreases=params.maxDegIncreases; 
    end

    if(isfield(params,'useMotzkinNull'))
        useMotzkinNull=params.useMotzkinNull;
    end
end

lRx1=sRRef(1:2);
lRx1Dot=sRRef(3:4);
lRx=sRx(1:2,:);
lRxDot=sRx(3:4,:);

%Variables are ordered [tx,ty,r1,r2,r3];
polyCoeffMats=cell(5,1);

%%Add the equation for the r1 term.
l1x=lRx1(1);
l1y=lRx1(2);

coeffs=zeros(3,3,3,3,3);
coeffs(0+1,0+1,2+1,0+1,0+1)=1;
coeffs(2+1,0+1,0+1,0+1,0+1)=-1;
coeffs(0+1,2+1,0+1,0+1,0+1)=-1;
coeffs(1+1,0+1,0+1,0+1,0+1)=2*l1x;
coeffs(0+1,1+1,0+1,0+1,0+1)=2*l1y;
coeffs(0+1,0+1,0+1,0+1,0+1)=-l1x^2-l1y^2;
polyCoeffMats{1}=coeffs;

%Add the equations for the r2 and r3 terms.
for curR=2:3
    ljx=lRx(1,curR-1);
    ljy=lRx(2,curR-1);
    coeffs=zeros(3,3,3,3,3);
    if(curR==2)
        coeffs(0+1,0+1,0+1,2+1,0+1)=1;
    else
        coeffs(0+1,0+1,0+1,0+1,2+1)=1;
    end

    coeffs(2+1,0+1,0+1,0+1,0+1)=-1;
    coeffs(0+1,2+1,0+1,0+1,0+1)=-1;
    coeffs(1+1,0+1,0+1,0+1,0+1)=2*ljx;
    coeffs(0+1,1+1,0+1,0+1,0+1)=2*ljy;
    coeffs(0+1,0+1,0+1,0+1,0+1)=-ljx^2-ljy^2;

    polyCoeffMats{curR}=coeffs;
end

%Now, we set the equations based on the measurements.
l1xDot=lRx1Dot(1);
l1yDot=lRx1Dot(2);
for curR=1:2
    lj=lRx(:,curR);
    ljDot=lRxDot(:,curR);
    ljxDot=lRxDot(1,curR);
    ljyDot=lRxDot(2,curR);

    f1j=fRat(curR);

    coeffs=zeros(2,2,2,2,2);

    %The coefficient for the r1 term.
    coeffs(0+1,0+1,1+1,0+1,0+1)=-f1j*dot(lj,ljDot);

    %The coefficient for the rj term.
    rjCoeff=dot(lRx1,lRx1Dot);
    if(curR==1)
        coeffs(0+1,0+1,0+1,1+1,0+1)=rjCoeff;
    else
        coeffs(0+1,0+1,0+1,0+1,1+1)=rjCoeff;
    end

    %r1*rj coefficient
    r1rjCoeff=c*(f1j-1);
    if(curR==1)
        coeffs(0+1,0+1,1+1,1+1,0+1)=r1rjCoeff;
    else
        coeffs(0+1,0+1,1+1,0+1,1+1)=r1rjCoeff;
    end

    %r1*tx coefficient
    coeffs(1+1,0+1,1+1,0+1,0+1)=f1j*ljxDot;

    %r1*ty coefficient
    coeffs(0+1,1+1,1+1,0+1,0+1)=f1j*ljyDot;

    %rj*tx coefficient
    if(curR==1)
        coeffs(1+1,0+1,0+1,1+1,0+1)=-l1xDot;
    else
        coeffs(1+1,0+1,0+1,0+1,1+1)=-l1xDot;
    end

    %rj*ty coefficient
    if(curR==1)
        coeffs(0+1,1+1,0+1,1+1,0+1)=-l1yDot;
    else
        coeffs(0+1,1+1,0+1,0+1,1+1)=-l1yDot;
    end

    polyCoeffMats{3+curR}=coeffs;
end

[theRoots,exitCode]=polyRootsMultiDim(polyCoeffMats,maxDegIncreases,useMotzkinNull);

%Eliminate complex solutions.
sel=sum((abs(imag(theRoots))<AbsTol) | abs(imag(theRoots))<RelTol*abs(real(theRoots)),1)~=0;
zCart=real(theRoots(:,sel));
zCart=zCart(1:2,:);

%Now, due to squaring in the equations, the sign of the solutions will
%generally be incorrect. Thus, we shall discard any solutions where the
%sign of the computed range rate measurements does not match the sign of
%the actual range rate measurements. However, this poses problems for range
%rate values near zero. Thus, as opposed to comparing the sign, we use the
%same AbsTol and RelTol as for determining whether something is complex to
%determine whether we should discard a solution.
numSol=size(zCart,2);
sel=true(numSol,1);
for curSol=1:numSol
    diff=zCart(:,curSol)-lRx1;
    rrRef=-lRx1Dot'*diff/norm(diff);
    
    for curEq=1:2
        diff=zCart(:,curSol)-lRx(:,curEq);
        rrCur=-lRxDot(:,curEq)'*diff/norm(diff);
        
        fRatCur=(1-rrRef/c)/(1-rrCur/c);
        
        absDiff=abs(fRatCur-fRat(curEq));

        if(~(absDiff<AbsTol || absDiff<RelTol*abs(fRat(curEq))))
            sel(curSol)=false;
            break;
        end
    end
end
zCart=zCart(:,sel);

end