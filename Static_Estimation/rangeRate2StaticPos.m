function [zCart,exitCode]=rangeRate2StaticPos(rr,sRx,AbsTol,RelTol,params)
%%RANGERATE2STATICPOS Given the minimum number of range rate measurements
%          for observability (2 in 2D and 3 in 3D) determine the location
%          of a stationary emitter. The receivers must be moving for the
%          problem to be observable. For example, the receiver could be a
%          drone that is flying around and the emitter a cell phone (which
%          would have a known broadcast frequency) held by a stationary
%          person talking on it.
%
%INPUTS: rr A 2X1 (for 2D) or 3X1 (for 3D) vector of range rate
%           measurements. Units are typically meters per second.
%       sRx A 4X2 or 6X3 matrix of the stacked position and velocity of the
%           receiver for each measurement. For example, in 3D
%           lRx(:,i)=[x;y;z;xDot;yDot;zDot]. For observability, none of the
%           receivers can be stationary.
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
%           absDiff<RelTol*abs(TDOA(curEq))) is true. This criterion is
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
%OUTPUTS: zCart The 2XnumSol or 3XnumSol set of 2D or 3D real solutions for
%               the location of the target.
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
%EXAMPLE 1:
%This is a 2D exmaple.
% uTrue=[1e3;5e3];
% s=[500, 1100;
%     2500, 2500]; 
% sDot=[300, 300;
%       0, 0];
% stateRx=[s;sDot];
% %Generate measurements of observed frequencies
% rr=zeros(2,1);
% for curMeas=1:2
%    rr(curMeas)=-sDot(:,curMeas)'*(uTrue-s(:,curMeas))/norm(uTrue-s(:,curMeas));
% end
% zCart=rangeRate2StaticPos(rr,stateRx)
%One will obtain the correct solution plus one other solution due to
%ambiguity.
%
%EXAMPLE in 3D:
%This is a 3D example.
% uTrue=[1e3;5e3;2e3];
% s=[500, 1100,5000;
%    3000, 1000,-1000;
%    0,    2000,1000];
% sDot=[300, 100,0;
%       0,   100,300;
%       0,   0,  0];
% stateRx=[s;sDot];
% %Generate measurements of observed frequencies
% rr=zeros(3,1);
% for curMeas=1:3
%    rr(curMeas)=-sDot(:,curMeas)'*(uTrue-s(:,curMeas))/norm(uTrue-s(:,curMeas));
% end
% zCart=rangeRate2StaticPos(rr,stateRx)
% %One will obtain the correct solution plus three other solutions (due to
% %ambiguity).
%
%REFERENCES:
%[1] D. F. Crouse, "General Multivariate Polynomial Target Localization and
%    Initial Estimation," Journal of Advances in Information Fusion, vol.
%    13, no. 1, pp. 68-91, Jun. 2018.
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(AbsTol))
    AbsTol=1e-9;
end

if(nargin<4||isempty(RelTol))
   RelTol=1e-7; 
end

if(nargin<5)
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

numDim=size(sRx,1)/2;

lList=sRx(1:numDim,:);
lDotList=sRx((numDim+1):end,:);

switch(numDim)
    case 2%2D
        polyCoeffMats=cell(2,1);
        for curMeas=1:2
            l=lList(:,curMeas);
            lDot=lDotList(:,curMeas);
            rDot=rr(curMeas);

            lTilde=2*(lDot*(l'*lDot)-rDot^2*l);
            cTilde=rDot^2*norm(l)^2-(l'*lDot)^2;

            coeffs=zeros(3,3);
            coeffs(2+1,0+1)=(rDot^2-lDot(1)^2);
            coeffs(0+1,2+1)=(rDot^2-lDot(2)^2);
            
            coeffs(1+1,1+1)=-2*lDot(1)*lDot(2);
            
            coeffs(1+1,0+1)=lTilde(1);
            coeffs(0+1,1+1)=lTilde(2);
            
            coeffs(0+1,0+1)=cTilde;
            polyCoeffMats{curMeas}=coeffs;
        end
    case 3%3D
        polyCoeffMats=cell(3,1);
        for curMeas=1:3
            l=lList(:,curMeas);
            lDot=lDotList(:,curMeas);
            rDot=rr(curMeas);

            lTilde=2*(lDot*(l'*lDot)-rDot^2*l);
            cTilde=rDot^2*norm(l)^2-(l'*lDot)^2;

            coeffs=zeros(3,3,3);
            coeffs(2+1,0+1,0+1)=(rDot^2-lDot(1)^2);
            coeffs(0+1,2+1,0+1)=(rDot^2-lDot(2)^2);
            coeffs(0+1,0+1,2+1)=(rDot^2-lDot(3)^2);
            
            coeffs(1+1,1+1,0+1)=-2*lDot(1)*lDot(2);
            coeffs(1+1,0+1,1+1)=-2*lDot(1)*lDot(3);
            coeffs(0+1,1+1,1+1)=-2*lDot(2)*lDot(3);
            
            coeffs(1+1,0+1,0+1)=lTilde(1);
            coeffs(0+1,1+1,0+1)=lTilde(2);
            coeffs(0+1,0+1,1+1)=lTilde(3);
            
            coeffs(0+1,0+1,0+1)=cTilde;
            polyCoeffMats{curMeas}=coeffs/1e3;
        end
    otherwise
        error('Invalid Dimensionality')
end
[theRoots,exitCode]=polyRootsMultiDim(polyCoeffMats,maxDegIncreases,useMotzkinNull);

%Eliminate complex solutions.
sel=sum((abs(imag(theRoots))<AbsTol) | abs(imag(theRoots))<RelTol*abs(real(theRoots)),1)~=0;
zCart=real(theRoots(:,sel));

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
    for curEq=1:numDim
        diff=zCart(:,curSol)-lList(:,curEq);
        rrComp=-lDotList(:,curEq)'*diff/norm(diff);

        absDiff=abs(rrComp-rr(curEq));

        if(~(absDiff<AbsTol || absDiff<RelTol*abs(rr(curEq))))
            sel(curSol)=false;
            break;
        end
    end
end
zCart=zCart(:,sel);

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
