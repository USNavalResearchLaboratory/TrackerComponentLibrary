function [xList,uList]=RungeKCurvedAtTimes(xInit,uInit,times,aDyn,uFunc,deltaTMax,order,solutionChoice)
%%RUNGEKCURVEDATTIMES Perform multiple steps of Runge-Kutta propagation
%                     with a possibly variable time-interval between steps
%                     given that the axes of the local coordinate system
%                     change as the target moves. The position components
%                     of the state are kept in the global coordinate system
%                     whereas the rest of the components are kept in the
%                     local coordinate system. This Runge-Kutta method is
%                     used to integrate flat-Earth models on a curved
%                     Earth.
%
%INPUTS: xInit The xDimX1 initial value of the state over which integration
%              is being performed. The first three components of xInit must
%              be the position in global Cartesian coordinates. The next 3
%              components must be the target velocity in local (flat)
%              Cartesian coordinates. The other components are arbitrary in
%              the local coordinate system.
%        uInit If the orthonormal basis vectors evolve according to a
%              differential equation, as is the case when moving along
%              geodesics on the surface of the Earth, then this is a 3X3
%              matrix of vectors specifying the local coordinate axes.
%              uInit(:,i) corresponds to the ith position component in
%              xInit. That is xInit(i) for i from 1 to 3. On the other
%              hand, if the basis vectors are deterministically known at
%              all locations as a function of x and t, then an empty matrix
%              should be passed for uInit.
%        times The times at which state estimates are desired. times(1) is
%              the time of xInit.
%         aDyn aDyn(xVal,curT) returns the derivative of the state taken at
%              time curT. If this parameter is omitted, then
%              aDyn=@(x,t)aPoly(x,3); will be used. This assumes that the
%              state is in 3D and a simple, position, velocity, etc. motion
%              model is used.
%        uFunc If uInit is not empty, then uFunc is a function handle of
%              the form uDyn(u,x,t) that returns the derivatives of the
%              basis vectors, where u is the current set of basis vectors,
%              as is the case when traveling a geodesic curve. If uInit is
%              an empty matrix, then uFunc should be a function handle of
%              the form uBasis(x,t), where x is the state and t is the
%              time, that returns the local basis vectors in the global
%              coordinate system. If omitted, and uInit is not empty, then
%              it is assumed that one is traveling geodesic-style
%              trajectories on the WGS-84 reference ellipsoid and the
%              function uDotEllipsoid is used as the derivative function.
%              If omitted and uInit is empty, then it is assumed that one
%              is traveling constant-heading trajectories (rhumb-lines) on
%              the WGS-84 reference ellipsoid and the function
%              @(x,t)getENUAxes(ellips2Cart(x(1:3))) is used.
%    deltaTMax The maximum allowable step size for propagation of the
%              state. If this parameter is omitted, then Inf will be used,
%              meaning that the steps are between the individual times.
%        order The order of the Runge-Kutta method. If this parameter is
%              omitted, then the default order of 4 is used. Order can
%              range from 1 to 7.
% solutionChoice When multiple formulae are implemented, this selects which
%              one to use. Otherwise, this parameter is not used.
%              Currently, only the order=4 method has multiple solutions
%              implemented in which case omitting this parameter or setting
%              it to zero used the Dormand and Prince Algorithm, and
%              setting it to 1 uses the Fehlberg algorithm.
%
%OUTPUTS: xList The target state at the times given in times.
%         uList This is the 3*3*numTimes local coordinate system axes at
%               the times given in times.
%
%This function maps an arbitrary continuous-time, deterministic flat-Earth
%dynamic model to curvature of the WGS-84 ellipsoid. The algorithm is
%described in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(solutionChoice))
    solutionChoice=0;
end

if(nargin<7||isempty(order))
    order=4;
end

if(nargin<6||isempty(deltaTMax))
    deltaTMax=Inf;
end

if(nargin<5||isempty(uFunc))
    %If the basis vectors are deterministic
    if(isempty(uInit))
        uFunc=@(x,t)getENUAxes(Cart2Ellipse(x(1:3)));
    else%If the basis vectors evolve.
        uFunc=@(u,x,t)uDotEllipsoid(u,x);
    end
end

if(nargin<4||isempty(aDyn))
    aDyn=@(x,t)aPoly(x,3);
end

xDim=length(xInit);
numTimes=length(times);

if(isempty(uInit))
    %If the basis vectors are deterministically given in all locations.
    aDynGlob=@(x,t)aDynFlat2CurvedUDet(x,t,aDyn,uFunc);
    xList=RungeKAtTimes(xInit,times,aDynGlob,deltaTMax,order,solutionChoice);
    
    if(nargout>1)
        uList=zeros(3,3,numTimes);
        for curTime=1:numTimes
            uList(:,:,curTime)=uFunc(xList(:,curTime),times(curTime));
        end
    end
else
    %We have to perform integration over a stacked state and basis vector
    %set.
    xStacked=[xInit;uInit(:)];
    aDynGlob=@(x,t)aDynFlat2CurvedUDyn(x,t,aDyn,uFunc);

    xStackedList=RungeKAtTimes(xStacked,times,aDynGlob,deltaTMax,order,solutionChoice);

    xList=reshape(xStackedList(1:xDim,:),xDim,numTimes);
    uList=reshape(xStackedList((xDim+1):end,:),3,3,numTimes);
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
