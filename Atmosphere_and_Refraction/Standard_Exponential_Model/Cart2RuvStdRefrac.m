function [z,uTx,uTarRx,uTarTx]=Cart2RuvStdRefrac(zC,useHalfRange,zTx,zRx,M,Ns,includeW,ce,rE,spherCent)
%%CART2RUVSTDREFRAC Convert points in Cartesian coordinates (either the
%         global system or a local system at the receiver) into local
%         bistatic r-u-v coordinates of the receiver, accounting for how a
%         standard exponential atmospheric model warps the measurements.
%         r-u-v coordinates consist of a bistatic range and direction
%         cosines at the receiver. The "direction cosines" u and v are just
%         the x and y coordinates of a unit vector from the receiver to
%         the target in the coordinate system at the receiver (Here, to the
%         refraction-corrupted apparent direction). This basically assumes
%         that the boresight direction of the receiver is the z axis.
%         Assuming the target is in front of the receiver, the third unit
%         vector coordinate is not needed. However, with the includeW
%         option, it can be provided, resulting in r-u-v-w coordinates. This
%         function is not suitable for computing refraction between
%         satellites, grazing the Earth's atmosphere. The algorithm might
%         have an error if the ray path goes too far underground.
%
%INPUT: zC  A 3XN matrix of Cartesian points (target locations) in global
%           [x;y;z] Cartesian coordinates.
% useHalfRange A boolean value specifying whether the bistatic range value
%           should be divided by two. This normally comes up when operating
%           in monostatic mode, so that the range reported is a one-way
%           range. The default if this parameter is not provided (or an
%           empty matrix is provided) is false.
%       zTx The 3X1 [x;y;z] location vector of the transmitter in global
%           Cartesian coordinates.  If this parameter is omitted or an
%           empty matrix is passed, then the receiver is placed at the origin.
%       zRx The 3X1 [x;y;z] location vector of the receiver in global
%           Cartesian coordinates. If this parameter is omitted or an empty
%           matrix is passed, then the receiver is placed at the origin.
%         M A 3X3 rotation matrix to go from the alignment of the global
%           coordinate system to the local alignment of the receiver. The z
%           vector of the local coordinate system of the receiver is the
%           pointing direction of the receiver. If this matrix is omitted,
%           then the identity matrix is used.
%        Ns The atmospheric refractivity reduced to the reference sphere.
%           Note that the refractivity is (n-1)*1e6, where n is the index
%           of refraction. The function reduceStdRefrac2Spher can be used
%           to reduce a refractivity to the surface of a reference
%           ellipsoid. This function does not allow different
%           refractivities to be used as the transmitter and receiver. If
%           this parameter is omitted or an empty matrix is passed, a
%           default value of 313 is used.
%  includeW An optional boolean value indicating whether a third direction
%           cosine component should be included. The u and v direction
%           cosines are two parts of a 3D unit vector. Generally, one might
%           assume that the target is in front of the sensor, so the third
%           component would be positive and is not needed. However, the
%           third component can be included if ambiguity exists. The
%           default if this parameter is omitted or an empty matrix is
%           passed is false.
%        ce The optional decay constant of the exponential model. The
%           refractivity N at height h is N=Ns*exp(-ce*(h-h0)) where h0 is
%           the reference height (in this function, the height of the
%           reference ellipsoid surface is used). ce is related to the
%           change in refractivity at an elevation of 1km based on the
%           refractivity at sea level as
%           ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
%           where the change in refractivity for a change in elevation of
%           1km is DeltaN=-multConst*exp(expConst*Ns); In [1], standard
%           values for the two constants are expConst=0.005577; and
%           multConst=7.32; If ce is omitted or an empty matrix is passed,
%           the value based on the standard model is used.
% rE,spherCent The radius of the Earth to use for the spherical Earth
%           approximation used in the model and also the offset between the
%           global model and the local spherical model. It is assumed that
%           zC,zTx,and zRx are all given in the global model and will need
%           to be transformed to the local model to the used. If rE is
%           omitted or an empty matrix is passed, then the default of
%           [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx)) is
%           used. The defaults here mean that a WGS-84 reference ellipsoid
%           is approximated by the local osculating sphere.
%
%OUTPUTS: z The 3XN (or 4XN if includeW is true) matrix of location vectors
%           of the points in bistatic [r;u;v] coordinates. If
%           useHalfRange=true, then the r component is half the bistatic
%           range (half the round-trip range for a monostatic scenario).
%       uTx A 3XN set of unit vectors pointing from the transmitter to the
%           refraction-corrupted position (of the target as seen by the
%           transmitter. This is in the global coordinate system.
%    uTarRx A 3XN set of unit vectors pointing from the target to the
%           refraction-corrupted position of the receiver as seen by the
%           target. This is in the global coordinate system.
%    uTarTx A 3XN set of unit vectors pointing from the target to the
%           refraction-corrupted position of the transmitter as seen by the
%           target. This is in the global coordinate system.
%
%This function implements the refraction algorithm for the basic
%exponential atmosphere as described in [1] for the bistatic case. If the
%target is collocated with the transmitter or the receiver, then NaNs will
%be returned for some values. The basic exponential refraction model is in
%[2].
%
%The model is parameterized in terms of a height above a sphere. The Earth
%is more of an ellipsoid than a sphere. Thus, we use local spherical
%approximations about the transmitter and the receiver. That is, for
%computing the refraction from the transmitter to the target, we use the
%distance from the center of the Earth to the surface of the reference
%ellipsoid at the transmitter as the radius of an approximately spherical
%Earth. Similarly, the distance from the center of the Earth to the
%receiver is used in the approximation for the path from the target to the
%receiver.
%
%The algorithm in [1] performs ray tracing by solving a boundary value
%problem. here, the bvp5c function in Matlab is used to solve the problem.
%
%For paths that are nearly vertical, it is approximated that there is no
%bending in angle and an explicit solution to the integral over the index
%of refraction in the vertical direction is used to obtain the range.
%
%EXAMPLE:
%Here, we have two radars and one target near Hawaii.
% latLonRx=deg2rad([20.269202;-155.852051]);
% AltRx=0;
% latLonTx=deg2rad([20.724568;-155.978394]);
% AltTx=0;
% latLonTar=deg2rad([20.835390;-155.313721]);
% AltTar=8e3;%8km target altitude.
% %Convert locations to Cartesian.
% zRx=ellips2Cart([latLonRx;AltRx]);
% zTx=ellips2Cart([latLonTx;AltTx]);
% zTar=ellips2Cart([latLonTar;AltTar]);
% 
% %The receiver faces 45 degrees East of North and 15 degrees up from the
% %local ellipsoidal level.
% M=findRFTransParam([latLonRx;AltRx],deg2rad(45),deg2rad(15));
% Ns=350;%Assumed refractivity at the sea surface.
% useHalfRange=false;
% includeW=true;%Include third dimension of unit vector.
% [z,uTx,uTarRx,uTarTx]=Cart2RuvStdRefrac(zTar,useHalfRange,zTx,zRx,M,Ns,includeW);
% zNoRefrac=Cart2Ruv(zTar,useHalfRange,zTx,zRx,M,includeW);
% z(1)-zNoRefrac(1)%Bistatic range difference of 31.0813 meters
% %Direction difference of 0.0927 degrees
% rad2deg(angBetweenVecs(z(2:end),zNoRefrac(2:end)))
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%[2] B. R. Bean and G. D. Thayer, CRPL Exponential Reference Atmosphere.
%    Washington, D.C.: U. S. Department of Commerce, National Bureau of
%    Standards, Oct. 1959. [Online]. Available:
%    http://digicoll.manoa.hawaii.edu/techreports/PDF/NBS4.pdf
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(zC,2);

%The assumes refractivity at sea level to use if none is provided.
if(nargin<7||isempty(includeW))
    includeW=false; 
end

if(nargin<6||isempty(Ns))
    Ns=313;
end

if(nargin<5||isempty(M))
    M=eye(3); 
end

if(nargin<4||isempty(zRx))
    zRx=zeros(3,1); 
end

if(nargin<3||isempty(zTx))
    zTx=zeros(3,1);
end

if(nargin<2||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<8||isempty(ce))
    expConst=0.005577;
    multConst=7.32;

    %The change in refractivity at an elevation of 1km based on the
    %refractivity on the surface of the Earth.
    DeltaN=-multConst*exp(expConst*Ns);
    ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
end

if(nargin<9||isempty(rE))
    %Use the radius of the Earth that is the radius of the osculating
    %sphere at the location of the observed. This will be the radius used
    %in the local spherical Earth approximation for computing atmospheric
    %refraction. This uses the WGS-84 reference ellipsoid.
    [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx));
end

%Adjust all the Cartesian values based on the osculating sphere.
zC=zC-spherCent;
zTx=zTx-spherCent;
zRx=zRx-spherCent;

%Allocate space
if(includeW)
    z=zeros(4,numMeas);
else
    z=zeros(3,numMeas);
end
uTx=zeros(3,numMeas);
uTarTx=zeros(3,numMeas);
uTarRx=zeros(3,numMeas);

if(any(zRx~=zTx))%If the scenario is bistatic
    for curMeas=1:numMeas
        [r2,uArrive,uTarRx(:,curMeas)]=atmosRefracMeas(zRx,zC(:,curMeas),Ns,ce,rE);
        [r1,uTx(:,curMeas),uTarTx(:,curMeas)]=atmosRefracMeas(zTx,zC(:,curMeas),Ns,ce,rE);

        r=(r1+r2);
        u=M*uArrive;
        
        if(useHalfRange)
            r=r/2; 
        end
        
        zCur=[r;u];
    
        if(includeW)
            z(:,curMeas)=zCur;
        else
            z(:,curMeas)=zCur(1:3);
        end
    end
else%The scenario is monostatic.
    for curMeas=1:numMeas
        [range,uArrive,uTarRx(:,curMeas)]=atmosRefracMeas(zTx,zC(:,curMeas),Ns,ce,rE);
        uTx=uArrive;
        r=2*range;%Round-trip range.
        u=M*uArrive;

        if(useHalfRange)
            r=r/2; 
        end

        zCur=[r;u];

        if(includeW)
            z(:,curMeas)=zCur;
        else
            z(:,curMeas)=zCur(1:3);
        end

        uTx(:,curMeas)=uArrive;
        uTarTx(:,curMeas)=uTarRx(:,curMeas);
    end
end
end

function [range,uArrive,uDepart]=atmosRefracMeas(xObs,xObj,Ns,ce,rE)
%%ATMOSREFRACMEAS  Given the location of an observer and an object in the
%                  atmosphere of the Earth, find the delay and angle of
%                  arrival of a signal from the object to the observer,
%                  accounting for basic,standard  refraction. A
%                  low-fidelity exponential atmospheric model is used. This
%                  function is not suitable for computing refraction
%                  between satellites, grazing the Earth's atmosphere. The
%                  algorithm might have an error if the raypath
%                  goes too far underground.
%
%INPUTS: xObs The Cartesian location of the observer in ECEF coordinates in
%             meters as [x;y;z].
%        xObj The Cartesian location of the object being observed in
%             ECEF coordinates in meters as [x;y;z].
%          Ns The refractivity as reduced to at sea level.
% expConst,multConst The parameters of the refractivity model such that
%             increasing the height by 1km, the model is that the
%             refractivity changes by deltaN=-multConst*exp(expConst*N).
%             deltaN cannot be negative.
%          rE The radius of the Earth to use in the spherical Earth
%             approximation. An osculating sphere near a sensor is
%             suggested.
%
%OUTPUTS: range  The apparent one-way range (in meters) of the signal from
%                the transmitter to the target and back to the receiver.
%        uArrive A unit vector in ECEF coordinates pointing in the apparent
%                direction of the signal the observer received (as seen by
%                the observer).
%        uDepart The direction of the signal departing the object that
%                arrives at the observer. Put another way, if the observer
%                were to transmit a signal to the object, this is the
%                apparent direction of the observer as seen by the object.
%
%This function implements the refraction algorithm for the basic
%exponential atmosphere as described in [1] for the monostatic case. If the
%points are collocated, then NaNs are returned for uArrive and uDepart.
%
%The function will fail for paths that go too deep into the Earth.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%We need the conversion from the 3D coordinate system of the observer and
%object into the 2D coordinate system used for raytracing. The 2D
%coordinate system has the center of the Earth as its origin and the x-y
%axes are in the plane of the vector from the observer to the target. One
%vector common to both coordinate systems in the local up vector, which
%will be the local y axis. The second vector common to both will be the
%local x vector, which will be the projection of xObj-xObs onto the local
%tangent plane. Here, the vertical is the spherical model vertical. Since
%the precision of the model is low enough that the difference between the
%spherical and gravitational verticals shouldn't matter.
uENU=getENUAxes(Cart2Ellipse(xObs,[],rE,0));
uVertGlobal=uENU(:,3);
uVertLocal=[0;1;0];

vec2TarGlobal=xObj-xObs;

%The projection of the xObj-xObs vector into the local tangent plane can be
%obtained by subtracting the component of the vector that is orthogonal to
%the plane.
uHorizGlobalOrig=vec2TarGlobal-dot(vec2TarGlobal,uENU(:,3))*uENU(:,3);

uHorizGlobal=uHorizGlobalOrig/norm(uHorizGlobalOrig);
uHorizLocal=[1;0;0];

%Find the rotation matrix from the global coordinate system into the local
%coordinate system.
ECEF2LocalRot=findTransParam([uVertLocal,uHorizLocal],[uVertGlobal,uHorizGlobal]);

%The third (z) coordinate in the local system should be zero after this
%transformation.
vec2TarLocal=ECEF2LocalRot*vec2TarGlobal;

%The location of the receiver in the local 2D coordinate system.
x0Init=0;
y0Init=norm(xObs);

%The location of the target in the local 2D coordinate system.
x1Init=vec2TarLocal(1);
y1Init=vec2TarLocal(2)+y0Init;

%If the two points are nearly vertical, then the ray tracing algorithm will
%fail. For nearly vertical points, the bending due to refraction in the
%model should be negligible, so we can perform an integral in the y
%direction to solve for the excess range instead of having to solve the
%more complicated general bistatic problem. To deal with x1Init not being
%exactly zero, we actually go to a full range of norm([x1Init;y1Init])
if(norm(uHorizGlobalOrig)<1e-3)
    %Integrating the index of refraction from y0Init to
    %norm([x1Init;y1Init]) at a constant x yields the following measured
    %range. It is not just the geometric range, because the index of
    %refraction is not a constant 1.
    yMax=norm([x1Init;y1Init]);
   
    range=((exp(ce*(rE-y0Init))-exp(ce*(rE-yMax)))*Ns)/(1e6*ce)+yMax-y0Init;

    uArrive=vec2TarGlobal/norm(vec2TarGlobal);
    uDepart=-uArrive;
    return;
end

%Now, set up the boundary-value problem to determine the path taken by
%light between the target and the receiver.

%The initial guess is just the linear solution. The solver requires a fixed
%number of steps. 20 is probably sufficient for things near the Earth. that
%is, up to distances of, say 400km. We can scale the number of steps as 20
%for every 400 kilometers with a minimum of, say 10.
%Things outside of the atmosphere should use the astronomical refraction
%routines.
numSteps=max(20,ceil(20*norm(vec2TarLocal)/400e3));
x=linspace(x0Init,x1Init,numSteps);
slope=(y1Init-y0Init)/(x1Init-x0Init);
b=y1Init-slope*x1Init;%The y-intercept.
%The initial estimate of the solution.
solInit=bvpinit(x,@(x)[x*slope+b;slope]);

%Now, solve the differential equation.
oldOpts=bvpset();
newOpts=bvpset(oldOpts,'RelTol',1e-8,'AbsTol',1e-8,'FJacobian',@(x,y)odefunJacob(x,y,Ns,rE,ce),'BCJacobian',@bcfunJacob);%Increase the accuracy.
sol=bvp5c(@(x,y)expDiffEq(x,y,Ns,rE,ce),@(y0,y1)bcfun(y0,y1,y0Init,y1Init),solInit,newOpts);

%Get the refraction-corrupted range measurement for a signal traveling from
%the object to the observer. 
range=integral(@(x)pathFun2D(x,sol,Ns,rE,ce),x0Init,x1Init,'AbsTol',eps(1),'RelTol',1e-15);

if(nargout>1)
    %Get the angle of arrival for a signal traveling from the object to the
    %observer. The angle is determined by the slope at the initial point.
    thetaOrig=atan(sol.y(2,1));
    uLocal=[cos(thetaOrig);sin(thetaOrig);0];

    %The inverse rotation is given the transpose of the rotation
    %matrix. This is the apparent direction of the object as seen by the
    %observer.
    uArrive=ECEF2LocalRot'*uLocal;

    thetaEnd=atan(sol.y(2,end));
    uLocal=[cos(thetaEnd);sin(thetaEnd);0];
    %This is the apparent direction of the observer as seen by the object.
    uDepart=-ECEF2LocalRot'*uLocal;
end
end

function val=pathFun2D(x,sol,Ns,rE,ce)
    %This function is used to integrate the time taken
    y=deval(x,sol);
    val=(1+NRefracExp(x,y(1,:),Ns,rE,ce)).*sqrt(1+y(2,:).^2);
end

function res=bcfun(y0,y1,y0Init,y1Init)
    %The residue to define the boundary condition for the numeric
    %differential equation solver as applied to the 2D exponential
    %atmospheric model.

    res=[y0(1)-y0Init;
         y1(1)-y1Init];
end

function [dbcy0,dbcy1]=bcfunJacob(~,~)
    %The Jacobians of the boundary conditions for raytracing the 2D
    %exponential atmospheric refraction model.
    dbcy0=[1 0
           0 0];
    dbcy1=[0 0
           1 0];
end

function J=odefunJacob(x,y,Ns,rE,ce)
    %The Jacobian of the differential equation for raytracing the 2D
    %exponential atmospheric model.
    expVal=NRefracExp(x,y(1),Ns,rE,ce);

    J=zeros(2,2);
    J(1,2)=1;
    J(2,1)=ce*(1+y(2)^2)*(-expVal)*(ce*y(1)*(x*y(2)-y(1))*sqrt(x^2+y(1)^2)+x*(x+y(1)*y(2))*(expVal+1))/((x^2+y(1)^2)^(3/2)*(expVal+1)^2);
    J(2,2)=ce*(x-2*y(1)*y(2)+3*x*y(2)^2)*expVal/((expVal+1)*sqrt(x^2+y(1)^2));
end

function dxdy=expDiffEq(x,y,Ns,rE,ce)
    %Find the refractivity at location (x,y).
    expVal=NRefracExp(x,y(1),Ns,rE,ce);

    dxdy=[y(2)
          ce*(1+y(2)^2)*(x*y(2)-y(1))*expVal/((expVal+1)*sqrt(x^2+y(1)^2))];
end

function nRefrac=NRefracExp(x,y,Ns,rE,ce)
    %The refractivity. This is 10^6*(index of refraction-1)
    nRefrac=1e-6*Ns*exp(-ce*(sqrt(x.^2+y.^2)-rE));
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
