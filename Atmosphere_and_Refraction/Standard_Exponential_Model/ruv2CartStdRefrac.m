function zCart=ruv2CartStdRefrac(zRUVBiased,useHalfRange,zTx,zRx,M,Ns,ce,rE,spherCent,xMax)
%%RUV2CARTSTDREFRAC Convert points in refraction-corrupted bistatic r-u-v
%         (or r-u-v-w) coordinates into Cartesian coordinates using a
%         standard exponential atmospheric model.  If r-u-v coordinates are
%         used, the target is assumed to be in front of the receiver (local
%         z coordinate is positive). r-u-v coordinates consist of a
%         bistatic range and direction cosines at the receiver. The
%         "direction cosines" u and v are just the x and y coordinates of a
%         unit vector from the receiver to the target in the coordinate
%         system at the receiver. This basically assumes that the boresight
%         direction of the receiver is the z axis. Assuming the target is
%         in front of the receiver, the third unit vector coordinate is not
%         needed. However, r-u-v-w coordinates include the third component.
%         For monostatic coordinate conversion where the range is a one-way
%         range, set useHalfRange=true and zTx and zRx to the same value.
%         This function is not suitable for computing refraction between
%         satellites, grazing the Earth's atmosphere. The algorithm might
%         have an error if the ray path goes too far underground.
%
%INPUTS: z A 3XN matrix of refraction-corrupted measurements with elements
%          [r;u;v], where r is the bistatic range from the transmitter to
%          the target to the receiver, and u and v are direction cosines.
%          Each u,v pair should have a magnitude less than or equal to
%          one. If the magnitude is greater than one, then the pair is
%          normalized before conversion to avoid imaginary results.
%          Alternatively, one can pass a 4XN matrix of [r;u;v;w] vectors
%          where [u;v;w] form a full unit vector in the receiver's local
%          3D Cartesian coordinates.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is false.
%      zTx The 3X1 [x;y;z] location vector of the transmitter in global
%          Cartesian coordinates.  If this parameter is omitted or an
%          empty matrix is passed, then the receiver is placed at the origin.
%      zRx The 3X1 [x;y;z] location vector of the receiver in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to the local alignment fo the receiver. The z
%          vector of the local coordinate system of the receiver is the
%          pointing direction of the receiver. If this matrix is omitted,
%          then the identity matrix is used.
%       Ns The atmospheric refractivity reduced to the reference sphere.
%          Note that the refractivity is (n-1)*1e6, where n is the index of
%          refraction. The function reduceStdRefrac2Spher can be used to
%          reduce a refractivity to the surface of a reference ellipsoid.
%          This function does not allow different refractivities to be
%          used as the transmitter and receiver. If this parameter is
%          omitted or an empty matrix is passed, a default value of 313 is
%          used.
%       ce The optional decay constant of the exponential model. The
%          refractivity N at height h is N=Ns*exp(-ce*(h-h0)) where h0 is
%          the reference height (in this function, the height of the
%          reference ellipsoid surface is used). ce is related to the
%          change in refractivity at an elevation of 1km based on the
%          refractivity at sea level as
%          ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
%          where the change in refractivity for a change in elevation of
%          1km is DeltaN=-multConst*exp(expConst*Ns); In [1], standard
%          values for the two constants are expConst=0.005577; and
%          multConst=7.32; If ce is omitted or an empty matrix is passed,
%          the value based on the standard model is used.
% rE,spherCent The radius of the Earth to use for the spherical Earth
%           approximation used in the model and also the offset between the
%           global model and the local spherical model. It is assumed that
%           zC,zTx,and zRx are all given in the global model and will need
%           to be transformed to the local model to the used. If rE is
%           omitted or an empty matrix is passed, then the default of
%           [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx)) is
%           used. The defaults here mean that a WGS-84 reference ellipsoid
%           is approximated by the local osculating sphere.
%     xMax This function traces the ray to a maximum displacement in the
%          local tangent plane prior (or vertically for nearly vertical
%          directions) to attempting to find a particular range. This is an
%          optional parameter specifying the maximum distance to trace the
%          ray. The default if this parameter is omitted or an empty matrix
%          is passed is 1000e3 (1000km).
%
%This function implements the algorithm of [1]. If the target is collocated
%with the transmitter or the receiver, then the algorithm will fail. The
%basic exponential refraction model is in [2].
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
%The algorithm of [1] performs ray tracing by solving an initial value
%problem and performing a line search.
%
%For paths that are nearly vertical, it is approximated that there is no
%bending in angle and an explicit solution to the integral over the index
%of refraction in the vertical direction is used to obtain the range offset
%in the vertical direction.
%
%EXAMPLE:
%In an example near Hawaii, we will convert a position into bistatic r-u-v
%coordinates using Cart2RuvStdRefrac and then we will convert it back using
%this function.
% latLonRx=deg2rad([20.269202;-155.852051]);%Degrees converted to radians.
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
% zRUVBiased=Cart2RuvStdRefrac(zTar,useHalfRange,zTx,zRx,M,Ns,includeW);
% %Now that we have the bistatic r-u-v value, we will convert it back to
% %Cartesian 
% zTarCart=ruv2CartStdRefrac(zRUVBiased,useHalfRange,zTx,zRx,M,Ns);
% norm(zTarCart-zTar)
% %One will see that zTarCart is less than 1mm away from zTar
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%[1] B. R. Bean and G. D. Thayer, CRPL Exponential Reference Atmosphere.
%    Washington, D.C.: U. S. Department of Commerce, National Bureau of
%    Standards, Oct. 1959. [Online]. Available:
%    http://digicoll.manoa.hawaii.edu/techreports/PDF/NBS4.pdf
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10||isempty(xMax))
    xMax=1000e3;%1000 kilometer assumed maximum x displacement.
end

if(nargin<8||isempty(rE))
    %Use the radius of the Earth that is the radius of the osculating
    %sphere at the location of the observed. This will be the radius used
    %in the local spherical Earth approximation for computing atmospheric
    %refraction. This uses the WGS-84 reference ellipsoid.
    [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx));
end

%Transform the sensor locations and use spherCent=0 in the following
%computation. In the end, we will have to transform the target location
%back to the global sphere.
zRx=zRx-spherCent;
zTx=zTx-spherCent;

numMeas=size(zRUVBiased,2);

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

if(nargin<7||isempty(ce))
    expConst=0.005577;
    multConst=7.32;

    %The change in refractivity at an elevation of 1km based on the
    %refractivity on the surface of the Earth.
    DeltaN=-multConst*exp(expConst*Ns);
    ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
end

if(useHalfRange)
    zRUVBiased(1,:)=2*zRUVBiased(1,:);
    useHalfRange=false;
end

%If it is monostatic, then the extra boundary value problem for the
%target-transmitter path need not be solved.
if(all(zTx==zRx))
    isMonostatic=true;
else
    isMonostatic=false;
end

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
uENU=getENUAxes(Cart2Ellipse(zRx,[],rE,0));
uVertGlobal=uENU(:,3);
uVertLocal=[0;1;0];

zCart=zeros(3,numMeas);

for curMeas=1:numMeas
    %The apparent position in global coordinates.
    tCartBiased=ruv2Cart(zRUVBiased(:,curMeas),useHalfRange,zTx,zRx,M);

    vec2TarGlobal=tCartBiased-zRx;

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
    y0Init=norm(zRx);

    %The location of the target in the local 2D coordinate system.
    x1Init=vec2TarLocal(1);
    y1Init=vec2TarLocal(2)+y0Init;

    %The derivative of y at x=0 is just the slope of the line in the apparent
    %direction of the target.
    y0Dot=(y1Init-y0Init)/(x1Init-x0Init);

    if(norm(uHorizGlobalOrig)<1e-3)
        oldOpts=optimset();
        newOpts=optimset(oldOpts,'TolX',1e-8);

        %If the vector passed only has uv and not a complete unit vector, add
        %the third element.
        if(size(zRUVBiased,1)==3)
            %The real command just tries to add some robustness if the 
            %magnitudes of u and v are too large.
            zMeas=[zRUVBiased(:,curMeas);real(sqrt(1-sum(zRUVBiased(2:3,curMeas).^2)))];
        else
            zMeas=zRUVBiased(:,curMeas);
        end

        yTrue=fminbnd(@(y)rangeCostVertical(y,y0Init,isMonostatic,zTx,zMeas,Ns,rE,ce),y0Init,y0Init+xMax,newOpts);

        %Convert back into ECEF.
        zCart(:,curMeas)=ECEF2LocalRot'*[0;yTrue-y0Init;0]+zRx+spherCent;
        return;
    end

    %We now have the initial conditions and can solve the initial value
    %problem.
    xSpan=[x0Init;xMax];%The range of values over which it will be solved.
    initialCond=[y0Init;y0Dot];

    %Up the accuracy
    oldOpts=odeset();
    newOpts=odeset(oldOpts,'RelTol',1e-12,'AbsTol',1e-12,'Jacobian',@odefunJacob);
    sol=ode45(@(x,y)expDiffEq(x,y,Ns,rE,ce),xSpan,initialCond,newOpts);
    %The solution sol is the monostatic traced path. We must find the value of
    %x such that the desired range is acquired after integrating over the path.

    oldOpts=optimset();
    newOpts=optimset(oldOpts,'TolX',1e-8);
    xTrue=fminbnd(@(x)rangeCost(x,isMonostatic,sol,zTx,zRx,ECEF2LocalRot,y0Init,Ns,rE,ce,zRUVBiased(1,curMeas)),0,xMax,newOpts);

    yTrue=deval(sol,xTrue);
    yTrue=yTrue(1);%yTrue(1) is the position.

    %Convert back into ECEF.
    zCart(:,curMeas)=ECEF2LocalRot'*[xTrue;yTrue-y0Init;0]+zRx+spherCent;
end
end

function val=rangeCost(x,isMonostatic,sol,zTx,zRx,ECEF2LocalRot,y0Init,Ns,rE,ce,biasedRange)
%The standard bistatic cost function for determining the range along the
%line of sight when the line of sight is not (almost) directly above the
%receiver.

    %Evaluate the range integral from the receiver to the target.
    rRx=integral(@(x)pathFun2D(x,sol,Ns,rE,ce),0,x,'RelTol',1e-10,'AbsTol',1e-10);

    if(isMonostatic==true)
        rTx=rRx;%No need to do the extra ray tracing in the monostatic case.
    else%It is bistatic
        %Get the location of the target at this point and solve the
        %boundary value problem to get the range from the target location
        %hypothesis to the transmitter.
        y=deval(sol,x);
        
        vecLocal=[x;y(1)-y0Init;0];
        %Convert back into ECEF.
        tarLocGlobalCur=ECEF2LocalRot'*vecLocal+zRx;
        
        %Now, solve for the apparent one-way range between the target and 
        %the transmitter.
        z=Cart2RuvStdRefrac(tarLocGlobalCur,true,zTx,zTx,[],Ns,[],ce,rE,0);
        rTx=z(1);
    end
    %Return the squared difference between the integrated value and the
    %actual value.
    val=biasedRange-rTx-rRx;
    val=val^2;
end

function val=rangeCostVertical(yMax,y0Init,isMonostatic,zTx,zMeas,Ns,rE,ce)
    rRx=((exp(ce*(rE-y0Init))-exp(ce*(rE-yMax)))*Ns)/(1e6*ce)+yMax-y0Init;    

    if(isMonostatic==true)
        rTx=rRx;%No need to do the extra ray tracing in the monostatic case.
    else%It is bistatic
        %Get the location of the target at this point and solve the
        %boundary value problem to get the range from the target location
        %hypothesis to the transmitter.
        tarLocGlobalCur=zMeas(2:end)*rRx;
        
        %Now, solve for the apparent one-way range between the target and 
        %the transmitter.
        z=Cart2RuvStdRefrac(tarLocGlobalCur,true,zTx,zTx,[],Ns,[],ce,rE,0);
        rTx=z(1);
    end
    %Return the squared difference between the integrated value and the
    %actual value.
    val=zMeas(1)-rTx-rRx;
    val=val^2;
end


function val=pathFun2D(x,sol,Ns,rE,ce)
    %This function is used to integrate the time taken
    y=deval(x,sol);
    val=(1+NRefracExp(x,y(1,:),Ns,rE,ce)).*sqrt(1+y(2,:).^2);
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

function [nRefrac,ce]=NRefracExp(x,y,Ns,rE,ce)
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
