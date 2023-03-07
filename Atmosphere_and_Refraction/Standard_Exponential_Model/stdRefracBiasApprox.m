function [deltaROneWay,deltaTheta]=stdRefracBiasApprox(L,thetaEl,radarHeight,Ns,ce,rE,algorithm)
%%STDREFRACTBIASAPPROX Approximate the offsets in range and in elevation
%  angle accounting for how a standard exponential atmospheric model warps
%  a measurement. This only accounts for one-way range and measured
%  elevation. Both algorithms use a spherical Earth approximation.
%
%INPUTS: L The straight-line distance from the sensor to the target. This
%          is typically in meters.
%  thetaEl The elevation angle above the local tangent plane (with respect
%          to a spherical Earth) at the sensor from the receiver to the
%          target in radians.
% radarHeight The height of the radar above the surface of the reference
%          sphere.
%       Ns The atmospheric refractivity reduced to the reference sphere.
%          Note that the refractivity is (n-1)*1e6, where n is the index
%          of refraction. The function reduceStdRefrac2Spher can be used
%          to reduce a refractivity to the surface of a reference
%          ellipsoid. This function does not allow different
%          refractivities to be used as the transmitter and receiver. If
%          this parameter is omitted or an empty matrix is passed, a
%          default value of 313 is used.
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
%       rE The radius of the Earth to use for the spherical Earth
%          approximation. This is typically in meters. If this is omitted
%          or an empty matrix is passed, then the default of
%          Constants.WGS84MeanRadius is used.
% algorithm Choose the algorithm to use. Possible values are:
%          0 Use the fast approximation of [1]. This is not very accurate
%            (even assuming that the exponential atmospheric model is the
%            truth), but it it fast. This method does not work for
%            elevation angles over 49 degrees.
%          1 (The default if omitted or an empty matrix is passed) Use the
%            method described in [2], which requires solving a boundary
%            value problem and performing numerical integration.
%
%OUTPUTS: deltaROneWay The one-way range bias from the sensor to the
%                      target.
%           deltaTheta The elevation angle bias (with elevation measured
%                      above the local tangent plane in the spherical Earth
%                      approximation) of the measurement.
%
%The use of an exponential atmospheric model is often not the most advanced
%refraction model, but it is simple.
%
%Above 10 degrees elevation angle an approximation for the difference of
%the error functions in [1] is used if algorthm 0 is chosen. However, the
%approximation given in the paper does not work. Instead, the error
%function is treated as 1-the complementary error function and the first
%term in an asymptotic expansion (a Poincaré expansion) of the
%complementary error function is used. Thus
%erfc(x)=(approx)1/(x*sqrt(pi))*exp(-x^2)
%
%EXAMPLE 1:
%Here, we consider the difference in the bias when using each algorithm.
%This also shows how to use the osculating sphere approximation including
%the spherical center offset:
% latLonRx=deg2rad([20.269202;-155.852051]);
% AltRx=100;
% latLonTar=deg2rad([20.835390;-155.313721]);
% AltTar=13e3;%13km target altitude.
% %Convert locations to Cartesian.
% zRx=ellips2Cart([latLonRx;AltRx]);
% zTar=ellips2Cart([latLonTar;AltTar]);
% [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx));
% ce=1.593321419785181e-04;
% Ns=350;%Assumed refractivity at the sea surface.
% 
% %Convert to Cartesian coordinates with respect to to osculating sphere.
% zRxLocal=zRx-spherCent;
% zTarLocal=zTar-spherCent;
% vec2Tar=zTarLocal-zRxLocal;
% 
% uENU=getENUAxes(Cart2Ellipse(zRxLocal,[],rE,0));
% %Get the local offset of the target in local coordinates, where "up" is one 
% %coordinate axis and the other is in the local tangent plane. 
% yTar=uENU(:,3)'*vec2Tar;%How far "up" it goes.
% %How far in the tangent plane it goes.
% xTar=sqrt((uENU(:,1)'*vec2Tar)^2+(uENU(:,2)'*vec2Tar)^2);
% thetaEl=atan(yTar/xTar);
% L=norm(vec2Tar);
% radarHeight=norm(zRxLocal)-rE;
% [deltaROneWay0,deltaTheta0]=stdRefracBiasApprox(L,thetaEl,radarHeight,Ns,ce,rE,0)
% [deltaROneWay1,deltaTheta1]=stdRefracBiasApprox(L,thetaEl,radarHeight,Ns,ce,rE,1)
%One can see that the two algorithms lead to results that are generally in
%the vicinity of each other.
%
%EXAMPLE 2:
%This just demonstrates that the results of this function with algorithm 1
%are the same as what is implied by the output of Cart2RuvStdRefrac (within
%finite precision limits.
% latLonRx=deg2rad([20.269202;-155.852051]);
% AltRx=0;
% latLonTar=deg2rad([20.835390;-155.313721]);
% AltTar=8e3;%8km target altitude.
% %Convert locations to Cartesian.
% zRx=ellips2Cart([latLonRx;AltRx]);
% zTar=ellips2Cart([latLonTar;AltTar]);
% [rE,spherCent]=osculatingSpher4LatLon(Cart2Ellipse(zRx));
% ce=1.593321419785181e-04;
% Ns=350;%Assumed refractivity at the sea surface.
% 
% %The receiver faces 45 degrees East of North and 15 degrees up from the
% %local ellipsoidal level.
% M=findRFTransParam([latLonRx;AltRx],deg2rad(45),deg2rad(15));
% 
% %Convert to Cartesian coordinates with respect to to osculating sphere.
% zRxLocal=zRx-spherCent;
% zTarLocal=zTar-spherCent;
% vec2Tar=zTarLocal-zRxLocal;
% 
% uENU=getENUAxes(Cart2Ellipse(zRxLocal,[],rE,0));
% %Get the local offset of the target in local coordinates, where "up" is one 
% %coordinate axis and the other is in the local tangent plane. 
% yTar=uENU(:,3)'*vec2Tar;%How far "up" it goes.
% %How far in the tangent plane it goes.
% xTar=sqrt((uENU(:,1)'*vec2Tar)^2+(uENU(:,2)'*vec2Tar)^2);
% thetaEl=atan(yTar/xTar);
% L=norm(vec2Tar);
% radarHeight=norm(zRxLocal)-rE;
% [deltaROneWay,deltaTheta]=stdRefracBiasApprox(L,thetaEl,radarHeight,Ns,ce,rE,1);
% 
% useHalfRange=true;
% includeW=true;%Include third dimension of unit vector.
% z=Cart2RuvStdRefrac(zTar,useHalfRange,zRx,zRx,M,Ns,includeW,ce,rE,spherCent);
% zNoRefrac=Cart2Ruv(zTar,useHalfRange,zRx,zRx,M,includeW);
% rDiffTraced=(z(1)-zNoRefrac(1));
% deltaThetaTraced=angBetweenVecs(z(2:end),zNoRefrac(2:end));
% 
% RelDiffR=(deltaROneWay-rDiffTraced)./rDiffTraced
% RelDiffTheta=(deltaTheta-deltaThetaTraced)./deltaThetaTraced
%
%REFERENCES:
%[1] J. C. Kerce, W. D. Blair, and G. C. Brown, "Modeling refraction errors
%    for simulation studies of multisensor target tracking," in Proceedings
%    of the Thirty-Sixth Southeastern Symposium on System Theory, Atlanta,
%    GA, 16 Mar. 2004, pp. 97-101.
%[2] D. F. Crouse, "Basic tracking using 3D monostatic and bistatic
%    measurements in refractive environments," IEEE Aerospace and
%    Electronic Systems Magazine, vol. 29, no. 8, Part II, pp. 54-75, Aug.
%    2014.
%
%November 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(algorithm))
    algorithm=1;
end

if(nargin<6||isempty(rE))
    rE=Constants.WGS84MeanRadius;
end

if(nargin<5||isempty(ce))
    expConst=0.005577;
    multConst=7.32;

    %The change in refractivity at an elevation of 1km based on the
    %refractivity on the surface of the Earth.
    DeltaN=-multConst*exp(expConst*Ns);
    ce=log(Ns/(Ns+DeltaN))/1000;%Units of inverse meters.
end

if(nargin<4||isempty(Ns))
    Ns=313;
end

RRadar=rE+radarHeight;

if(algorithm==0)
    alpha=Ns/1e6;
    beta=ce;
    
    if(thetaEl>deg2rad(49))
        error('This algorithm does not work for angles >49 degrees.')
    end

    FVal=F(beta,L,thetaEl,RRadar,rE);
    deltaROneWay=alpha*FVal;
    deltaTheta=alpha*beta*cos(thetaEl)*FVal;
else
    slope=tan(thetaEl);
    x0Init=0;
    x1Init=L*cos(thetaEl);
    y0Init=RRadar;
    y1Init=RRadar+sin(thetaEl)*L;

    %If the two points are nearly vertical, then the ray tracing algorithm will
    %fail. For nearly vertical points, the bending due to refraction in the
    %model should be negligible, so we can perform an integral in the y
    %direction to solve for the excess range instead of having to solve the
    %more complicated general bistatic problem.
    if(abs(pi/2-thetaEl)<1e-3)
        yMax=RRadar+L;
        deltaROneWay=((exp(ce*(rE-y0Init))-exp(ce*(rE-yMax)))*Ns)/(1e6*ce);
        deltaTheta=0;
        return;
    end

    %The initial guess is just the linear solution. The solver requires a fixed
    %number of steps. 20 is probably sufficient for things near the Earth. that
    %is, up to distances of, say 400km. We can scale the number of steps as 20
    %for every 400 kilometers with a minum of, say 10.
    %Things outside of the atmosphere should use the astronomical refraction
    %routines.
    numSteps=max(20,ceil(20*L/400e3));
    x=linspace(x0Init,x1Init,numSteps);
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
    deltaROneWay=range-L;
    deltaTheta=atan(sol.y(2,1))-thetaEl;
end
end

function val=F(beta,L,thetaEl,RRadar,rE)
%% F This implements the F function from [1].
%
%REFERENCES:
%[1] J. C. Kerce, W. D. Blair, and G. C. Brown, "Modeling refraction errors
%    for simulation studies of multisensor target tracking," in Proceedings
%    of the Thirty-Sixth Southeastern Symposium on System Theory, Atlanta,
%    GA, 16 Mar. 2004, pp. 97-101.

    lUpper=sqrt(beta/(2*RRadar))*cos(thetaEl)*(L+RRadar*sec(thetaEl)*tan(thetaEl));
    lLower=sqrt((beta*RRadar)/2)*tan(thetaEl);
    
    if(thetaEl<=deg2rad(10))
        erfDiff=(erf(lUpper)-erf(lLower));
    else
        %This is corrected from [1]. This comes from an asymptotic
        %expansion of the complementary error function.
        erfDiff=1/(lLower*sqrt(pi))*exp(-lLower^2)-1/(lUpper*sqrt(pi))*exp(-lUpper^2) ;
    end

    val=sqrt(pi*RRadar/(2*beta))*cos(thetaEl)*exp(-beta*(RRadar-rE))*...
        exp(beta*RRadar*tan(thetaEl)^2/2)*erfDiff;
end

function val=pathFun2D(x,sol,Ns,rE,ce)
    %This function is used to integrate the time taken.
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
