function [aDeriv,AJacob]=aJ1Gravity(xState,includeCoriolis,GM,omega)
%%AJ1GRAVITY This is the J1 dynamic model for simple ballistic targets near
%            the Earth, neglecting atmospheric drag. It is assumed that the
%            global origin is at the center of the Earth. If one is
%            tracking in an Earth-centered Earth-fixed (ECEF) coordinate
%            system, then the Coriolis effect should be included; in an
%            Earth-centered inertial (ECI) coordinate system, it should be
%            omitted. It is assumed that the rotation axis of the Earth
%            coincides with the z axis.
%
%INPUTS: xState A 6X1 target state consisting of 3D Cartesian
%               [position;velocity]. The position must be in an
%               Earth-centered coordinate system, such as the international
%               terrestrial reference system (ITRS). If an Earth-centered
%               inertial system, such as the geocentric celestial reference
%               system (GCRS), this function does not take into
%               account the rotation between the ITRS and the GCRS.
% includeCoriolis An optional parameter indicating whether the Coriolis
%               and centrifugal forces (due to the Earth's rotation) should
%               be taken into account. This is needed when tracking in an
%               ECEF coordinate system. The default if this parameter is
%               omitted is true.
%            GM The universal gravitational constant times the mass of the
%               Earth. Units are typically m^3/s^2. If this parameter is
%               omitted or an empty matrix is passed, then the default of
%               Constants.EGM2008GM is used.
%         omega If the Coriolis effect and the centrifugal acceleration due
%               to the Earth's rotation should be taken into account, this
%               is the rotation rate of the Earth, usually in radians per
%               second. If this parameter is omitted or an empty matrix is
%               passed, then omega=Constants.EGM2008EarthRotationRate is
%               used.
%
%OUTPUTS: aDeriv The 6X1 derivative of the state xState with respect to
%                time. The elements are [velocity;acceleration].
%         AJacob The Jacobian of aDeriv, which can be useful in extended
%                Kalman filters. This is the derivative of each component
%                of aDeriv (selected by row) with respect to the elements
%                of the state [x,y,z,xDot,yDot,zDot] selected by column.
%
%The J1 model is simply Newton's law for a point mass. It is the J2 model
%without the oblateness term. The J2 model is derived in Appendix F.4 of
%[1]. The Coriolis terms are derived in Appendix A of [1]. 
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(includeCoriolis))
   includeCoriolis=true; 
end

if(nargin<3||isempty(GM))
   %Universal gravitational constant times the mass of the Earth
    GM=Constants.EGM2008GM; 
end

if(nargin<4||isempty(omega))
   omega=Constants.EGM2008EarthRotationRate;
end

r=xState(1:3);%Position
x=r(1);
y=r(2);
z=r(3);

v=xState(4:6);%Velocity
rMag2=dot(r,r);
rMag=sqrt(rMag2);
rMag3=rMag*rMag2;

%Newton's law for a point or sphere.
aNewton=-(GM/rMag3)*r;

if(includeCoriolis)
    %The rotation vector for the Earth is
    %Omega=[0;0;1]*omega;
    %and the Coriolis and centrifugal forces to add to accel are
    %aCoriolis=-2*bsxfun(@cross,Omega,v);
    %aCentrifugal=-bsxfun(@cross,Omega,bsxfun(@cross,Omega,r));

    %However, it is faster to add the components directly:
    aCoriol=[2*omega*v(2,:)+omega^2*x;-2*omega*v(1,:)+omega^2*y;0];
else
    aCoriol=[0;0;0];
end

a=aNewton+aCoriol;

%The time-derivative of the state vector under the J1 dynamic model.
aDeriv=[v;a];

if(nargout>1)
    rMag5=rMag3*rMag2;
    x2=x*x;
    y2=y*y;
    z2=z*z;

    axx=(GM*(2*x2-y2-z2))/rMag5;
    ayy=-(GM*(x2-2*y2+z2))/rMag5;
    azz=-(GM*(x2+y2-2*z2))/rMag5;
    axy=(3*GM*x*y)/rMag5;
    axz=(3*GM*x*z)/rMag5;
    ayz=(3*GM*y*z)/rMag5;
    
    %If Coriolis and centrifugal acceleration terms should be added due to
    %the rotation of the Earth.
    if(includeCoriolis)
        daXCorioldx=omega^2;
        daXCorioldvy=2*omega;
        daYCorioldy=omega^2;
        daYCorioldvx=-2*omega;
    else
        daXCorioldx=0;
        daXCorioldvy=0;
        daYCorioldy=0;
        daYCorioldvx=0;
    end

    AJacob=[0,                  0,              0,      1,              0,              0;
            0,                  0,              0,      0,              1,              0;
            0,                  0,              0,      0,              0,              1;
            axx+daXCorioldx,    axy,            axz,    0,              daXCorioldvy,   0;
            axy,                ayy+daYCorioldy,ayz,    daYCorioldvx,   0,              0;
            axz,                ayz,            azz,    0,              0,              0];
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
