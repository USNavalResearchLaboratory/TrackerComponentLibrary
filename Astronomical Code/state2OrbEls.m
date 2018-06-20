function orbEls=state2OrbEls(stateVec,elType,GM,epsVal)
%%STATE2ORBELS Convert a state vector consisting of at least position and
%              velocity into a set of orbital elements. Such elements are
%              for a simple two-body problem where the satellite (object)
%              in motion has negligible mass. The elements can be found
%              with respect to a given epoch time.
%
%INPUTS: stateVec A 6XnumVec matrix of numVec state vectors consisting of
%                 position and velocity in a Cartesian (quasi)-inertial
%                 coordinate system where the gravitating body is at the
%                 origin. Extra state elements are ignored. The units are
%                 assumed to be meters and meters per second.
%          elType A value indicating the type of orbital elements. Possible
%                 values are:
%                 0 (The default if omitted) The elements are Gooding's
%                   universal orbital elements.
%                 2 The elements are direct equinoctial orbital elements.
%                 3 The elements are retrograde equnoctial orbital
%                   elements.
%              GM Optionally, the universal gravitation constant times the
%                 mass of the body about which the object with the given
%                 state vector is orbiting. If this parameter is omitted,
%                 the value in Constants.WGS84GMWithAtmosphere is used. The
%                 units are m^3/sec^2.
%          epsVal If universal orbital elements are chosen, this is a
%                 precision bound used for equality comparisons when
%                 determining whether a trajectory is rectilinear,
%                 parabolic or circular. If omitted, a default value of
%                 eps is used.
%
%OUTPUTS: orbEls An 6XnumVec set of vectors of orbital elements, the format
%                of which depends on the elType parameter. The format of
%                the elements is discussed below.
%
%Gooding's universal orbital elements are consist of (in order)
%alpha=GM/a where a is the semi-major axis in meters
%q=a*(1-e)  where e is the eccentricity (unitless). This is the perifocal
%           distance.
%i          inclination in radians
%Omega      longitude (right ascension) of the ascending node in radians
%omega      argument of periapsis/perigee in radians
%tau        the time at pericenter (seconds)
%
%The equinoctial orbital elements (for both direct and retrograde elements)
%are
%a      semi-major axis in meters
%h      first eccentricity vector component (unitless)
%k      second eccentricity vector component (unitless)
%p      first ascending node vector component (unitless)
%q      second ascending node vector component (unitless)
%lambda mean longitude in radians.
%
%The algorithm using Gooding's universal orbital elements is very robust to
%numerical errors. It is implemented using the algorithm of [1], where the
%elements are also described.
%
%Equinoctial orbital elements are discussed in Section 2 of [2]. and in
%[3]. The diffference between direct and retrograde elements is essentially
%a matter of handedness.
%
%Orbital elements are discussed in general in Chapter 2.2 of [4].
%
%The inverse of this function is orbEls2State, which can also predict the
%state forward in time. Adding a time interval in seconds to the tau
%argument also predicts the trajectory forward in terms of universal
%orbital elements.
%
%REFERENCES:
%[1] R. H. Gooding, "On universal elements, and conversion procedures to
%    and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%    pp. 283-298, 1988.
%[2] D. A. Danielson, C. P. Sagovac, B. Neta, and L. W. Early,
%    "Semianalytic satellite theory," Mathematics Department, Naval
%    Postgraduate School, Monterey, CA, Tech. Rep., 1995. [Online].
%    Available: http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix= html&identifier=ADA531136
%[3] R. A. Broucke and P. J. Cefola, "On the equinoctial orbit elements,"
%    Celestial Mechanics, vol. 5, no. 3, pp. 303-310, 1972.
%[4] O. Montenbruck and E. Gill, Satellite Orbits: Models, Methods
%    Applications. Berlin: Springer, 2000.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    epsVal=eps;
end

if(nargin<3)
   GM=Constants.WGS84GMWithAtmosphere;
end

if(nargin<2)
    elType=0;
end

numVec=size(stateVec,2);
orbEls=zeros(6,numVec);

switch(elType)
    case 0
        for curVec=1:numVec
            orbEls(:,curVec)=state2OrbElsUniv(stateVec(1:6,curVec),GM,epsVal);
        end
    case 1
        for curVec=1:numVec
            orbEls(:,curVec)=state2OrbElsEquinoctial(stateVec(1:6,curVec),1,GM);
        end
    case 2
        for curVec=1:numVec
            orbEls(:,curVec)=state2OrbElsEquinoctial(stateVec(1:6,curVec),-1,GM);
        end
    otherwise
        error('Invalid element type provided.')
end

end

function elEqui=state2OrbElsEquinoctial(state,I,GM)
%%STATE2ORBELSEQUINOCTIAL An algorithm to convert from a target state to a
%                         set of equinoctial orbital elements.
%
%The input I is the retrograde factor ond is +1 for direct equinoctial
%elements and -1 for retrograde equinoctial elements.
%
%The algorithm is taken from Section 2.1.5 of
%D. A. Danielson, C. P. Sagovac, B. Neta, and L. W. Early, "Semianalytic
%satellite theory," Mathematics Department, Naval Postgraduate School,
%Monterey, CA, Tech. Rep., 1995. [Online]. Available:
%http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix= html&identifier=ADA531136
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

%position
r=state(1:3);
%velocity
rDot=state(4:6);

%Equation 1
a=1/(2/norm(r)-norm(rDot)^2/GM);%Semi-major axis.
%Equation 2
w=cross(r,rDot)/norm(cross(r,rDot));

%Equation 3
p=w(1)/(1+I*w(3));%First ascending node vector component.
q=-w(2)/(1+I*w(3));%Second ascending node vector component.

%Equation 4
e=-r/norm(r)+cross(rDot,cross(r,rDot))/GM;

%Equation 2.1.4-1
coeff=(1/(1+p^2+q^2));
f=coeff*[1-p^2+q^2;
         2*p*q;
         -2*I*p];
g=coeff*[2*I*p*q;
         (1+p^2-q^2)*I;
         2*q];
     
%Equation 5
h=dot(e,g);%First eccentricity vector component.
k=dot(e,f);%Second eccentricity vector component.

%Equation 6
X=dot(r,f);
Y=dot(r,g);

%Equation 2.1.4-4
b=1/(1+sqrt(1-h^2-k^2));
%Equation 7
denom=a*sqrt(1-h^2-k^2);
sinF=h+((1-h^2*b)*Y-h*k*b*X)/denom;
cosF=k+((1-k^2*b)*X-h*k*b*Y)/denom;
F=atan2(sinF,cosF);

lambda=F+h*cosF-k*sinF;%Mean longitude

elEqui=[a;h;k;p;q;lambda];

end


function elUniv=state2OrbElsUniv(state,GM,epsVal)
%%STATE2ORBELSUNIV An implementation of Gooding's conversion from a target
%                  state to universal orbital elements. This implements the
%                  PV3ELS function.
%
%The implementation is taken from Appendix D of
%R. H. Gooding, "On universal elements, and conversion procedures to
%and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%pp. 283-298, 1988.
%with minor changes so that equality comparisons are replaced with
%comparisons within a region epsVal.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

rVec=state(1:3);
vVec=state(4:6);
x=rVec(1);
y=rVec(2);
z=rVec(3);

xyMag=norm(rVec(1:2));
r=norm(rVec);

%The projection of the velocity onto the radial component.
VR=dot(rVec,vVec)/r;%Radial velocity

%The angular momentum vector.
hVec=cross(rVec,vVec);

if(norm(hVec)^2<epsVal)%Rectilinear orbit
    %This works for an axial orbit as well as for a general rectilinear
    %orbit.
    i=pi/2;%inclination
    %For axial or general rectilinear orbits, atan2 will provide the
    %correct result. When axial and both are zero, it will return 0.
    Omega=atan2(y,x);%Longitude of the ascending node.
    %Angle from assumed reference direction (argument of latitude)
    u=atan2(z,xyMag);
    VT=0;%Transverse velocity
else%Non degenerate orbit
    b=cross(hVec,rVec);
    
    %r2h is r^2*h as explained before the beginning of Section 5.
    r2h=cross(rVec,b);
    W=dot(r2h(1:2),r2h(1:2));
    
    i=atan2(sqrt(W),r2h(3));%inclination
    
    if(W<epsVal)%If the orbit is in the reference plane
        Omega=0;%Longitude of the ascending node.
        u=atan2(y*sign(r2h(3)),x);%Angle from assumed reference direction
    else%General orbit
        Omega=atan2(r2h(1),-r2h(2));%Longitude of the ascending node.
        %Angle from assumed reference direction (argument of latitude)
        u=atan2(norm(r2h)*z,r^2*b(3));
    end
    VT=norm(r2h)/r^3;%Transverse velocity
end

[alpha,q,omega,tau]=PV2ELS(VR,VT,r,u,GM,epsVal);

elUniv=[alpha;
            q;
            i;
        Omega;
        omega;
          tau];
end

function [alpha,q,omega,tau]=PV2ELS(VR,VT,r,u,GM,epsVal)
%%PV2ELS  An implementation of the PV2ELS subroutine of Gooding's algorithm
%        for universal orbital element conversion. This is a
%        two-dimensional conversion suroutine.
%
%INPUTS VR, VT Radial and transerver velocity. Note VT>=0.
%       r      Radial distance.
%       u      Angle from assumed reference direction.
%       GM     Universal gravitational constant times the mass of the
%              gravitating body.
%       epsVal A small number for determining equality within finite
%              precision bounds.
%
%The algorithm is taken from Appendix C of
%R. H. Gooding, "On universal elements, and conversion procedures to
%and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%pp. 283-298, 1988.
%with minor changes so that equality comparisons are not performed.  
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

sw=0.25;

VMag2=VR^2+VT^2;
alpha=2*GM/r-VMag2;

d=r*VR;
h=r*VT;
p=h^2;

esq1=p*alpha;
es=d*sqrt(abs(alpha));
ec=r*VMag2-GM;
if(alpha>0)%One formula superior for the ellipse
    e=sqrt(ec^2+es^2);
else%Different formula superior for the hyperbola
    e=sqrt(GM^2-esq1);
end

q=p/(GM+e);
if(abs(alpha)<=epsVal)%Parabola
    tau=d*(2*q+r)/(3*GM);
    v=2*atan2(VR,VT);%The true anomaly
elseif(abs(e)<=epsVal)%Circle
    tau=0;
    v=0;%The true anomaly
else%Ellipse or hyperbola
    e1=alpha*q;
    
    if(alpha>0)%Ellipse
        eh=atan2(es,ec);
        if(GM*eh^2/6+e1>GM*sw)%General case
            em=GM*eh-es;
            ecesq=GM*ec-e^2;
        else%For e1 and eh both near zero
            em=GM*GoodingSinDiffFun(e1/GM,eh);
            ecesq=(esq1*ec^2-e^2*es^2)/(e^2+GM*ec);
        end
    else%Hyperbola
        eh=asinh(es/e);
        
        if(GM*eh^2/6-e1>GM*sw)%General case
            em=es-GM*eh;
            ecesq=e^2-GM*ec;
        else%For e1 and eh both near zero
            em=e*GoodingHyperSinDiffFun(-e1/e,es/e);
            ecesq=-(esq1*ec^2+e^2*es^2)/(e^2+GM*ec);
        end
    end
    %Still ellipse or hyperbola
    en=abs(alpha)^(3/2);
    tau=em/en;
    v=atan2(es*h*sqrt(abs(alpha)),ecesq);%The true anomaly
end
%All orbits
omega=u-v;%The argument of periapsis

%Adjust revolutions, if necessary, so that omega remains in the range
%-pi to pi. The second condition in the if-statement takes care of possible
%parabolic case; this is just for the elliptical case.
if(alpha>0&&abs(alpha)>epsVal)
    adj=2*pi*fix(abs(omega/(2*pi))+(1/2))*sign(omega);
    omega=omega-adj;
    tau=tau+adj/en;
end

end


function x=GoodingSinDiffFun(e,EE)
%%GOODINGSINDIFFFUN Evaluate the function EE-(1-e)*sin(EE) using Gooding's
%                   EMKEP procedure for when e and EE are close
%                   to (1,0). This is supposed to be more accurate than
%                   just directly evaluating the functions, unless EE is
%                   large as it is then supposed to worsen rounding errors.
%
%
%The algorithm is the EMKPL algorithm taken from Appendix C of
%A. W. Odell and R. H. Gooding, "Procedure for solving Kepler's
%equation," Celestial Mechanics, vol. 38, no. 4, pp. 307-334, Apr. 1986.
%modified to solve EE-(1-e)*sin(EE) instead of EE-e*sin(EE).
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

x=e*sin(EE);
EE2=-EE^2;
term=EE;
d=0;
while(1)
    d=d+2;
    term=term*EE2/(d*(d+1));
    x0=x;
    x=x-term;
    if(x==x0)
        break;
    end
end

end


function x=GoodingHyperSinDiffFun(g1,s)
%%GOODINGHYPERSINDIFFFUN Evaluate the function s-(1-g1)*asinh(s) when
%                        (g1,s) is close to (0,0) using Gooding's method.
%                        This is supposed to have a higher precision than
%                        just explicitly evaluating the function.
%
%The algorithm is the SHMKEP function taken from Appendix B of 
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

g=1-g1;
t=s/(1+sqrt(1+s^2));
x=s*(g1+g*t^2);
term=2*g*t;
twoI1=1;
%Iterate until convergence or until a maximum number of iterations is
%reached.
maxIter=64;
for curIter=1:maxIter
    twoI1=twoI1+2;
    term=term*t^2;
    x0=x;
    x=x-term/twoI1;
    if(x==x0)
        break;
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

