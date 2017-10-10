function stateVec=orbEls2State(orbEls,deltaT,elType,GM,epsVal)
%%ORBELS2STATE Obtain a velocity and position from a set of orbital
%              elements at a desired time in seconds since the epoch time.
%
%INPUTS: orbEls A 6XnumVec vector of orbital numVec elements to convert.
%               The elements are either Gooding's universal orbital
%               elements or equinoctial orbital elements, as specified by
%               the parameter elType. The components of the elements
%               depends on the type and are explained below.
%        deltaT The time in seconds forward that the trajectory should be
%               propagated from the epoch time of the orbital elements. If
%               this parameter is omitted, a value of 0 is used. For
%               example, if one wishes to propagate the state forward an
%               hour from the time associated with the elements, one can
%               use deltaT=60*60.
%        elType A value indicating the type of orbital elements. Possible
%               values are:
%               0 (The default if omitted) The elements are Gooding's
%                 universal orbital elements.
%               1 The elements are direct equinoctial orbital elements.
%               2 The elements are retrograde equnoctial orbital
%                 elements.
%            GM An optional value of the universal gravitational constant
%               times the mass of the Earth. If omitted, the value
%               Constants.WGS84GMWithAtmosphere is used. The units are
%               m^3/sec^2.
%        epsVal If universal orbital elements are chosen, this is a
%               precision bound used on the alpha term (GM divided by the
%               semi-major axis) for determining whether the trajectory is
%               parabolic. If omitted, a default value of eps is used.
%
%OUTPUTS: stateVec A 6XnumVec matrix of numVec state vectors consisting of
%                 position and velocity in a Cartesian inertial coordinate
%                 system where the gravitating body is at the origin. Units
%                 are meters, and then meters per second.
%
%Gooding's universal orbital elements are consist of (in order)
%alpha=GM/a where a is the semi-major axis in meters
%q=a*(1-e)  where e is the eccentricity (unitless). This is the perifocal
%           distance.
%i          inclination in radians
%Omega      longitude of the ascending node in radians
%omega      argument of periapsis in radians
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
%Equinoctial orbital elements are discussed in Section 2 of [2] and in [3].
%The diffference between direct and retrograde elements is essentially a
%matter of handedness.
%
%Orbital elements are discussed in general in Chapter 2.2 of [4].
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

if(nargin<5)
    epsVal=eps;
end

if(nargin<4)
   GM=Constants.WGS84GMWithAtmosphere; 
end

if(nargin<3)
    elType=0;
end

if(nargin<2)
    deltaT=0;
end

numVecs=size(orbEls,2);
stateVec=zeros(6,numVecs);

switch(elType)
    case 0
        for curVec=1:numVecs
            stateVec(:,curVec)=orbEls2StateUniv(orbEls(:,curVec),deltaT,GM,epsVal);
        end
    case 1
        for curVec=1:numVecs
            stateVec(:,curVec)=orbEls2StateEquinoctial(orbEls(:,curVec),deltaT,1,GM);
        end
    case 2
        for curVec=1:numVecs
            stateVec(:,curVec)=orbEls2StateEquinoctial(orbEls(:,curVec),deltaT,-1,GM);
        end
    otherwise
        error('Invalid element type provided.')
end
end


function state=orbEls2StateUniv(theEls,deltaT,GM,epsVal)
%%ORBELS2STATEUNIV  An implementation of Gooding's conversion from
%                   universal orbital elements to state vectors. This
%                   implements the ELS3PV function.
%
%The implementation is taken from Appendix B of
%R. H. Gooding, "On universal elements, and conversion procedures to
%and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%pp. 283-298, 1988.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

%If we are given classical elements...
alpha=theEls(1);%alpha=GM/a where a is the semi-major axis in meters
q=theEls(2);%q=a*(1-e), where e=eccentricity. q is the perifocal distance
i=theEls(3);%Inclination
Omega=theEls(4);%Longitude of the ascending node
omega=theEls(5);%Argument of periapsis
tau=theEls(6)+deltaT;%The time at pericenter.

[r,u,VR,VT]=ELS2PV(alpha,q,omega,tau,GM,epsVal);

c=cos(u);
s=sin(u);
x1=r*c;
y1=r*s;
x2=VR*c-VT*s;
y2=VR*s+VT*c;

c=cos(i);
s=sin(i);
z=y1*s;
y1=y1*c;
zDot=y2*s;
y2=y2*c;
c=cos(Omega);
s=sin(Omega);
x=x1*c-y1*s;
y=x1*s+y1*c;
xDot=x2*c-y2*s;
yDot=x2*s+y2*c;

state=[x;y;z;xDot;yDot;zDot];
end

function state=orbEls2StateEquinoctial(theEls,deltaT,I,GM)
%%ORBELS2STATEEQUINOCTIAL An algorithm to convert from equinoctial orbital
%                         elements to a target state vector.
%
%The input I is the retrograde factor ond is +1 for direct equinoctial
%elements and -1 for retrograde equinoctial elements.
%
%The algorithm is taken from Section 2.1.4 of
%D. A. Danielson, C. P. Sagovac, B. Neta, and L. W. Early, "Semianalytic
%satellite theory," Mathematics Department, Naval Postgraduate School,
%Monterey, CA, Tech. Rep., 1995. [Online]. Available:
%http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix= html&identifier=ADA531136
%
%The equation for propagating the elements forward in time is taken from
%R. A. Broucke and P. J. Cefola, "On the equinoctial orbit elements,"
%Celestial Mechanics, vol. 5, no. 3, pp. 303-310, 1972.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

a=theEls(1);%Semi-major axis
h=theEls(2);%First eccentricity vector component
k=theEls(3);%Second eccentricity vector component
p=theEls(4);%First ascending node vector component
q=theEls(5);%Second ascending node vector component
lambda=theEls(6);%Mean longitude (at epoch)

n=sqrt(GM/a^3);%The mean motion
%The propagation of the mean longitude is from Equation 15 in Brouke and
%Cefola.
lambda=lambda+n*deltaT;

%Section 2.1.4 conversion from equinoctial elements to position and
%velocity.
coeff=(1/(1+p^2+q^2));
f=coeff*[1-p^2+q^2;
         2*p*q;
         -2*I*p];
g=coeff*[2*I*p*q;
         (1+p^2-q^2)*I;
         2*q];

%Solve the Equinoctial form of Kepler's Equation. to get F.
F=solveEquinoctialKeplersEq(lambda,h,k);

n=sqrt(GM/a^3);
b=1/(1+sqrt(1-h^2-k^2));
denom=1-h*sin(F)-k*cos(F);
sinL=((1-k^2*b)*sin(F)+h*k*b*cos(F)-h)/denom;
cosL=((1-h^2*b)*cos(F)+h*k*b*sin(F)-k)/denom;

r=a*(1-h^2-k^2)/(1+h*sinL+k*cosL);
X=r*cosL;
Y=r*sinL;
XDot=-n*a*(h+sinL)/sqrt(1-h^2-k^2);
YDot=n*a*(k+cosL)/sqrt(1-h^2-k^2);

r=X*f+Y*g;
rDot=XDot*f+YDot*g;

state=[r;rDot];

end


function [r,u,VR,VT]=ELS2PV(alpha,q,omega,tau,GM,epsVal)
%%ELS2PV An implementation of the ELS2PV subroutine of Gooding's algorithm
%        for universal orbital element conversion. This is a
%        two-dimensional conversion suroutine.
%
%The algorithm is taken from Appendix A of
%R. H. Gooding, "On universal elements, and conversion procedures to
%and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%pp. 283-298, 1988.
%with minor changes so that alpha is compared to an epsilon value rather
%than zero for determining the parabolic case.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(abs(alpha)<=epsVal) %Parabolic case, GM cannot be zero.
    d=solveCubicEqGoodin(0.5/GM,q,1.5*GM*tau);
    r=q+(1/2)*d^2/GM;
    h=sqrt(2*GM*q);
    v=2*atan2(d,h);
else%Ellipse or hyperbola
    e1=alpha*q;
    e=GM-e1;
    ep1=GM+e;
    h=sqrt(q*ep1);
    
    em=tau*abs(alpha)^(3/2);
    if(alpha>0)%ellipse, GM cannot be zero
        ee2=(1/2)*solveKeplersEq(em/GM,e1/GM,2);
        s2=sin(ee2);
        c2=cos(ee2);
        r=q+2*e*s2^2/alpha;
        d=2*e*s2*c2/sqrt(abs(alpha));
        v=2*atan2(ep1*s2,h*sqrt(abs(alpha))*c2);
        emv=em/GM-v;
        v=v+4*pi*fix(abs(emv/(4*pi))+(1/2))*sign(emv);
    else%hyperbola
        s=solveKeplersEq(em/e,-e1/e,3);
        c=sqrt(1+s^2);
        s2=s^2/(c+1);
        r=q-e*s2/alpha;
        d=e*s/sqrt(abs(alpha));
        v=atan2(s*h*sqrt(abs(alpha)),-GM*s2-e1);
    end
end

%All orbits
u=omega+v;
VR=d/r;
VT=h/r;
end

function x=solveCubicEqGoodin(a,b,c)
%SOLVECUBICEQGOODING  This solves the equation a*x^3+3*b*x-2*c=0, where
%                     a>=0 and b^3+a*c^2>=0 for the real root of x.
%                     Additionally, if a and b are both zero, then zero is
%                     generated in lieu of an indeterminate solution.
%
%INPUTS: a,b,c The coefficients in the equation a*x^3+3*b*x-2*c=0 where
%              a>=0 and b^3+a*c^2>=0.
%
%OUTPUTS: x The one real solution to a*x^3+3*b*x-2*c=0.
%
%In the report
%R. H. Gooding and A. W. Odell, "The hyperbolic Kepler's equation,
%and the elliptic equations revisited," Royal Aerospace Executive,
%Procurement Executive, Ministry of Defence, Farnborough, Hants, United
%Kingdom, Tech. Rep. 369, Jul. 1989.
%the source of the algorithm is listed as being in 
%R. H. Gooding, "Solution of the hyperbolic Kepler's equation," Royal
%Aerospace Executive, Procurement Executive, Ministry of Defence,
%Farnborough, Hants, United Kingdom, Tech. Rep. 87042, 1987.
%However, the 1987 report appears to be unobtainable. Nonetheless, the
%solution is given in Equation 4.47 in
%G. Beutler, Methods of Celestial Mechanics: Physical, Mathematical and
%Numerical Principles. Berlin: Springer, 2005, vol. 1.
%The modification to return zero instead of infinity is consistent with how
%Gooding used the algorithm.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(a==0&&b==0||c==0)
    x=0;
else
    d=nthroot(sqrt(a)*abs(c)+sqrt(b^3+a*c^2),3)^2;
    x=2*c/(d+b+b^2/d);
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
