function [vMinStartEllip,tMinEllip,tMinParab]=orbVelDet2PtMinEng(r1Vec,r2Vec,GM)
%%ORBVELDET2PTMINENG Determine the time of the minimum energy elliptical
%                    trajectory between two orbital points and the initial
%                    velocity needed to achieve that trajectory. This is
%                    for the two-body problem assuming a Keplerian
%                    gravitational model in a generic inertial coordinate
%                    system. Also determine the time of the minimum energy
%                    parabolic trajectory, which is less than the time of
%                    the elliptical trajectory. Note that these are
%                    zero-revolution solutions. That is, the path does not
%                    complete a full revolution of the planet between the
%                    points. Also, the shortest way between the points is
%                    used for the elliptical solution (there are two ways
%                    to go around an ellipse).
%
%INPUTS: r1Vec,r2Vec 3X1 position vectors in (quasi)-inertial Cartesian
%                    coordinates where the massive body is at the origin.
%                    For example, Earth-centered inertial. The units are
%                    assumed meters per second.
%                 GM An optional value of the universal gravitational
%                    constant times the mass of the Earth. If omitted, the
%                    value Constants.WGS84GMWithAtmosphere is used. The
%                    units are m^3/sec^2.
%
%OUTPUTS: vMinStartEllip The velocity vector at r1Vec for the minimum
%                        energy elliptical trajectory from r1Vec to r2Vec.
%              tMinEllip The time needed to traverse the minimum energy
%                        elliptical trajectory from r1Vec to r2Vec.
%              tMinParab The time needed to traverse the minimum energy
%                        parabolic trajectory between r1Vec and r2Vec.
%
%The algorithm is described step-by-step in Chapter 7.6.1 of [1] and with a
%more detailed derivation in Chapter 3.3 of [2].
%
%REFERENCES:
%[1] D. A. Vallado and W. D. McClain, Fundamentals of Astrodynamics and
%    Applications, 4th ed. Hawthorne, CA: Microcosm press, 2013.
%[2] R. H. Battin, Astronautical Guidance. New York: McGraw Hill, 1964.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   GM=Constants.WGS84GMWithAtmosphere;
end

r1=norm(r1Vec);
r2=norm(r2Vec);

deltaV=angBetweenVecs(r1Vec,r2Vec);
cosDeltaV=cos(deltaV);
sinDeltaV=sin(deltaV);

%Chord length between the points.
c=sqrt(r1^2+r2^2-2*r1*r2*cosDeltaV);
%Semiparameter
s=(r1+r2+c)/2;

%Semi-major axis of the minimum-energy ellipse.
aMin=s/2;
pMin=r1*r2/c*(1-cosDeltaV);
alphae=pi;

betae=2*asin(sqrt((s-c)/s));

%In Chapter 3.3 of Battin and Chapter 7.6.1 of Vallado.
tMinEllip=sqrt(aMin^3/GM)*(alphae-(betae-sin(betae)));

%Equation 3.24 in Battin and 7-39 in Vallado
tMinParab=(1/3)*sqrt(2/GM)*(s^(3/2)-(s-c)^(3/2));

vMinStartEllip=sqrt(GM*pMin)/(r1*r2*sinDeltaV)*(r2Vec-(1-r2/pMin*(1-cosDeltaV))*r1Vec);
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
