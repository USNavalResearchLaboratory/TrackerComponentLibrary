function f=altEllipsParam2Flattening(omega,a,C20Bar,GM)
%%ALTELLIPSPARAM2FLATTENING Some reference ellipsoids, such as that of the
%                           EGM2008, are defined using a second degree
%                           spherical harmonic gravitational coefficient
%                           and not a flattening factor. This function
%                           takes the parameters of the reference ellipsoid
%                           and returns a flattening factor, so that a more
%                           standard parameterization can be used.
%
%INPUTS: omega The average rotation rate of the planet, generally having
%              units of radians per second.
%            a The semi-major axis of the reference ellipsoid, generally
%              having units of meters.
%       C20Bar The fully normalized second degree zonal coefficient in a
%              spherical harmonic representation of the ellipse's
%              potential. In older models, a J2 term is used. The J2 term
%              is related to C20Bar by J2=-C20Bar*sqrt(5).
%           GM The universal gravitational constant times the mass of the
%              planet, generally having units of m^3/s^2.
%
%OUTPUTS:  f The flattening factor of the reference ellipsoid. 
%
%The conversion is an iterative method taken from [1].
%
%REFERENCES:
%[1] H. Moritz, "Geodetic reference system 1980," Bulletin géodésique,
%    vol. 54, no. 3, pp. 395-405, 1980.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

C2=C20Bar*sqrt(5);
J2=-C2;

e=0.8;
eOld=Inf;
numIter=0;
while(abs(eOld-e)>10^3*eps(e))
    eTilde=e/sqrt(1-e^2);%The second eccentricity.
    q0=(1/2)*((1+3/eTilde^2)*atan(eTilde)-3/eTilde);
    eOld=e;
    e=sqrt(3*J2+(4/15)*(omega^2*a^3/GM)*(e^3/(2*q0)));
    numIter=numIter+1;
end

f=1-sqrt(1-e^2);
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
