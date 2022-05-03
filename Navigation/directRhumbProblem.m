function latLonEnd=directRhumbProblem(latLonStart,azimuth,dist,algorithm,a,f)
%%DIRECTRHUMBPROBLEM Given a starting location on the surface of the
%                 reference ellipsoid as well as a constant heading a
%                 distance to travel, determine the location of a vessel
%                 traveling at the constant heading after traveling the
%                 distance. A constant heading trajectory is a rhumb line
%                 (loxodrome) and is usually not the shortest path between
%                 two points.
%
%INPUTS: latLonStart A 2X1 vector of the starting ellipsoidal latitude
%                    (North) and longitude (East) in radians. This cannot
%                    be a pole (cannot be +/-pi/2).
%            azimuth The constant heading that is to be traveled in
%                    radians East of North.
%               dist The NX1 or 1XN set of distances on the ellipsoid-
%                    approximated Earth starting from latLonStart where one
%                    wishes to find the stopping point. The distances must
%                    be less than the rhumb distance to reach the pole or
%                    invalid results will be obtained.
%          algorithm An optional parameter specifying the algorithm to use.
%                    Possible values are:
%                    0 Use the algorithm of [1] coupled with a formula
%                      using isometric latitudes, which are described in
%                      Chapter 3 of [2] to get the azimuth angle.
%                    1 Use the explicit approximation of [4].
%                  a The semi-major axis of the reference ellipsoid. If
%                    this argument is omitted or an empty matrix is passed,
%                    the value in Constants.WGS84SemiMajorAxis is used.
%                  f The flattening factor of the reference ellipsoid. If
%                    this argument is omitted or an empty matrix is passed,
%                    the value in Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonEnd The ellipsoidal latitude and longitude that one will be
%                    after starting at latLonStart and traveling a distance
%                    of dist at a constant heading of azimuth.
%
%Algorithm 0 is taken from [1]. However, a formula using isometric
%latitudes, which are described in Chapter 3 of [2] to get the azimuth
%angle was used, because it is simpler. The formula is also explicitly
%mentioned in Equation 2 of [3]. However, the expression for computing the
%distance from [3] is only for a sphere, not for an ellipsoid, which is why
%the Carlton-Wippern distance computation using an incomplete elliptic
%integral of the second kind is preferred.
%
%Both techniques have singularities with azimuth values near +/-pi/2, and
%different approximations are used. Algorithm 0 is slower than 1 but is
%more accurate.
%
%EXAMPLE:
%We show that after solving
%the indirect rhumb problem, the trajectory obtained leads to the
%desired enpoint when given to directRhumbProb. We also show the error
%incurred by using algorithm 1, which is faster.
% latStart=degMinSec2Rad(37,47.5);
% lonStart=degMinSec2Rad(-122,-27.8);
% latEnd=degMinSec2Rad(-33,-51.7);
% lonEnd=degMinSec2Rad(151,12.7);
% latLonStart=[latStart;lonStart];
% latLonEnd=[latEnd;lonEnd];
% 
% [azStart,dist]=indirectRhumbProblem(latLonStart,latLonEnd);
% latLonPoint0=directRhumbProblem(latLonStart,azStart,dist,0);
% latLonPoint1=directRhumbProblem(latLonStart,azStart,dist,1);
% max(abs(wrapRange(latLonPoint0-latLonEnd,-pi,pi)))
% max(abs(wrapRange(latLonPoint1-latLonEnd,-pi,pi)))
%One will see that the value from algorithm 0 is pretty much exact, whereas
%the value from algorithm 1, the approximationm, is close but still notably
%different.
%
%REFERENCES:
%[1] K. C. Carlton-Wippern, "On loxodromic navigation," Journal of 
%    Navigation, vol. 45, no. 2, pp. 292-297, May 1992.
%[2] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[3] J. Alexander, "Loxodromes: A rhumb way to go," Mathematics Magazine,
%    vol. 77, no. 5, pp. 349-356, Dec. 2004.
%[4] G. H. Kaplan, "Practical Sailing Formulas for Rhumb-Line Tracks on an
%    Oblate Earth," Navigation: Journal of the Institute of Navigation,
%    vol. 42, no. 2, pp. 313-326 Summer 1995.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(algorithm))
    algorithm=0;
end

switch(algorithm)
    case 0
        dist=dist(:);
        N=length(dist);
        
        latLonEnd=zeros(2,N);
        
        for k=1:N
            latLonEnd(:,k)=directRhumbProbCarlton(latLonStart,azimuth,dist(k),a,f);
        end
    case 1
        latLonEnd=directRhumbProblemKaplan(latLonStart,azimuth,dist,a,f);
    otherwise
        error('Unknown algorithm specified.')
end
end

function latLonEnd=directRhumbProbCarlton(latLonStart,azimuth,dist,a,f)
%%DURECTRHUMBPROBCARLTON Solve the direct rhumb problem using the iterative
%             approximation in [1]. A formula using isometric latitudes,
%             which are described in Chapter 3 of [2] to get the azimuth
%             angle was used, because it is simpler than the formula of [1].
%
%REFERENCES:
%[1] K. C. Carlton-Wippern, "On loxodromic navigation," Journal of 
%    Navigation, vol. 45, no. 2, pp. 292-297, May 1992.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

latStart=latLonStart(1);
lonStart=latLonStart(2);

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

%Convert the ellipsoidal latitudes to a reduced co-latitude. A co-
%latitude is pi/2 minus the latitude.
nu1=pi/2-ellipsLat2ReducedLat(latStart,f);

%Make sure that the azimuth is in the range of -pi to pi.
azimuth=wrapRange(azimuth,-pi,pi,false);

%If the azimuth is too close to due East/West, then use the solution for
%due East/West to avoid numerical precision issues. Otherwise, use the
%solution for other headings in general.
if(abs(abs(azimuth)-pi/2)>1e-8)
    %The solution for general headings.

    %The following sets up and performs the iteration from Equation 13,
    %which inverts the incomplete elliptic function of the second kind
    %using Newton's method.
    nu2=nu1;
    if(-pi/2<azimuth&&azimuth<pi/2)
        signVal=1;
    else
        signVal=-1;
    end

    %It seems to converge in under 6 iterations
    numIter=6;
    for curIter=1:numIter
        num=abs(ellipIntInc2Kind(nu2,e^2)-ellipIntInc2Kind(nu1,e^2))-dist*abs(cos(azimuth))/a;
        denom=sqrt(1-e^2*sin(nu2)^2);

        nu2=nu2+signVal*num/denom;
    end

    %Convert the reduced co-latitude of the ending point to an ellipsoidal
    %latitude. 
    latEnd=reducedLat2EllipsLat(pi/2-nu2,f);

    %Next, we want to find the ending longitude. This could be done using
    %equation 11 in the Carlton-Wippern paper. However, we are going to do
    %it using isometric latitudes, since it is simpler.
    lonEnd=lonStart+tan(azimuth)*(ellipsLat2IsoLat(latEnd,f)-ellipsLat2IsoLat(latStart,f));

    %However, if the trajectory is too close to 0/ 180 degrees, then the
    %ellipsLat2IsoLat function might not be finite, in which case the
    %difference will be a NaN, rather than 0. Thus, we will check for that
    %case. In such an instance, the ending longitude is just set to the
    %starting latitude.
    if(~isfinite(lonEnd))
       lonEnd=lonStart; 
    end
else
    %The solution for East-West headings. This inverts equation 14b to get
    %the ending longitude.

    %If it is going East (positive direction)
    if(abs(azimuth-pi/2)<abs(azimuth+pi/2))
        signVal=1;
    else
        signVal=-1;
    end
    lonEnd=lonStart+signVal*dist/(a*abs(sin(nu1)));

    %If one tries to have an East-West trajectory at the pole, then the
    %denominator above will be 0. This keeps it from returning a NaN.
    if(~isfinite(lonEnd))
        lonEnd=lonStart; 
    end

    latEnd=latStart;
end
%Make sure that the longitude is within the proper range.
latLonEnd=[latEnd;wrapRange(lonEnd,-pi,pi,false)];

end

function latLonEnd=directRhumbProblemKaplan(latLonStart,azimuth,dist,a,f)
%%DIRECTRHUMBPROBKAPLAN Solve the direct rhumb problem using the explcit
%                       approximation in [1].
%
%This is a non-iterative approximation that is significantly better than a
%spherical Earth approximation.
%
%REFERENCES:
%[1] G. H. Kaplan, "Practical Sailing Formulas for Rhumb-Line Tracks on an
%    Oblate Earth," Navigation: Journal of the Institute of Navigation,
%    vol. 42, no. 2, pp. 313-326 Summer 1995.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<4||isempty(f))
    a=Constants.WGS84SemiMajorAxis;
end

phi0=latLonStart(1);
lambda0=latLonStart(2);

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

eSinPhi0=e*sin(phi0);
if(abs(abs(azimuth)-pi/2)>1e-8)
    %If it is not nearly a +/- 90 degree heading.
    
    M0=1/sqrt(1-eSinPhi0^2)^3;
    M=a*(1-e^2)*M0;%Meridian radius of curvature.
    phiPrime=phi0+dist*cos(azimuth)/M;

    %Equation 13 in [1].
    phi=phi0+M0*((1-(3/4)*e^2)*(phiPrime-phi0)+(3/8)*e^2*(sin(2*phiPrime)-sin(2*phi0)));

    eSinPhi=e*sin(phi);
    %Equation 16 in [1].
    lambda=lambda0+tan(azimuth)*(log(tan(pi/4+phi/2))+...
           (e/2)*log((1-eSinPhi)./(1+eSinPhi))-...
           log(tan(pi/4+phi0/2))-(e/2)*log((1-eSinPhi0)/(1+eSinPhi0)));
else
    M=a*(1-e^2)/sqrt(1-eSinPhi0^2)^3;
    phiPrime=phi0+dist*cos(azimuth)/M;
    
    %Prime vertical radius of curvature.
    N0=a/sqrt(1-eSinPhi0^2);
        
    %These expressions are in Equation 17 in [1].
    phi=phiPrime;
    lambda=lambda0+dist*sin(azimuth)./(N0*cos((1/2)*(phi0+phiPrime)));
    
    %The zero azimuth, case where the denominatot above is zero.
    sel=~isfinite(lambda);
    lambda(sel)=lambda0(sel);
end

latLonEnd=[phi;lambda];
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
