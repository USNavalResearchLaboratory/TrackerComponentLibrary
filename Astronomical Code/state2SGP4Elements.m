function [SGP4Elements,didConverge]=state2SGP4Elements(stateVec,dragTerm,dragType,TTEpoch1,TTEpoch2,opsMode,gravityModel)
%%STATE2SGP4ELEMENTS Convert the state of a target orbiting the Earth given
%                    in terms of position and velocity in True Equator Mean
%                    Equinox (TEME) of date coordinate system into orbital
%                    elements used by the Simplified General Perturbations
%                    4 (SGP4)/ Simplified Deep Space Perturbations 4 (SDP4)
%                    propagator given by the propagateOrbitSGP4 function.
%                    This only work for elliptical orbits, not for orbits
%                    that are on escape trajectories (parabolic,
%                    hyperbolic).
%
%INPUTS: stateVec A 6X1 state vector consisting of 3D position and velocity
%                 components in the TEME coordinate system in meters and
%                 meters per second.
%        dragTerm The drag term associated with the target. This can either
%                 be the B* term used in the SGP4 propagator (which has
%                 units of inverse Earth radii) or the ballistic
%                 coefficient BC, defined as BC=Cd*Area/mass where Cd is
%                 the drag coefficient having units of kg/m^2. If this
%                 parameter is omitted, zero drag is assumed.
%        dragType Specified the type of the drag term. Possible values are
%                 0 (The default if omitted) The drag term is a ballistic
%                    coefficient.
%                 1 The drag term is the B* term used in the SGP4
%                   propagator.
% TTEpoch1, TTEpoch2 The epoch time of the SGP4Elements given in
%                 terrestrial time (TT) as a two-part Julian date. These
%                 are only required if 2*pi/SGP4Elements(6)>=225, which is
%                 when the deep space propagator is used. Otherwise, these
%                 inputs can be omitted or empty matrices can be passed.
%         opsMode An optional parameter specifying the orbital propagation
%                 model. Possible values are
%                 0 (The default if omitted) Use a model that is supposed
%                   to be similar to what the Air Force Space Command
%                   (AFSPC) uses when they publish orbital elements.
%                 1 Use a model that is supposed to be a slight improvement
%                   over what the AFSPC uses.
%    gravityModel An optional parameter specifying which gravitational
%                 model should be used. Since the gravitational models stem
%                 from ellipsoidal Earth models, the choice of the model
%                 also affects the radius of the Earth term used in the
%                 algorithm. Possible values are
%                 0 (The default if omitted). Use a gravitational model
%                   based on the WGS-72 reference ellipsoid. This appears
%                   to be what the AFSPC uses in their TLE sets.
%                 1 Use a gravitational model based on the WGS-84 reference
%                   ellipsoid.
%
%OUTPUTS: SGP4Elements The 7X1 vector of SGP4 orbital elements that
%                      correspond to the target state vector or an empty
%                      matrix if the algorithm was unable to produce an
%                      appropriate set of elements. Note that if the 
%                      algorithm does not converge, SGP4Elements might not
%                      be empty, but will be inaccurate. Also, a failure of
%                      the algorithm to produce elements does not mean that
%                      such elements do not exist. The seven elements are
%                      1) eccentricity
%                      2) inclination  (radians)
%                      3) argument of perigee/periapsis (radians)
%                      4) right ascension of the ascending node (radians)
%                      5) mean anomaly (radians)
%                      6) mean motion (radians per second [TT])
%                      7) BSTAR drag term. This pseudo ballistic
%                        coefficient has units of inverse Earth radii.
%                        Normally, a ballistic coefficient BC is mass per
%                        unit area divided by a drag coefficient. The BSTAR
%                        drag term is a scaled version of the normal
%                        ballistic coefficient. That is,
%                        BSTAR=BC*rho0/2 where rho0=2.461e-5 kg/m^2. The
%                        ballistic coefficient BC is Cd*Area/mass where Cd
%                        is the drag coefficient of the object.
%         didConverge A boolean value indicating whether the algorith
%                     converged to an accurate solution.
%
%The algorithm first converts the state vector into the set of orbital
%elements described in [1] via the state2OrbEls function. These orbital
%elements are then manipulated to obtain the orbital elements for the SGP4
%propagator under a simple Keplerian dynamic model. Newton's method is then
%used to refine the initial estimate for the SGP4 propagator models. The
%basic idea behind the technique is from [2], though the code provided
%therein is not used.
%
%The reference value of the atmospheric density in kg/m^2 used in the SGP4
%propagator as given in Appendix B of [3].
%
%The relationship between the semi-major axis and the mean motion is taken
%from [4]. That paper also uses a Newton's method approach to find SGP4
%orbital elements.
%
%The algorithm will not work if the state vector is a non-orbital
%trajectory. The algorithm is sensitive to the initialization value
%produced by Gooding's algorithm and can fail if the initialization is not
%sufficiently accurate.
%
%Note that this is NOT based on the official SGP4 orbital propagator used
%by the U.S. Air Force and cannot be assumed to be as reliable or produce
%identical results to the official propagator. Information on obtaining the
%U.S. Air Force's official propagator is given at
%http://www.afspc.af.mil/units/ASDA/
%
%REFERENCES:
%[1] R. H. Gooding, "On universal elements, and conversion procedures to
%    and from position and velocity," Celestial mechanics, vol. 44, no. 3,
%    pp. 283-298, 1988.
%[2] D. E. Andersen, "Computing mean orbital elements from a state vector,"
%    Master's thesis, Air Force Institute of Technology, Wright- Patterson
%    Air Force Base, OH, Dec. 1994.
%[3] D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, "Revisiting
%    spacetrack report # 3: Rev 2," in Proceedings of the AIAA/AAS
%    Astrodynamics Specialist Conference and Exhibit, Keystone, CO, 21-24
%    Aug. 2006.
%[4] D. A. Vallado and P. Crawford, "SGP4 orbit determination," in
%    Proceedings of the AIAA/AAS Astrodynamics Specialist Conference and
%    Exhibit, Honolulu, HI, 18-21 Aug. 2008.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    dragTerm=0;
end

if(nargin<3)
    dragType=0;
end

if(nargin<4)
    TTEpoch1=[];
    TTEpoch2=[];
end

if(nargin<6)
    opsMode=0;
end

if(nargin<7)
    gravityModel=0;
end

switch(dragType)
    case 0%The ballistic coefficient is given.


        rho0=2.461e-5;%kg/m^2

        %The BStar drag term in units of inverse Earth radii. The formula
        %is taken from http://celestrak.com/columns/v04n03/ as the units in
        %the aforementioned Vallado paper do not appear to be correct.
        BStarDrag=dragTerm*rho0/2;
    case 1%The BSTAR drag term is given directly.
        BStarDrag=dragTerm;
    otherwise
        error('Invalid drag type specified.')
end

if(gravityModel==0)
    GM=398600.8e9;%WGS-72 value used in the SGP4 code in km^3/s^2
else
    GM=398600.5e9;%WGS-84 value used in the SGP4 code in km^3/s^2
end


%Convert the state to Gooding's orbital elements under a simple Keplerian
%dynamic model.
orbEls=state2OrbEls(stateVec,0,GM);

%Get the components of Gooding's universal elements
alpha=orbEls(1);%alpha=GM/a where a is the semi-major axis in meters
q=orbEls(2);%q=a*(1-e)  where e is the eccentricity (unitless). This is the
%perifocal distance.
i=orbEls(3);%inclination in radians
Omega=orbEls(4);%longitude (right ascension) of the ascending node in radians
omega=orbEls(5);%argument of periapsis/perigee in radians
tau=orbEls(6);%the time at pericenter (seconds)

els2Est=zeros(6,1);%Allocate space.

e=(1-q*alpha/GM);
els2Est(1)=e;%eccentricity
els2Est(2)=i;%inclination
els2Est(3)=wrapRange(omega,0,2*pi);%argument of periapsis/perigee
els2Est(4)=wrapRange(Omega,0,2*pi);%right ascension of the ascending node
%The mean motion n is sqrt(GM/a^3). This is from Equation 6 in the Vallado/
%Crawford SGP4 orbit determination paper.
a=q/(1-e);
n=sqrt(GM/a^3);%mean motion in radians per second.
els2Est(5)=wrapRange(n*tau,0,2*pi);%The mean anomaly.
els2Est(6)=n;%The mean motion in radians per second

%The initial estimates of the SGP4 orbital elements are based on a
%Keplerian model. Now, use Newton's method to estimate the value under the
%SGP4 model. Numerical differentiation is used for the derivatives

%A function handle to convert the current estimate into a state vector.
SGP4Els2StateFunc=@(els2Est)propagateOrbitSGP4([els2Est;BStarDrag],0,TTEpoch1,TTEpoch2,opsMode,gravityModel);

didConverge=false;
for numIter=1:100
%Get the state for the current estimates of the orbital elements
    [xStateEst,errorState]=propagateOrbitSGP4([els2Est;BStarDrag],0,TTEpoch1,TTEpoch2,opsMode,gravityModel);

    %If any type of error occurred (other than indicating that the
    %satellite decayed), then stop.
    if(errorState~=0&&errorState~=6)
        SGP4Elements=[];
        return;
    end
    
    %Check for an invalid SGP4 estimate.
    if(any(~isfinite(xStateEst)))
        SGP4Elements=[];
        return;
    end
    
    %Check for convergence.
    %The position vector will never be zero (center of the Earth).
    relPosErr=abs(xStateEst(1:3)-stateVec(1:3))/norm(stateVec(1:3));
    %Except for the extrema of rectilinear orbits (in the TEME coordinate
    %system), which are probably not of interest, the velocity vector will
    %never be zero.
    relVelErr=abs(xStateEst(4:6)-stateVec(4:6))/norm(stateVec(4:6));
    if(all(relPosErr<=1e-15)&&all(relVelErr<=1e-15))
        didConverge=true;
        break;
    end
    
%Numerically compute the Jacobian matrix using numerical differentiation.
%This requires some sort of epsilon offset to each of the terms. A value of
%1e-5 times the value of the term, limited to certain ad-hoc minimum values
%is used.
    epsilon=max(1e-5*abs(els2Est),[1e-12*ones(5,1);1e-15]);
    J=numDiff(els2Est,SGP4Els2StateFunc,6,2,epsilon);
    
    %If a problem finding the numeric derivative arose.
    if(isempty(J))
        SGP4Elements=[];
        return;
    end
    
%Update using a step of Newton's method. This is essentially just
%els2Est=els2Est+inv(J)*(stateVec-xStateEst);
%However, we will use singular value decomposition as it is less prone to
%numerical problems than the inv function and it avoids Matlab producing
%a warning if the conditioning is bad.
    [U,S,V]=svd(J);
%We want a matrix whose diagonals are the inverse of the diagonals of S.
%However, if any of the singular values are zero/ near zero, then there
%will e bad results. By setting the inverses of zero singular values to
%zero, we are essentially taking a matrix inverse that is providing a
%least squares fit when the output is not in the domain of the singular
%matrix.
    s=diag(S);
    sel=s/s(1)<eps;
    sInv=1./s;
    sInv(sel)=0;
 
    %Otherwise, update the estimate.
    els2Est=els2Est+V*diag(sInv)*(U'*(stateVec-xStateEst));
    
    %If the eccentricity is negative, make it zero or the iteration will
    %fail.
    els2Est(1)=min(max(els2Est(1),0),1);
    %Keep all of the other values in the minimum range. Otherwise,
    %divergence is likely.
    els2Est(2)=wrapRange(els2Est(2),0,pi,true);
    els2Est(3)=wrapRange(els2Est(3),0,2*pi);
    els2Est(4)=wrapRange(els2Est(4),0,2*pi);
    els2Est(5)=wrapRange(els2Est(5),0,2*pi);
end

SGP4Elements=zeros(7,1);%Allocate space for the return value.
SGP4Elements(1:6)=els2Est;
SGP4Elements(7)=BStarDrag;
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
