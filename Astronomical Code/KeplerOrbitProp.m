function statePred=KeplerOrbitProp(stateOrig,deltaT,GM,epsVal)
%%KEPLERORBITPROP Propagate a state consisting of position and velocity
%                 forward while under a simple, two-mody ballistic,
%                 Keplerian dynamic model.
%
%INPUTS: stateOrig A 6XN matrix of N states to propagate forward in time by
%                  the same interval. A state consists of poistion and
%                  velocity in a (quasi)-inertial coordinate system with
%                  the mass producing body at the origin. Units are meters
%                  and meter per second.
%           deltaT The time duration in seconds by which the state should
%                  be predicted. Positive for forward, negative for
%                  backward.
%               GM An optional value of the universal gravitational
%                  constant times the mass of the Earth. If omitted, the
%                  value Constants.WGS84GMWithAtmosphere is used. The
%                  units are m^3/sec^2.
%           epsVal An optional parameter that is used to determine the
%                  cutoff for when e is zero, one, when the angular
%                  momentum is zero, and when the specific energy is zero.
%                  The default if omitted is 1e-15.
%
%OUTPUTS: statePred The state predicted forward deltaT seconds. Note that
%                   when predicting rectinlinear motion models forward, the
%                   bahaviour as a target passes the origin (the
%                   mass-producing point object) is not well defined in
%                   terms of the physics.
%
%The algorithm is an implementation of [1].
%
%REFERENCES:
%[1] D. Condurache and V. Martinusi "A complete closed form vectorial
%    solution to the Kepler problem," vol. 42, no. 5, Oct. 2007, pp.
%    465-476.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    epsVal=1e-15;
end

if(nargin<3)
    GM=Constants.WGS84GMWithAtmosphere; 
end

numStates=size(stateOrig,2);
statePred=zeros(6,numStates);
for curState=1:numStates
    r0=stateOrig(1:3,curState);
    v0=stateOrig(4:6,curState);

    %Specific angular momentum, Equation 2
    Omega=cross(r0,v0);

    %Specific energy, Equation 3
    h=(1/2)*norm(v0)^2-GM/norm(r0);

    %Eccentricity vector, Equation 4
    eVec=cross(v0,Omega)/GM-r0/norm(r0);
    e=norm(eVec);

    %The mean motion, Equation 24
    n=(2*abs(h))^(3/2)/GM;
    if(norm(Omega)>epsVal)%Non-rectilinear motion, Section 4.1
        if(abs(e)<=epsVal)%e=0 case, Circular trajectory
            %Equation 42
            r=cos(n*deltaT)*r0+sin(n*deltaT)*v0/n;
            %Equation 44
            v=-n*sin(n*deltaT)*r0+cos(n*deltaT)*v0;
        elseif(abs(e-1)<=epsVal)%Parabolic trajectory
            %Equation 50
            tau0=dot(r0,v0)/GM;

            %The semi-latus rectum, mentioned after Equation 45
            p=norm(Omega)^2/GM;

            %Equation 51, taking t0=0.
            tP=0-(tau0/2)*(p+(GM/3)*tau0^2);

            deltaTPt=tP-deltaT;
            %Equation 52
            tauT=1/GM^(1/3)*(nthroot(3*deltaTPt+sqrt(9*deltaTPt^2+p^3/GM),3)+...
                nthroot(3*deltaTPt-sqrt(9*deltaTPt^2+p^3/GM),3));

            %Equation 45
            r=(1/2)*(p-GM*tauT^2)*eVec+tauT*cross(Omega,eVec);

            %Equation 47
            v=(2/(p+GM*tauT^2))*(-GM*tauT*eVec+cross(Omega,eVec));
        elseif(e<1)%Elliptical trajectory
            %Equation 24
            a=GM/(2*e*abs(h))*eVec;
            b=1/(e*sqrt(2*abs(h)))*cross(Omega,eVec);

            %Equation 32, for the elliptical case
            cosE0=(1/e)*(1-n*norm(r0)/sqrt(2*abs(h)));
            sinE0=n/(2*abs(h)*e)*dot(r0,v0);
            E0=atan2(sinE0,cosE0);
            %The solution to Equation 34
            E=solveKeplersEq(n*deltaT+E0-e*sin(E0),e);

            %Equation 26 for the propagated position
            r=a*(cos(E)-e)+b*sin(E);

            %Equation 29 for the propagated velocity
            v=(n/(1-e*cos(E)))*(-a*sin(E)+b*cos(E));
        else%e>1 || h>0, the hyperbolic trajectory
            %Equation 24
            a=GM/(2*e*abs(h))*eVec;
            b=1/(e*sqrt(2*norm(h)))*cross(Omega,eVec);

            %Equation 62, for the hyperbolic case
            E0=asinh(n*dot(r0,v0)/(2*e*h));
            %The solution to equation 64
            E=solveKeplersEq(n*deltaT+e*sinh(E0)-E0,e);

            %Equation 57
            r=a*(e-cosh(E))+b*sinh(E);

            %Equation 60
            v=(n/(e*cosh(E)-1))*(-sinh(E)*a+cosh(E)*b);
        end
    else%Omega=0, rectlinear trajectory, Section 4.2
        if(abs(h)<=epsVal)%Zero specific energy
            %Equation 84
            tau0=dot(r0,v0)/GM;
            %Equation 85 taking t0=0;
            tP=0-(GM/6)*tau0^3;
            tauT=nthroot(6*(deltaT-tP)/GM,3);

            %Equation 80
            r=(GM/2)*tauT^2*r0/norm(r0);
            %Equation 82
            v=(2/tauT)*r0/norm(r0);
        elseif(h<0)%Negative specific energy
            %From before Equation 70
            a=GM/(2*norm(h));

            %Equation 73
            cosE0=1-n*norm(r0)/sqrt(2*abs(h));
            sinE0=(n/(2*abs(h)))*dot(r0,v0);
            E0=atan2(sinE0,cosE0);

            %Equation 74 for t0=0
            tP=0-(1/n)*(E0-sin(E0));

            %Solving Equation 72
            E=solveKeplersEq(n*(deltaT-tP),e);

            %Equation 69
            r=a*(1-cos(E))*r0/norm(r0);
            %Equation 71
            v=(n*a*sin(E)/(1-cos(E)))*(r0/norm(r0));
        else%h>0, positive specific energy
            temp=n*dot(r0,v0)/(2*h);
            %Equation 96 using t0=0.
            tP=0-(1/n)*(temp-asinh(temp));
            %Solving Equation 94
            E=solveKeplersEq(n*(deltaT-tP),1,1);

            %Given before Equation 92.
            a=GM/(2*h);

            %Equation 91
            r=a*(cosh(E)-1)*r0/norm(r0);

            %Equation 93
            v=n*a*sinh(E)/(cosh(E)-1)*r0/norm(r0);
        end
    end
    statePred(:,curState)=[r;v];
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
