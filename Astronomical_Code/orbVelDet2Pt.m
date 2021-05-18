function [V1Comp,V2Comp,alpha]=orbVelDet2Pt(r1Val,r2Val,deltaT,angParam,isECEF,longWay,GM,EarthRotRate)
%%ORBVELDET2PT Solve Lambert's problem. That is, given two position
%            measurements of a ballistic (orbital) object and the time
%            difference of the measurements, determine the velocities of
%            the object at each of the two points assuming a Keplerian
%            gravitational model with just two bodies (the Earth (or Sun or
%            other planet), and the object (which is not massive). This
%            function can be useful to help start a track on a ballistic
%            target and it can be useful in designing transfer orbits. The
%            function can work in a generic inertial coordinate system, or
%            it can be used with a generic Earth-centered Earth-fixed
%            (ECEF) coordinate system. That is, rather than being
%            rigorously-defined in terms of precession, nutation, Earth
%            orientation parameters, etc, the inputs are given in a generic
%            ECEF coordinate system that rotates at a rate of EarthRotRate
%            radians per second about the z-axis with respect to a generic
%            Earth centered inertial (ECI) coordinate system.
%
%INPUTS: The function can be called either as
%        orbVelDet2Pt(r1Vec,r2Vec,deltaT,m,longWay,GM), where one passes
%        the vectors r1Vec,r2Vec to the points and m is the number of
%        complete revolutions between the points, or one can use 
%        orbVelDet2Pt(r1,r2,deltaT,theta,longWay,GM), where r1 and r2 are
%        scalar values. The type of input is inferred from the dimensions
%        of r1Val.
%        r1Val,r2Val Depending on the calling method these can be are a 3XN
%                    set of 3D position vectors, for N instances of
%                    Lambert's that are to be solved, in (quasi)-inertial
%                    Cartesian coordinates  where the massive body is at
%                    the origin. For example, a generic ECI system (or a
%                    set of rotating ECEF coordinates if isECEF=true). The
%                    units are assumed meters. Alternatively, they can just
%                    be a 1XN set of ranges to the target if non-Cartesian
%                    outputs are desired.
%             deltaT The positive time in seconds between when r1Val and
%                    r2Val are measured.
%           angParam If r1Val and r2Val are just ranges, then angParam is
%                    the angle between the two full position vectors.
%                    angParam can be grater than 2*pi to indicate that a
%                    full orbital revolution was made between the points.
%                    It can also be greater in magnitude than pi/2 if the
%                    long way around the ellipse between the points is
%                    taken. On the other hand, if r1Val and r2Val are full
%                    3D vectors, then this is m, the number of complete
%                    orbital revolutions made for the transfer. For m>0,
%                    there might not always be a solution to Lambert's
%                    problem or there might be two solutions. If omitted, a
%                    value of 0 is used.
%             isECEF An optional boolean value indicating whether the r1Val
%                    and r2Val are given in a generic ECEF coordinate
%                    system and thus the outputs V1Comp and V2Comp are
%                    given in ECEF coordinates. The default if omitted is
%                    false. This must be false if r1Val and r2Val are just
%                    ranges. If true, then r1Val and r2Val must be given in
%                    units of meters and deltaT in units of seconds.
%            longWay This is only used if r1Val and r2Val are vectors (An
%                    empty matrix can be passed otherwise). longWay is an
%                    An optional boolean value indicating whether the
%                    orbital path goes the long way around during the time.
%                    For example, if r1Val and r2Val point 1 degree apart,
%                    the target could have moved 1 degree in angle during
%                    the period deltaT (longWay=false) or it could have
%                    moved 359 degrees (longWay=true). The default if
%                    omitted is false.
%                 GM An optional value of the universal gravitational
%                    constant times the mass of the Earth. If omitted, the
%                    value Constants.EGM2008GM is used. The units are
%                    m^3/sec^2.
%      EarthRotRate The rotation rate of the Earth in radians per second.
%                   This parameter is only used if isECEF=true. If this
%                   parameter is omitted or an empty matrix is passed, then
%                   Constants.EGM2008EarthRotationRate is used.
%
%OUTPUTS: V1Comp If r1Val and r2Val scalar ranges, then V1Comp is a
%                2XnumSol vector containing the radial and transverse
%                components of the velocity at the starting point assuming
%                a Keplerian dynamic model (i.e. no drag, higher-order
%                gravitational terms, etc.) V1Comp=[VR1;VT1];. If multiple
%                values are passed for r2Val and r2Val, then V1Comp is a
%                cell array of the solutions for all of the passed values.
%                A unit vector in the transverse direction would be 
%                cross(cross(r1Vec,r2Vec)/norm(cross(r1Vec,r2Vec)),r1Vec/norm(r1Vec));
%                if r1Vec and r2Vec are full 3D vectors to the starting and
%                ending points. On the other hand, r1Val and r2Val are
%                vectors, then V1Comp is a full 3X1 Cartesian unit vector
%                for the velocity of the target at the starting point. If a
%                solution does not exist, then empty matrices will be
%                returned. When m=0, a solution should always exist unless
%                finite precision problems arise. When m>0, a solution
%                might not exist. If m>0, two solutions might exist.
%         V2Comp The same as V1Comp, except at the ending point.
%          alpha A numSolX1 vector containing values of GM/a, where a is
%                the semi-major axis of the Keplerian orbit. If multiple
%                values are passed for r1Vec and r2Vec, then alpha is a
%                cell array holding the solutions for all of the values.
%                alpha can be useful when converting to orbital elements
%                and the trajectory is close to parabolic as the value
%                found here is not as susceptible to finite precision
%                errors.
%
%Much of the algorithm is taken from the report [1] and the paper [2],
%which are extensions of [3]. Note that [2] does not provide as good an
%explanation of how the algorithm works as [1].
%
%As an example for testing the algorithm, one can use the example
%parameters from Example 7-5 in Chapter 7.6 of [4], which are
% r1Vec=10^3*[15945.34;0;0];
% r2Vec=10^3*[12214.83899;10249.46731;0];
% m=0;
% deltaT=76*60;
% [VStart3D,VEnd3D]=orbVelDet2Pt(r1Vec,r2Vec,deltaT,m)
%For a vector solution and then
% r1=norm(r1Vec);
% r2=norm(r2Vec);
% deltaTheta=angBetweenVecs(r1Vec,r2Vec)+2*pi*m;
% [VStart2D,VEnd2D]=orbVelDet2Pt(r1,r2,deltaT,deltaTheta)
%for a solution in terms of components.
%
%The correction for ECEF coordinates just comes from the basic physics of
%converting between rotating and non-rotating coordinate systems, which is
%mentioned in Appendix A of [5].
%
%REFERENCES:
%[1] R. H. Gooding, "On the solution of Lambert's orbital boundary-value
%    problem," Royal Aerospace Executive, Procurement Executive, Ministry
%    of Defence, Farnborough, Hants, United Kingdom, Tech. Rep. 88027,
%    Apr. 1988.
%[2] R. H. Gooding, "A procedure for the solution of Lambert's orbital
%    boundary-value problem," Celetial Mechanics and Dynamical Astronomy,
%    vol. 48, no. 2, pp. 145-165, 1990.
%[3] E. R. Lancaster and R. C. Blanchard, "A unified form of Lambert's
%    theorem," National Aeronautics and Space Administration, Goddard
%    Space Flight Center, Greenbelt, MD, Tech. Rep. TN D-5368, Sep. 1969.
%[4] D. A. Vallado and W. D. McClain, Fundamentals of Astrodynamics
%    and Applications, 4th ed. Hawthorne, CA: Microcosm press, 2013.
%[5] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(r1Val,2);%The number of problems to solve.

if(nargin<8||isempty(EarthRotRate))
    EarthRotRate=Constants.EGM2008EarthRotationRate;
end

if(nargin<7||isempty(GM))
    GM=Constants.EGM2008GM;
end

if(nargin<6||isempty(longWay))
    longWay=false;
end

if(nargin<5||isempty(isECEF))
    isECEF=false;
end

if(nargin<4||isempty(angParam))
    m=0;
else
    m=angParam;
end
    
%If one value is provided, then return that directly.
if(N==1)
    [V1Comp,V2Comp,alpha]=solve1LambertProblem(r1Val,r2Val,deltaT,m,isECEF,longWay,GM,EarthRotRate);
    return;
end

%If multiple values are provided, then return all the results in a cell
%array.
V1Comp=cell(N,1);
V2Comp=cell(N,1);
alpha=cell(N,1);
for curProb=1:N
    [V1Comp{curProb},V2Comp{curProb},alpha{curProb}]=solve1LambertProblem(r1Val(:,curProb),r2Val(:,curProb),deltaT,m,isECEF,longWay,GM,EarthRotRate);
end
end

function [V1Comp,V2Comp,alpha]=solve1LambertProblem(r1Val,r2Val,deltaT,m,isECEF,longWay,GM,omega)
rDim=size(r1Val,1);
usingVecInputs=(rDim==3);
%The error tolerance for the convergence of x.
epsVal=eps;

if(~usingVecInputs)%If ranges were passed
    r1=r1Val;
    r2=r2Val;
    theta=m;
    
    if(isECEF)
        error('the ECEF option can only be used when using vector inputs for r1Val and r2Val')
    end
else
    r1=norm(r1Val);
    r2=norm(r2Val);
    
    if(isECEF)%Rotate r2 into the ECI position
        %omega is the rotation rate of the Earth.
        Omega=[0;0;1]*omega;%The z-axis is the axis of rotation.
        theta=deltaT*omega;%The angle rotated between detections.
        ECEF2ECI2=axisAng2RotMat([0;0;1],theta);%Rotation matrix.
        %Rotate into ECI coordinates consistent with the initial location
        %of r1Val.
        r2Val=ECEF2ECI2*r2Val;
    end
    
    %The short-way angle between the vectors.
    theta=angBetweenVecs(r1Val,r2Val)+m*2*pi;
    
    %If the long direction around the orbit is desired, adjust the angle and
    %make it negative so that the final velocity has the correct sign.
    if(longWay==true)
        theta=theta-2*pi;
    end
end

[V1Comp,V2Comp,alpha]=VALAMB(GM,r1,r2,theta,deltaT,epsVal);

%If full vector inputs are given, then provide full vector outputs.
if(usingVecInputs)
    %The number of solutions found.
    numSol=size(V1Comp,2);
    V1=zeros(3,numSol);
    V2=zeros(3,numSol);
    for curSol=1:numSol
        %To get the 3D vectors, we need to get the radial and tangent vectors.
        uR1=r1Val/r1;%Unit radial vector for r1.
        uR2=r2Val/r2;%Unit radial vector for r2.

        %A vector normal to the plane in which r1Vec and r2vec reside.
        rNormal=cross(r1Val,r2Val);
        uNormal=rNormal/norm(rNormal);

        %Find the tangent vectors.
        uT1 = cross(uNormal,uR1);%Unit vectors in the tangential-directions
        uT2 = cross(uNormal,uR2);

        %If any of the values in the unit vectors are not finite, then the
        %vectors are parallel (the trajectory is radial) and the tangential
        %components of V1Comp and V2Comp should be zero anyway, so just put in
        %the proper values. Note, however, that the results are correct for
        %parts of a 1D trajectory, which can be viewed as a collapsed ellipse.
        %Thus, the velocity at 2 might be after PASSING point 2 and then going
        %back (up/down) on the way back unless deltaT is the minimum energy
        %time.
        if(any(~isfinite(uT1)))
            V1(:,curSol)=V1Comp(1,curSol)*uR1;
        else%If the trajectory is not radial.
            V1(:,curSol)=V1Comp(1,curSol)*uR1+V1Comp(2,curSol)*uT1;
        end
        
        if(any(~isfinite(uT2)))
            V2(:,curSol)=V2Comp(1,curSol)*uR2;
        else%If the trajectory is not radial.
            V2(:,curSol)=V2Comp(1,curSol)*uR2+V2Comp(2,curSol)*uT2;
        end
    end
    
    V1Comp=V1;
    V2Comp=V2;

    if(isECEF)
        %The solution must be rotated back into ECEF coordinates from ECI.
        V1Comp=V1Comp-cross(Omega,r1Val);
        V2Comp=ECEF2ECI2'*(V2Comp-cross(Omega,r2Val));
    end
end
end

function [V1Vec,V2Vec,alpha]=VALAMB(GM,r1,r2,theta,deltaT,epsVal)
%%VALAMB  This function is taken from Appendix E of Gooding's report. This
%         solves Lambert's problem in the plane in which the vectors
%         reside. Thus, the problem is solved given scalar distances r1 and
%         r2, as well as the angle between them, theta, which can be
%         greater than 2*pi if a full revolution is made between the
%         vectors.
%
%INPUTS: GM   The universal gravitational constant times the mass of the
%             gravitating body (Earth). The units are The units are
%             m^3/sec^2.
%       r1,r2 The scalar distances from the center of the massive body to
%             the object at two times (meters).
%       theta The angle in radians between the vectors that r1 and r2 are
%             associated with.
%      deltaT The time in seconds between observations.
%      epsVal A small value to determine whether things are zero within
%             numerical precision bounds. For example, one can use eps.
%
%OUTPUTS: V1Vec A 2XN vector of the radial and transverse components of the
%               starting velocit(ies).
%         V2Vec A 2XN vector of the radial and transverse components of the
%               ending velocit(ies).
%         alpha The orbital element alpha=GM/a associated with the
%               trajectories, where a if the semi-major axis of the orbit.
%
%The original function is given in Appendix E of
%R. H. Gooding, "On the solution of Lambert's orbital boundary-value
%problem," Royal Aerospace Executive, Procurement Executive, Ministry
%of Defence, Farnborough, Hants, United Kingdom, Tech. Rep. 88027,
%Apr. 1988.
%
%The name VALAMB, with the extra A reflects the fact that this function has
%been modified to also return the value A,which is required when this is
%used as a subroutine in the algorithm 
%R. H. Gooding, "A new procedure for orbit determination based on three
%lines of sight (angles only)," Royal Aerospace Executive, Procurement
%Executive, Ministry of Defence, Farnborough, Hants, United Kingdom,
%Tech. Rep. 93004, Apr. 1993.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

    m=fix(theta/(2*pi));
    thr2=theta/2-m*pi;
    deltaR=r1-r2;
    r1r2=r1*r2;
    r1r2th=4*r1r2*sin(thr2)^2;
    
    csq=deltaR^2+r1r2th;
    
    c=sqrt(csq);
    s=(r1+r2+c)/2;
    gms=sqrt(GM*s/2);%This is gamma in the text of the report
    
    qsqfm1=c/s;
    q=sqrt(r1r2)*cos(thr2)/s;
    
    %Force q and qsqfm1 to have a valid magnitudes given possible finite
    %precision errors.; they cannot exceed 1.
    if(abs(q)>1)
       q=sign(q); 
    end
    if(abs(qsqfm1)>1)
        qsqfm1=sign(qsqfm1);
    end
    
    if(abs(c)>epsVal)
        rho=deltaR/c;
        sigmaVal=r1r2th/csq;
    else
        rho=0;
        sigmaVal=1;
    end
    
    %Equation 2 in the Gooding report.
    T=4*gms*deltaT/s^2;
    
    [x1,x2,N]=xlamb(m,q,qsqfm1,T,epsVal);
    
    %If there are no solutions, then these will come back empty
    V1Vec=zeros(2,N);
    V2Vec=zeros(2,N);
    alpha=zeros(N,1);

    for i=1:N
        if(i==1)
            x=x1; 
        else
            x=x2;
        end
        
        [~,qzminx,qzplx,zplqx]=tlamb(m,q,qsqfm1,x,-1);
        VT2=gms*zplqx*sqrt(sigmaVal);
        VR1=gms*(qzminx-qzplx*rho)/r1;
        VT1=VT2/r1;
        VR2=-gms*(qzminx+qzplx*rho)/r2;
        VT2=VT2/r2;
        
        V1Vec(1,i)=VR1;
        V1Vec(2,i)=VT1;
        V2Vec(1,i)=VR2;
        V2Vec(2,i)=VT2;
        alpha(i)=(1-x^2)*2*GM/s;
    end
end

function [x,xpl,N]=xlamb(m,q,qsqfm1,Tin,epsVal)
%%XLAMB  This is a subroutine in Gooding's algorithm for solving Lambert's
%        problem. This subroutine finds the appropriate value(s) of x.
%
%The function is taken from Appendix D of
%R. H. Gooding, "On the solution of Lambert's orbital boundary-value
%problem," Royal Aerospace Executive, Procurement Executive, Ministry
%of Defence, Farnborough, Hants, United Kingdom, Tech. Rep. 88027,
%Apr. 1988.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Empirical values from Appendix B of Goodwin's paper/ Appendix D of the
%report that are not all explicitly defined in the text.
c0=1.7;
c1=0.5;
c2=0.03;
c3=0.15;
c41=1;
c42=0.24;

%From Equation 19 in the Gooding report, but divied by 2*pi
thr2=atan2(qsqfm1,2*q)/pi;

%These empty matrices are only returned if it does not converge.
x=[];
xpl=[];

skipTo3=false;
if(m==0)
%In the no revolution case, a solution for x always exists (given numerical
%precision limits). First, we have to find the sign of x. Gooding does this
%by finding an initial value of T that corresponds to x=0
    N=1;
    T0=tlamb(m,q,qsqfm1,0,0);
    TDiff=Tin-T0;
    if(TDiff<=0)%x>=0
        %Equation 11 in Gooding's paper
        x=T0*TDiff/(-4*Tin);
        %(-4 is the value of dT for x=0)
    else%x<0
        %Equation 13 in Gooding's paper
        x=-TDiff/(TDiff+4);
        %Equation 16 in Gooding's paper
        W=x+c0*sqrt(2*(1-thr2));
        if(W<0)
            %Equation 17 in Gooding's paper
            x=x-(-W).^(1/16).*x+sqrt(TDiff/(TDiff+1.5*T0));
        end
        
        w=4/(4+TDiff);
        %Equations 17 and 18 in Gooding's report, combined
        x=x*(1+x*(c1*w-c2*x*sqrt(w)));

        %There should always be a solution in the single revolution case.
        %This should take care of finite precision problems.
        x=max(-1,x);
    end
else
%If one more more revolutions occurs, first get TMin as a basis for a starter
%From between equations 19 and 20 in Gooding's report.
    xM=1/(1.5*(m+0.5)*pi);  
    if(thr2<0.5)
        %Equation 20 in Gooding's report
        xM=nthroot(2*thr2,8)*xM;
    elseif(thr2>0.5)
        %Equation 21 in Gooding's report.
        xM=(2-nthroot(2-2*thr2,8))*xM;
    end
    %Starter for TMin
    %This is an application of Halley's method.
    maxIter=12;
    didConverge=false;
    for i=1:maxIter
        [TMin,dT,d2T,d3T]=tlamb(m,q,qsqfm1,xM,3);
        %If the second derivative is zero. Gooding compared to zero; I used
        %a small value to help minimize finite precision errors.
        if(abs(d2T)<epsVal)
            didConverge=true;
            break;
        end
        xMold=xM;
        xM=xM-dT*d2T/(d2T^2-dT*d3T/2);
        
        if(abs(xMold/xM-1)<=epsVal)
            didConverge=true;
            break;
        end
    end
    
    %Break off and exit if tMin not located - should never happen.
    if(didConverge==false)
        N=0;%Returning zero indicates no solution found.
        return
    end
    
    %Find the corresponding TM.
    TDiffM=Tin-TMin;
    
    %Gooding compared equality to zero, I use a small number epsVal to help
    %minimize finite precision errors.
    if(abs(TDiffM)<epsVal)
        x=xM;
        N=1;
        %Exit if unique solution already from X(TMin);
        return;
    elseif(TDiffM<0)
        N=0;
        return;%Exit if no solution exists with this m
    else
        N=3;
        if(d2T==0)
            d2T=6*m*pi;
        end
        %Equation 23 in Gooding's report
        x=sqrt(TDiffM/(d2T/2+TDiffM/(1-xM)^2));

        %This is the patched started mentioned in Section 5.3 of Gooding's
        %report, but not spelt out until Appendix D.
        W=xM+x;
        W=W*4/(4+TDiffM)+(1-W)^2;
        x=x*(1-(1+m+c41*(thr2-1/2))/(1+c3*m)*x*(c1*W+c2*x*sqrt(W)))+xM;
        d2T2=d2T/2;
        if(x>=1)%Otherwise, no finite solution with x> xM
            N=1;
            skipTo3=true;
        end%Otherwise, no finite solution with x>xM
    end
end

%We now have an initial estimate of x. The true value of x will
%be found iteratively using Halley's method. The iteration for the next
%estimate of x in Halley's method is
%x=x-2*f*dfdx/(2*dfdx^2-f*d2fdx2)
%where f is a function of x and dfdx and d2fdx2 are respectively the first
%and second derivatives of the function. Whereas Newton's method requires
%the first derivative, Halley's method requires the second derivative to
%find the zero of the function f. Thus, if one wants to find an extrema of
%a function, then f would be the derivative of the function with respect to
%x, so one would need up to the third derivative of the function. Halley's
%method is described at
%Weisstein, Eric W. "Halley's Method." From MathWorld--A Wolfram Web
%Resource. http://mathworld.wolfram.com/HalleysMethod.html
%
%In the problem at hand, we want to minimize the difference between T,
%which is the true (normalized) time of flight needed to go between the
%points r1Vec and r2Vec, and the value obtained by the function
%findTAndDerivs, which is the time implied by the value x. The iterations
%inside the loop in the if(~skipTo3) statement do that.

while(1)%This is point 5 in Gooding's Fortran code.
    if(~skipTo3)
        %Gooding suggest 3 iterations. However, I noticed that it often
        %takes more to converge, though 60 should be much more than enough.
        maxIter=60;
        for i=1:maxIter
            [T,dT,d2T]=tlamb(m,q,qsqfm1,x,2);
            T=Tin-T;
            %Gooding compared dT to zero. Here, I compare the change to the
            %smallest possible change in x.
            deltaVal=T*dT/(dT^2+T*d2T/2);
            if(abs(deltaVal)>=eps(x))
               x=x+deltaVal;
            else
                break;
            end
        end

        %Return if only one solution, normally when m=0
        if(N~=3)
            return;
        end
        N=2;
        xpl=x;
    end
    
    %Second multi-rev starter
    T0=tlamb(m,q,qsqfm1,0,0);
    TDiff0=T0-TMin;
    TDiff=Tin-T0;
    if(TDiff<=0)
        %Equation 25 in Gooding's report.
        x=xM-sqrt(TDiffM/(d2T2-TDiffM*(d2T2/TDiff0-1/xM^2)));
    else
        x=-TDiff/(TDiff+4);

        %This is comparable to Equation 16 in Gooding's paper
        W=x+c0*sqrt(2*(1-thr2));
        if(W<0)
            x=x-(-W)^(1/16)*(x+sqrt(TDiff/(TDiff+1.5*T0)));
        end
        W=4/(4+TDiff);
        x=x*(1+(1+m+c42*(thr2-0.5))/(1+c3*m)*x*(c1*W-c2*x*sqrt(W)));

        if(x<=-1)
            N=N-1;%No finite solution with x<xM
            if(N==1)
                x=xpl;
            end
        end
    end

    skipTo3=false;
end

end

function [T,dT,d2T,d3T]=tlamb(m,q,qsqfm1,x,N)
%%TLAMB A subroutine from Gooding's report. This computes the
%       time-of-flight along with necessary derivatives.
%
%The function is taken from Appendix C of 
%R. H. Gooding, "On the solution of Lambert's orbital boundary-value
%problem," Royal Aerospace Executive, Procurement Executive, Ministry
%of Defence, Farnborough, Hants, United Kingdom, Tech. Rep. 88027,
%Apr. 1988.
%and has been modified to also return the variable A, if it is computed.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

LM1=(N==-1);
L1=(N>=1);
L2=(N>=2);
L3=(N==3);

qsq=q^2;
xsq=x^2;

%Equation B-1, Gooding's report
u=(1-x)*(1+x);

%This will only return empty if
T=[];

if(~LM1)
   dT=0;
   d2T=0;
   d3T=0;
   A=[];
end

sw=0.4;

%If an explicit solution can be used rather than a series representation.
%The value of 0.4 is mention after Equation B-17 in the report.
if(LM1||m>0||x<0||abs(u) > sw)
    %Equation B-2, Gooding's report
    y=sqrt(abs(u));
    %Equation B-3, Gooding's report
    z=sqrt(qsqfm1+qsq*xsq);
    
    qx=q*x;
    if(qx<=0.0)
        %Equation B-4, Gooding's report (alpha)
        A=z-q*x;
        %Equation B-6, Gooding's report (beta)
        B=q*z-x;
    end
    
    if(qx<0&&LM1)
        AA=qsqfm1/A;
        BB=qsqfm1*(qsq*u-xsq)/B;
    end
    
    if(qx==0&&LM1||qx>0)
        %Equation B-5, Gooding's report
        AA=z+q*x;
        %Equation B-7, Gooding's report
        BB=q*z+x;
    end
    
    if(qx>0)
        %From Equation B-10, Gooding's report
        A = qsqfm1/AA;
        %From Equation B-11, Gooding's report
        B = qsqfm1*(qsq*u-xsq)/BB;
    end

    if(LM1)
       dT=B;
       d2T=BB;
       d3T=AA;
    else
        if (qx*u >= 0)
            %Equation B-9, Gooding's report
            g=x*z+q*u;
        else
            %Equation after B-11 Gooding's report
            g=(xsq-qsq*u)/(x*z-q*u);
        end
        
        %Equation B-8, Gooding's report
        f = A*y;

        if(x<=1)
            %Equation B-12a in Gooding's report for elliptical orbits.
            T = m*pi+atan2(f, g);
        else%If the path is hyperbolic.
            %This comes from the discussion after B-12b. However, rather 
            %than using a series representation, Matlab's implementation of
            %the logarithm function sufficies if we just make sure that the
            %argument is never less than zero due to numerical precision
            %limitations.
            T = log(max(0,f+g));
        end

        %Equation B-13, in Gooding's report.
        T=2*(T/y+B)/u;
        
        if(L1&&z~=0)
            %As noted in Appendix B of the report, z can only be zero if x
            %is zero and abs(q)=1. At such a point, the derivative is
            %discontinuous, so getting an indeterminant solution is
            %appropriate.
            qz=q/z;
            qz2=qz*qz;
            qz=qz*qz2;
            
            
            %Equation B-14, Gooding's report
            dT=(3*x*T-4*(A+qx*qsqfm1)/z)/u;

            if(L2)
                %Equation B-16, Gooding's report
                d2T=(3*T+5*x*dT+4*qz*qsqfm1)/u;
            end

            if(L3)
                %Equation B-17, Gooding's report
                d3T=(8*dT+7*x*d2T-12*qz*qz2*x*qsqfm1)/u;
            end
        end
    end
else
    %If a series solution to find T and its derivatives must be used.
    %The series is in Equation B-28 in the report. Derivatives are taken of
    %the series.

    %These are powers of u needed for the various derivatives (with respect
    %to u) of the sum.
    u0i=1;
    if(L1)
        u1i=1;
    end
    if(L2)
        u2i=1;
    end
    if(L3)
        u3i=1;
    end
    
    term=4;
    tq=q*qsqfm1;
    i=0;
    
    %Get the first term in the b series.
    if (q < 0.5)
        %For b0 from Equation B-23, Gooding's report.
        tqsum=1-q*qsq;
    else
        %Equation B-24, Gooding's report.
        tqsum=(1/(1+q)+q)*qsqfm1;
    end
    
    ttmold=term/3;
    T=ttmold*tqsum;
    told=inf;%So that it will enter the loop.
    while(i<N||T~=told)%loop until convergence.
        i=i+1;
        p=i;
        u0i=u0i*u;
        if(L1&&i>1)
            u1i=u1i*u;
        end
        if(L2&&i>2)
            u2i=u2i*u;
        end
        if(L3&&i>3)
            u3i=u3i*u;
        end
        
        %Get the next a value using Equation B-27
        term=term*(p-(1/2))/p;
        tq=tq*qsq;
        %Get the next b value using Equation B-25
        tqsum=tqsum+tq;
        told=T;
        %Equation B-26, Gooding's report.
        tterm=term/(2*p+3);
        %Equation B-22, Gooding's report.
        tqterm=tterm*tqsum;
        %Add the nth term of the recursion in B-27 to the running sum for
        %x^2*T.
        T=T-u0i*((1.5*p+0.25)*tqterm/(p^2-0.25)-ttmold*tq);
        ttmold=tterm;
        tqterm=tqterm*p;
        
        %The derivatives are computed using derivatives of the running sum
        %in B-21.
        if(L1)
            dT=dT+tqterm*u1i;
        end
        if(L2)
            d2T=d2T+tqterm*u2i*(p-1);
        end
        if(L3)
            d3T=d3T+tqterm*u3i*(p-1)*(p-2);
        end
    end
    
    if(L3)
        %Equation B-33, Gooding's report
        d3T=8*x*(1.5*d2T-xsq*d3T);
    end
    if(L2)
        %Equation B-32, Gooding's report
        d2T=2*(2*xsq*d2T-dT);
    end
    
    if(L1)
        %Equation B-31, Gooding's report
        dT=-2*x*dT;
    end
    
    T=T/xsq;%Divide out an x^2
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
