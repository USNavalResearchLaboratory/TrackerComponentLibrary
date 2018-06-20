function val=aPolarCoordTurn2D(x,t,turnType,tauTurn,tauLinAccel)
%APOLARCOORDTURN2D The continuous-time drift function for a 2D coordinated
%           turn model with the velocity expressed in terms of a heading
%           angle and a speed rather than in terms of a Cartesian velocity
%           vector. The turn rate can be specified in terms of a turn rate
%           in radians per second, or in terms of a transversal
%           acceleration. Additionally, a linear acceleration can be given.
%           The turn rate and linear acceleration can optionally have time
%           constants associated with them, like in the Singer model,
%           modelling a tendancy to eventually want to return to
%           non-accelerating, straight-line motion.
%
%INPUTS: x   The target state for 2D motion where the velocity is given in
%            terms of heading and speed components. If there is no linear
%            acceleration (acceleration along the direction of motion),
%            then x can either be x=[x;y;h;v;omega], where h is the heading
%            in terms of radians counterclockwise from the x-axis, v is the
%            speed, and omega is the turn rate (the derivative of h with
%            respect to time) or  x=[x;y;h;v;at] where at is the the
%            transversal acceleration, which is orthogonal to
%            the velocity and is defined such that positive values of at
%            map to positive values of omega. If there is a linear
%            acceleration, then the target state is either
%            x=[x;y;h;v;omega;al] where omega is the turn rate and al
%            is the linear acceleration or the target state is
%            x=[x;y;h;v;at;al] if the turn is expressed in terms of a
%            transversal acceleration. The dimensionality of the state is
%            used to determine whether a linear acceleration component is
%            present. The linear acceleration component changes the speed.
%            That means that it is the derivative of the speed.
%        t  An unused time component so that aPolar2DCV can be used with
%           Runge-Kutta methods that expect the function to take two
%           parameters.
% turnType  A string specifying whether the turn is given in terms of a
%           turn rate in radians per second or a transversal acceleration
%           in m/s^2. Possible values are
%           'TurnRate'   The turn is specified in terms of a turn rate
%                       (The default if this parameter is omitted).
%           'TransAccel' The turn is specified in terms of a transversal
%                        acceleration.
%   tauTurn  The correlation time constant for the turn rate in seconds.
%            tau must be positive but does not have to be finite. If this
%            parameter is omitted, then tauTurn is set to infinity.
%tauLinAccel The correlation time constant for the linear acceleration (if
%            present) in seconds. This parameter is not needed if there is
%            no linear acceleration. If a linear acceleration is present
%            and this parameter is omitted, then tauLinAccel is set to
%            infinity.
%       
%OUTPUTS: val The time-derivative of the state under the polar coordinated
%             turn motion model in 2D. val has the same dimensionality as
%             x.
%
%The basic 2D polar coordinated turn model is described in [1]. It is also
%mentioned in [2], though no differential equations are given and a more
%detailed reference cited therein is a hard-to-get dissertation in French.
%The use of transversal acceleration is discussed in more detail in [3].
%Chapter 4.2.3 of [4] also provides an overview of turning models using
%polar coordinates.
%
%More information on turning modeling using trasverse and linear
%acceleration can be found in the comments to the function aCoordTurn2D.
%
%The corresponding diffusion matrix is given by the function
%DPolarCoordTurn2D. The corresponding discrete-time functions are
%FPolarCoordTurn2D and QPolarCoordTurn2D. However, note that the
%discrete-time functions are direct-discrete models and not discretizations
%of the continuous-time models as the propagated PDF does not remain
%Gaussian over time.
%
%REFERENCES:
%[1] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%[2] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[3] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design
%    of a multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%[4] S. Blackman and R. Popoli, Design and Analysis of Modern Tracking
%    Systems. Norwood, MA: Artech House, 1999.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Extract the heading and speed (polar velocity components).
theta=x(3);
v=x(4);

%The time constant for the linear acceleration.
if(nargin<5)
    tauLinAccel=Inf;
end

%The time constant for the turn.
if(nargin<4)
   tauTurn=Inf;
end

if(nargin<3)
    turnType='TurnRate';
end
    
switch(turnType)
    case 'TransAccel'%The turn is expressed in terms of a transverse
                     %acceleration.     
        at=x(5);%The transverse acceleration.
        v=sqrt(xDot^2+yDot^2);%The speed

        omega=at/v;
        if(~isfinite(omega))%Deal with zero velocity.
            omega=0;
        end
                     
        switch(length(x))
            case 5%There is no linear acceleration
            %This uses the relationship between turn rate, velocity and
            %transverse acceleration implied by Equation 2.A25 in
            %Blom's paper and puts it into Equations 61 and 67 in Li's
            %paper. The time-constant, which is not in Blom's paper,
            %was added using the same logic as Li's paper, but
            %substituting the appropriate terms in the expression for
            %the derivative of the turn rate.
                val=[v*cos(theta);%Position derivative
                     v*sin(theta);%Position derivative
                     omega;%Heading derivative
                     0;%Speed derivative
                     -(1/tauTurn)*at];%Transverse acceleration derivative
            case 6%There is a linear acceleration component.
                al=x(6);%The linear acceleration (vDot)
                
                val=[v*cos(theta);%Position derivative
                    v*sin(theta);%Position derivative
                    omega;%Heading derivative
                    al;%Speed derivative
                    -(1/tauTurn)*at;%Transverse acceleration derivative
                    -(1/tauLinAccel)*al];%Linear acceleration derivative
            otherwise
                error('The dimensionality of the state is neither 5 nor 6.');
        end
    case 'TurnRate'%The turn is expressed in terms of a turn rate.
        omega=x(5);%The turn rate
        
        switch(length(x))
            case 5%There is no linear acceleration
                %From Equations 61 and 67 in Li's paper.
                val=[v*cos(theta);%Position dervative
                     v*sin(theta);%Position derivative
                     omega;%Heading derivative
                     0;%Speed derivative
                     -(1/tauTurn)*omega];%Turn rate derivative
                
            case 6%There is a linear acceleration component.
                al=x(6);%The linear acceleration (vDot)
                
                %From Equations 61 and 67 in Li's paper with linear
                %acceleration added.
                val=[v*cos(theta);%Position dervative
                     v*sin(theta);%Position derivative
                     omega;%Heading derivative
                     al;%Speed derivative
                     -(1/tau)*omega;%Turn rate derivative
                     -(1/tauLinAccel)*al];%Linear acceleration derivative
            otherwise
                error('The dimensionality of the state is neither 5 nor 6.');
        end
    otherwise
        error('Unknown turn type specified.');
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

