function val= aCoordTurn2D(x,t,turnType,tauTurn,tauLinAccel)
%%ACOORDTURN2D The continuous-time drift function for a 2D coordinated turn
%              model with a Cartesian state. The turn rate can be specified
%              in terms of a turn rate in radians per second, or in terms
%              of a transversal acceleration. Additionally, a linear
%              acceleration can be given. The turn rate and linear
%              acceleration can optionally have time constants associated
%              with them, like in the Singer model, modelling a tendancy to
%              eventually want to return to non-accelerating, straight-line
%              motion.
%
%INPUTS:  x  The target state for 2D motion. If there is no linear
%            acceleration (acceleration along the direction of motion),
%            then x can either be x=[x;y;xdot;ydot;omega], where omega is
%            the turn rate estimate in radians per second
%            counterclockwise from the x-axis or x=[x;y;xdot;ydot;at] where
%            at is the the transversal acceleration, which is orthogonal to
%            the velocity and is defined such that positive values of at
%            map to positive values of omega. If there is a linear
%            acceleration, then the target state is either
%            x=[x;y;xdot;ydot;omega;al] where omega is the turn rate and al
%            is the linear acceleration or the target state is
%            x=[x;y;xdot;ydot;at;al] if the turn is expressed in terms of a
%            transversal acceleration. The dimensionality of the state is
%            used to determine whether a linear acceleration component is
%            present. The linear acceleration component changes the speed.
%            That means that it acts in the direction of the velocity
%            vector.
%          t An unused time component so that aCoordTurn2D can be used with
%            Runge-Kutta methods that expect the function to take two
%            parameters.
%  turnType  A string specifying whether the turn is given in terms of a
%            turn rate in radians per second or a transversal acceleration
%            in m/s^2. Possible values are
%            'TurnRate'   The turn is specified in terms of a turn rate
%                         (The default if this parameter is omitted).
%            'TransAccel' The turn is specified in terms of a transversal
%                         acceleration.
%   tauTurn  The correlation time constant for the turn rate in seconds.
%            tau must be positive but does not have to be finite. If this
%            parameter is omitted, then tauTurn is set to infinity.
%tauLinAccel The correlation time constant for the linear acceleration (if
%            present) in seconds. This parameter is not needed if there is
%            no linear acceleration. If a linear acceleration is present
%            and this parameter is omitted, then tauLinAccel is set to
%            infinity.
%
%OUTPUTS: val The time-derivative of the state under the 2D coordinated
%             turn model.
%
%The basic 2D coordinated turn model in Cartesian coordinates is described
%in Section VA of [1]. When the turn rate is something that must be
%estimated, it is assumed that the continuous-time turn rate model is
%omegaDot=-(1/tauTurn)*Omega+noise
%Note that the ordering of the state elements assumed by this function
%differs from the ordering of the state elements assumed in [1].
%
%The 2D coordinates turn model in Cartesian coordinates is also described
%in Chapter 4.2.3 of [2].
%
%The concept of using the transversal acceleration instead of the turn rate
%is not discussed in either of those references. It is, however, mentioned
%in [3], though no differential equations are given and a more detailed
%reference cited therein is a hard-to-get dissertation in French. The use
%of transversal acceleration is discussed in more detail in [4], though
%expressions are given when considering the 2D velocity are broken
%into components of heading and speed rather than in Cartesian space. The
%generalization to Cartesian space is not difficult and is done here.
%
%It can be shown that the magnitude of the transverse acceleration is equal
%to the speed times the turn rate. The model utilizing the transverse
%acceleration comes naturally from there. The relationship to the time
%constant in the transverse acceleration model comes from how it related to
%the derivative of the turn rate in [1].
%
%The optional linear acceleration provides a derivative of the speed of the
%target. The time constant operates in the same manner as the time constant
%for the turn rate. That is,
%alDot=-(1/tauLinAccel)*al+noise
%This is similar to how (total) acceleration decays in the Singer dynamic
%model, which was described in [5].
%
%The formulation in terms of transverse acceleration contains a singularity
%in the computation of the implied turn rate if the speed is zero. In such
%an instance, the turn rate is just set to zero. Similarly, singularities
%exist in applying the linear acceleration when the velocity of the target
%is zero. Values are zero are also substituted in such instances.
%
%The corresponding diffusion matrix is given by the function DCoordTurn2D.
%The corresponding discrete-time functions are FCoordTurn2D and
%QCoordTurn. However, note that the discrete-time functions are
%direct-discrete models and not discretizations of the continuous-time
%models as the propagated PDF does not remain Gaussian over time.
%
%REFERENCES:
%[1] X. R. Li and V. P. Jilkov, "Survey of maneuvering target tracking.
%    Part I: Dynamic models," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 39, no. 4, pp. 1333-1364, Oct. 2003.
%[2] S. Blackman and R. Popoli, Design and Analysis of Modern Tracking
%    Systems. Norwood, MA: Artech House, 1999.
%[3] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[4] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design
%    of a multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%[5] R. A. Singer,"Estimating optimal tracking filter performance for
%    manned maneuvering targets," IEEE Transactions on Aerospace and
%    Electronic Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Extract the components of velocity.
xDot=x(3);
yDot=x(4);

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
            case 5%There is no linear acceleration.
            %This uses the relationship between turn rate, velocity and
            %transverse acceleration implied by Equation 2.A25 in
            %Blom's paper and puts it into Equations 61 and 67 in Li's
            %paper. The time-constant, which is not in Blom's paper,
            %was added using the same logic as Li's paper, but
            %substituting the appropriate terms in the expression for
            %the derivative of the turn rate.
                val=[xDot;%Position derivative
                     yDot;%Position derivative
                     -omega*yDot;%Velocity derivative
                     omega*xDot;%Velocity derivative
                     -(1/tauTurn)*at];%Transverse acceleration derivative
            case 6%There is a linear acceleration component.
                al=x(6);
                %Get the linear acceleration vector. There will be NaNs if
                %the velocity is zero.
                linAccel=([xDot;yDot]/v)*al;

                %Deal with problems with a zero or nearly zero velocity.
                if(any(~isfinite(linAccel)))
                    linAccel(:)=0;
                end 
                
                val=[xDot;%Position derivative
                    yDot;%Position derivative
                    -omega*yDot+linAccel(1);%Velocity derivative
                    omega*xDot+linAccel(2);%Velocity derivative
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
                val=[xDot;%Position dervative
                     yDot;%Position derivative
                     -omega*yDot;%Velocity derivative
                     omega*xDot;%Velocity derivative
                     -(1/tauTurn)*omega];%Turn rate derivative
                
            case 6%There is a linear acceleration component.
                al=x(6);
                %Get the linear acceleration vector. There will be NaNs if
                %the velocity is zero.
                linAccel=([xDot;yDot]/v)*al;

                %Deal with problems with a zero or nearly zero velocity.
                if(any(~isfinite(linAccel)))
                    linAccel(:)=0;
                end
                
                %From Equations 61 and 67 in Li's paper with linear
                %acceleration added.
                val=[xDot;%Position dervative
                     yDot;%Position derivative
                     -omega*yDot+linAccel(1);%Velocity derivative
                     omega*xDot+linAccel(2);%Velocity derivative
                     -(1/tauTurn)*omega;%Turn rate derivative
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
