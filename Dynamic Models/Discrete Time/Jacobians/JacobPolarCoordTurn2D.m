function J=JacobPolarCoordTurn2D(T,x,turnType,discPoint,tauTurn,tauLinAccel)
%%JACOBCOORDTURN2D Evaluate the Jacobian (matrix of partial derivatives
%              with respect to the state elements) of a discrete-time state
%              transition function for a two-dimensional coordinated turn
%              model where the velocity is specified in terms of a heading
%              and a speed. The turn rate can be specified in terms of a
%              turn rate in radians per second, or in terms of a
%              transversal acceleration. Additionally, a linear
%              acceleration can be given. The turn rate and linear
%              acceleration can optionally have time constants associated
%              with them, like in the Singer model, modelling a tendancy to
%              eventually want to return to non-accelerating, straight-line
%              motion. This Jacobian is associated with the continuous-time
%              model fTransPolarCoordTurn2D.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%        x The  target state for 2D motion where the velocity is given in
%          terms of heading and speed components. If there is no linear
%          acceleration (acceleration along the direction of motion),
%          then x can either be x=[x;y;h;v;omega], where h is the heading
%          in terms of radians counterclockwise from the x-axis, v is the
%          speed, and omega is the turn rate (the derivative of h with
%          respect to time) or  x=[x;y;h;v;at] where at is the
%          transversal acceleration, which is orthogonal to the velocity
%          and is defined such that positive values of at map to positive
%          values of omega. If there is a linear acceleration, then the
%          target state is either x=[x;y;h;v;omega;al] where omega is the
%          turn rate and al is the linear acceleration or the target
%          state is x=[x;y;h;v;at;al] if the turn is expressed in terms
%          of a transversal acceleration. The dimensionality of the state
%          is used to determine whether a linear acceleration component
%          is present. The linear acceleration component changes the
%          speed. That means that it is the derivative of the speed.
% turnType A string specifying whether the turn is given in terms of a
%          turn rate in radians per second or a transversal acceleration
%          in m/s^2. Possible values are
%          'TurnRate'   The turn is specified in terms of a turn rate
%                       (The default if this parameter is omitted).
%          'TransAccel' The turn is specified in terms of a transversal
%                       acceleration.
% discPoint This optional parameter specified what value of the turn rate
%           is used for the discretized state prediction. The thre
%           possible values were suggested in Li's paper cited
%           below. Possible values are
%           0 (The default if omitted) use omega=x(5), or the equivalent
%             value when specifying the turn rate using a transverse
%             acceleration, from the non-predicted target state for
%             building the state transition matrix. \omega_k
%           1 Use the average value of the predicted omega (or the average
%             value of the transverse acceleration) over the interval T
%             for building the state transition matrix. \bar{\omega}
%           2 Use the approximate average value of the predicted omega
%             over the interval T for building the state transition
%             matrix. \bar{\omega} This is the suggestion of using half
%             the prior and prediction, as was given in Li's paper. When
%             given a transverse acceleration instead of a turn rate, half
%             of the prior and predicted accelerations is used.
%           3 Use the forward-predicted omega (or transverse acceleration)
%             for building the state transition matrix. \omega_{k+1}
%   tauTurn The correlation time constant for the turn rate in seconds.
%           tau must be positive but does not have to be finite. If this
%           parameter is omitted, then tauTurn is set to infinity.
% tauLinAccel The correlation time constant for the linear acceleration (if
%            present) in seconds. This parameter is not needed if there is
%            no linear acceleration. If a linear acceleration is present
%            and this parameter is omitted, then tauLinAccel is set to
%            infinity.
%
%OUTPUTS: J The 5X5 or 6X6 Jacobian matrix of the function
%           fTransPolarCoordTurn2D, where the partial derivatives are
%           ordered [df/dx(1), df/dx(2),...,df/dx(xDim)]. That is, column i
%           consists of partial derivatives with respect to element i of
%           the x vector.
%
%More information on the discrete-time turning model is given in the
%comments to the function fTransPolarCoordTurn2D.
%
%April 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The time constant for the linear acceleration.
if(nargin<6||isempty(tauLinAccel))
    tauLinAccel=Inf;
end

%The time constant for the turn.
if(nargin<5||isempty(tauTurn))
   tauTurn=Inf;
end

if(nargin<4||isempty(discPoint))
    discPoint=0;
end

if(nargin<3||isempty(turnType))
    turnType='TurnRate';
end

xDim=size(x,1);

beta=exp(-T/tauTurn);
switch(discPoint)
    case 0%Use the prior value of at/omega for the linearization.
        turnVal=x(5);
    case 1
        %Use the average value of at.omega over the prediction
        %interval for the linearization.
        turnVal0=x(5);
        turnVal=(tauTurn-beta*tauTurn+T*turnVal0)/T;

        %Assume tau=Inf or T=0 was used, in which case the
        %asymptotic average value of at/omega should be used.
        turnVal(~isfinite(tauTurn)||~isfinite(turnVal))=turnVal0(~isfinite(tauTurn)||~isfinite(turnVal));
    case 2%Use the simple mean of at/omega for the linearization.
        turnVal=(beta+1)*x(5)/T;
    case 3%Use the predicted value of at/omega for the linearization.
        turnVal=beta*x(5);
    otherwise
        error('Invalid value entered for discPoint');
end

v=x(4);%The speed

switch(turnType)
    case 'TransAccel'%The turn is expressed in terms of a transverse
        %acceleration.
        omega=turnVal./v;
        omega(~isfinite(omega))=0;%Deal with zero velocity
    case 'TurnRate'%The turn is expressed in terms of a turn rate.
        omega=turnVal;
    otherwise
        error('Unknown turn type specified.');
end
%We now have the omega term for the turn.

theta=x(3);%The heading (counterclockwise from the x axis).

sinHead=sin(theta);
cosHead=cos(theta);

sinVal=sin(omega*T);
cosVal=cos(omega*T);

sinRat=sinVal./omega;
sinRat(~isfinite(sinRat))=T;%The limit as omega goes to zero.
cosRat=(1-cosVal)/omega;
cosRat(~isfinite(cosRat))=0;%The limit as omega goes to zero.

switch(turnType)
    case 'TransAccel'
        at=x(5);
        J=[1,0,-v*(sinHead*sinRat+cosHead*cosRat),2*(cosHead*sinRat-sinHead*cosRat)-T*(cosHead*cosVal-sinHead*sinVal),cosHead*(T*cosVal/omega-sinRat/omega)-sinHead*(T*sinRat-cosRat/omega);
           0,1,v*(cosHead*sinRat-sinHead*cosRat), 2*(sinHead*sinRat+cosHead*cosRat)-T*(sinHead*cosVal+cosHead*sinVal),sinHead*(T*cosVal/omega-sinRat/omega)+cosHead*(T*sinRat-cosRat/omega);
           0,0,1,                                 -T*at/v^2,                                                          T/v;
           0,0,0,                                 1,                                                                  0;
           0,0,0,                                 0,                                                                  beta];
    case 'TurnRate'
        J=[1,0,-v*(sinHead*sinRat+cosHead*cosRat),cosHead*sinRat-sinHead*cosRat,v*(cosHead*(T*cosVal/omega-sinRat/omega)-sinHead*(T*sinRat-cosRat/omega));
           0,1,v*(cosHead*sinRat-sinHead*cosRat), sinHead*sinRat+cosHead*cosRat,v*(sinHead*(T*cosVal/omega-sinRat/omega)+cosHead*(T*sinRat-cosRat/omega));
           0,0,1,                                 0,                            T;
           0,0,0,                                 1,                            0;
           0,0,0,                                 0,                            beta];
end
if xDim==6
    betaLinAccel=exp(-T/tauLinAccel);
    J(4,6)=T;
    J(6,6)=betaLinAccel;
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
