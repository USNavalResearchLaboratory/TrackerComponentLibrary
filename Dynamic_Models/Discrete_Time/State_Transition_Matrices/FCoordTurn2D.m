function F=FCoordTurn2D(T,x,turnType,discPoint,tauTurn,tauLinAccel)
%%FCOORDTURN2D The state transition matrix for a two-dimensional
%              coordinated turn model with a Cartesian state. The turn
%              rate can be specified in terms of a turn rate in radians per
%              second, or in terms of a transversal acceleration.
%              Additionally, a linear acceleration can be given. The turn
%              rate and linear acceleration can optionally have time
%              constants associated with them, like in the Singer model,
%              modelling a tendency to eventually want to return to
%              non-accelerating, straight-line motion.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%        x The target state for 2D motion. If there is no linear
%          acceleration (acceleration along the direction of motion), then
%          x can either be x=[x;y;xdot;ydot;omega], where omega is the turn
%          rate estimate in radians per second counterclockwise from the
%          x-axis or x=[x;y;xdot;ydot;at] where at is the transversal
%          acceleration, which is orthogonal to the velocity and is defined
%          such that positive values of at map to positive values of omega.
%          If there is a linear acceleration, then the target state is
%          either x=[x;y;xdot;ydot;omega;al] where omega is the turn rate
%          and al is the linear acceleration or the target state is
%          x=[x;y;xdot;ydot;at;al] if the turn is expressed in terms of a
%          transversal acceleration. The dimensionality of the state is
%          used to determine whether a linear acceleration component is
%          present. The linear acceleration component changes the speed.
%          That means that it acts in the direction of the velocity
%          vector.
% turnType A string specifying whether the turn is given in terms of a turn
%          rate in radians per second or a transversal acceleration in
%          m/s^2. Possible values are
%          'TurnRate'   The turn is specified in terms of a turn rate
%                       (The default if this parameter is omitted).
%          'TransAccel' The turn is specified in terms of a transversal
%                       acceleration.
% discPoint This optional parameter specified what value of the turn rate
%          is used for the discretized state prediction. The three
%          possible values were suggested in [1]. Possible values are:
%          0 (The default if omitted) use omega=x(5), or the equivalent
%            value when specifying the turn rate using a transverse
%            acceleration, from the non-predicted target state for
%            building the state transition matrix. \omega_k
%          1 Use the average value of the predicted omega (or the average
%            value of the transverse acceleration) over the interval T for
%            building the state transition matrix. \bar{\omega}
%          2 Use the approximate average value of the predicted omega over
%            the interval T for building the state transition matrix,
%            \bar{\omega}. This is the suggestion of using half the prior
%            and prediction, as was given in [2]. When given a
%            transverse acceleration instead of a turn rate, half of the
%            prior and predicted accelerations is used.
%          3 Use the forward-predicted omega (or transverse acceleration)
%            for building the state transition matrix. \omega_{k+1}
%  tauTurn The correlation time constant for the turn rate in seconds.
%          tau must be positive but does not have to be finite. If this
%          parameter is omitted, then tauTurn is set to infinity.
% tauLinAccel The correlation time constant for the linear acceleration (if
%          present) in seconds. This parameter is not needed if there is no
%          linear acceleration. If a linear acceleration is present and
%          this parameter is omitted, then tauLinAccel is set to Inf.
%
%OUTPUTS: F The 5X5 (no linear accelration) or 6X6 (with linear
%           acceleration) state transition matrix under a possibly
%           linearly accelerating coordinated turn model where the velocity
%           is a Cartesian vector.
%
%The basic 2D coordinated turn model is described in Section VA of [1]. It
%is assumed that the continuous-time evoluation of the turn rate omega is
%omegaDot=-(1/tauTurn)*omega+noise, which discretizes to
%omega[k+1]=exp(-T/tauTurn)*omega[k]+noise.
%Analogously, the optional linear acceleration discretizes to
%al[k+1]=exp(-T/tauLinAccel)*al[k]+noise
%The F matrix, not including the optional linear velocity term, is taken
%from equation 73 of [1] for an unknown turn rate. The ordering of the
%elements in the state has been changed from the paper.
%
%The concept of using the transversal acceleration instead of the turn rate
%is not discussed in either of those references. It is, however, mentioned
%in [2], though no differential equations are given and a more detailed
%reference cited therein is a hard-to-get dissertation in French. The use
%of transversal acceleration is discussed in more detail in [3], though
%expressions are given when considering the 2D velocity are broken into
%components of heading and speed rather than in Cartesian space. The
%generalization to Cartesian space is not difficult and is done here.
%
%The average value of the predicted value of omega over the interval that
%is used when the turn rate is unknown and discPoint=1 is the solution to
%the integral: (1/T)*int_0^T(omega_0+exp(-t/tauTurn))*dt. The noise term is
%not present as it is assumed to be zero mean. The approximation of [1] is
%to use 0.5*(omega_{k+1}+omega_k)  --that is, the average of the prior and
%predicted omega values over the discrete prediction interval. Analogous
%things can be done in terms of the transverse acceleration.
%
%The optional linear acceleration provides a derivative of the speed of the
%target. The time constant operates in the same manner as the time constant
%for the turn rate. That is,
%alDot=-(1/tauLinAccel)*al+noise
%This is similar to how (total) acceleration decays in the Singer dynamic
%model, which was described in [4]. The linear acceleration time T is just
%added to the velocity in the direction of motion of the target in this
%discretization. Care is taken to avoid NaNs when the velocity is close to
%zero. Care is also taken to avoid NaNs when omega is close to zero.
%
%This state transition matrix goes with the process noise covariance matrix
%given by QCoordTurn.The corresponding continuous-time drift functions for
%are aCoordTurn2DOmega and aCoordTurn2DTrans with the diffusion matrix
%DCoordTurn2D. Note that this model is a direct-discrete-time model and is
%not just a discretization of the continuous-time model.
%
%REFERENCES:
%[1] X. R. Li and V. P. Jilkov, "Survey of maneuvering target tracking.
%    Part I: Dynamic models," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 39, no. 4, pp. 1333-1364, Oct. 2003.
%[2] P. Vacher, I. Barret, and M. Gauvrit, "Design of a tracking algorithm
%    for an advanced ATC system," in Multitarget-Multisensor Tracking:
%    Applications and Advances, Y. Bar-Shalom, Ed. Norwood, MA: Artech
%    House, 1992, vol. II, ch. 1.
%[3] H. A. P. Blom, R. A. Hogendoorn, and B. A. van Doorn, "Design of a
%    multisensor tracking system for advanced air traffic control," in
%    Multitarget-Multisensor Tracking: Applications and Advances, Y. Bar-
%    Shalom, Ed. Norwood, MA: Artech House, 1992, vol. II, ch. 2.
%[4] R.A.Singer,"Estimating optimal tracking filter performance for manned
%    maneuvering targets," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. AES-6, no. 4, pp. 473-483, Jul. 1970.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The time constant for the linear acceleration.
if(nargin<6)
    tauLinAccel=Inf;
end

%The time constant for the turn.
if(nargin<5)
   tauTurn=Inf;
end

if(nargin<4)
    discPoint=0;
end

if(nargin<3)
    turnType='TurnRate';
end

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
        if(~isfinite(tauTurn)||~isfinite(turnVal))
            turnVal=turnVal0;
        end
    case 2%Use the simple mean of at/omega for the linearization.
        turnVal=(beta+1)*x(5)/T;
    case 3%Use the predicted value of at/omega for the linearization.
        turnVal=beta*x(5);
    otherwise
        error('Invalid value entered for discPoint');
end

%Extract the components of the velocity.
xDot=x(3);
yDot=x(4);
v=sqrt(xDot^2+yDot^2);%The speed

switch(turnType)
    case 'TransAccel'%The turn is expressed in terms of a transverse
                     %acceleration.
        omega=turnVal/v;
        if(~isfinite(omega))%Deal with zero velocity.
            omega=0;
        end
    case 'TurnRate'%The turn is expressed in terms of a turn rate.
        omega=turnVal;
    otherwise
        error('Unknown turn type specified.');
end
%We now have the omega term for the turn.

switch(length(x))
    case 5%There is no linear acceleration.
        F=zeros(5,5);%Allocate space.
    case 6%There is a linear acceleration component.
        F=zeros(6,6);%Allocate space.
        betaLinAccel=exp(-T/tauLinAccel);
        %The decay term for the linear acceleration.
        F(6,6)=betaLinAccel;
    otherwise
        error('The dimensionality of the state is neither 5 nor 6.');
end

%The decay term for the turn rate/ transverse acceleration.
F(5,5)=beta;

sinVal=sin(omega*T);
cosVal=cos(omega*T);

sinRat=sinVal/omega;
if(~isfinite(sinRat))
    sinRat=T;%The limit as omega goes to zero.
end
cosRat=(1-cosVal)/omega;
if(~isfinite(cosRat))
    cosRat=0;%The limit as omega goes to zero.
end

F(1:4,1:4)=[1,   0, sinRat,  -cosRat;%The x-component.
            0,   1, cosRat,   sinRat;%The y-component.
            0,   0, cosVal,  -sinVal;%The x-velocity component.
            0,   0, sinVal,   cosVal];%The y-velocity component.

%If there is a linear acceleration, then it gets added in as an extra
%component to the acceleration terms of the F matrix.
if(length(x)==6)
    al=x(6);%The linear acceleration.

    %The change in velocity due to the linear acceleration is the following
    %times each of the velocity components.
    deltaVParam=al*T/v;

    %Deal with problems with a zero or nearly zero velocity.
    if(any(~isfinite(deltaVParam)))
        deltaVParam=0;
    end
    
    %The following two equations are because deltaV=[xDot;yDot]*al*T/v.
    %Thus, the effects can be put into the transition matrix by allowing
    %the xDot and yDot components of the transition matrix to multiply the
    %appropriate other components.
    F(3,3)=F(3,3)+deltaVParam;
    F(4,4)=F(4,4)+deltaVParam;
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
