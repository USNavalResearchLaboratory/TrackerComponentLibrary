function Q=QPolarCoordTurn2D(T,x,sigmaTurn2,param4,turnType,tauTurn,tauLinAccel)
%%QPOLARCOORDTURN2D Get the discrete-time process noise covariance matrix
%              for a two-dimensional coordinated turn model with a state
%              where the velocity is given in terms of a heading and a
%              speed. The turn rate can be specified in terms of a turn
%              rate in radians per second, or in terms of a transversal
%              acceleration. Additionally, a linear acceleration can be
%              given.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%        x The target state for 2D motion where the velocity is given in
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
%          turn rate and al is the linear acceleration or the target state
%          is x=[x;y;h;v;at;al] if the turn is expressed in terms of a
%          transversal acceleration. The dimensionality of the state is
%          used to determine whether a linear acceleration component is
%          present. The linear acceleration component changes the speed.
%          That means that it is the derivative of the speed.
% sigmaTurn2 If the turn is specified in terms of a turn rate in radians
%          per second, then this is the instantaneous variance of the turn
%          rate noise. If the turn is expressed in terms of a transverse
%          acceleration, then this is the instantaneous variance of the
%          transverse acceleration noise.
%   param4 If a linear acceleration is present this is the instantaneous
%          variance of the linear acceleration noise. Otherwise this is
%          the instantaneous variance of the velocity.
% turnType A string specifying whether the turn is given in terms of a
%          turn rate in radians per second or a transversal acceleration
%          in m/s^2. Possible values are
%          'TurnRate'   The turn is specified in terms of a turn rate
%                       (The default if this parameter is omitted).
%          'TransAccel' The turn is specified in terms of a transversal
%                       acceleration.
%  tauTurn The correlation time constant for the turn rate in seconds.
%          tau must be positive but does not have to be finite. If this
%          parameter is omitted, then tauTurn is set to infinity.
% tauLinAccel The correlation time constant for the linear acceleration (if
%          present) in seconds. This parameter is not needed if there is no
%          linear acceleration. If a linear acceleration is present and
%          this parameter is omitted, then tauLinAccel is set to Inf.
%
%OUTPUTS: Q The process noise covariance matrix under a direct discrete-
%           time coordinated turn model where the velocity is specified in
%           terms of a heading and speed. Note that Q is singular due to
%           the direct-discrete modeling of the process noise.
%
%Details of the overall discretized dynamic model are given in the comments
%to the function FPolarCoordTurn2D.
%
%When no linear acceleration component is present and the turn is expressed
%in terms of a turn rate, then the covariance matrix for the model is based
%on the unnumbered equation on the second page of [1], which is page 436 of
%the proceedings. The model essentially puts noise on the turn rate
%component and allows it to integrate into the heading component; it also
%puts a noise term on the speed. When a linear acceleration term is
%present, then one could extend the concept to putting a noise term on the
%linear acceleration and letting it integrate into the speed. In the
%instance where a transverse acceleration is given, one can approximate the
%conversion based on the fact that omega*v=at and scaling the covariance
%terms appropriately from the model used in Blackman.
%
%The paper [1] has a mysterious c_Q term, which he just set to 1. Upon
%inspection, the c_Q term appears to arise, because his process noise
%components are based off of a Singer (actually an integrated
%Ornstein-Uhlenbeck) model. Since the transition matrix
%function FPolarCoordTurn2D allows for the parameters tauTurn, and 
%tauLinAccel to specify the correlation-time of the noise for the turn
%component and for the linear acceleration, the parameters used by Blackman
%have been completely replaced with equivalent parameters in the context of
%the integrated Ornstein-Uhlenbeck model so that the correlation-time
%components could be used.
%
%Using the integrated Ornstein-Uhlenbeck model as a basis provides clues on
%how to set sigmaTurn2 and sigmaLin2.
%
%If qTurn is the value of the power spectral density of a continuous
%dynamic model driving the turn rate that has been tuned for particular
%parameters, using, for example, the function processNoiseSuggest
%with the 'PolyKalDirectDisc' model of order 1, then
%sigmaTurn2^2=qTurn/(2*tauTurn) is a good starting place to choose the
%noise for the model as per Section 8.2.2 of the book [2] when considering
%that the integrated Ornstein Uhlenbeck model has one level of integration
%less than the Singer model. Similar things can be done when the turn is
%represented in terms of a transverse acceleration and not a turn rate.
%
%Note that the discrete-time model in [1] does not integrate the errors
%into the position components. Thus, for large values of T, one would
%expect the discrete-time model to be particularly bad.
%
%The corresponding transition matrix is given by the function
%FPolarCoordTurn2D. The corresponding continuous-time functions are
%aPolarCoordTurn2DOmega, aPolarCoordTurn2DTrans, and DPolarCoordTurn2D.
%However, note that the discrete-time functions with unknown noise is a
%direct-discrete model and not a discretization of the continuous-time
%model.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[2] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7)
    tauLinAccel=Inf;
end

if(nargin<6)
    tauTurn=Inf;
end

if(nargin<5)
    turnType='TurnRate';
end


switch(turnType)
    case 'TransAccel'%The turn is expressed in terms of a transverse
                     %acceleration.
    %If the turn rate is expressed in terms of a transverse acceleration,
    %then the equivalent standard deviation for the instantaneous turn rate
    %must be obtained so that the Ornstein-Uhlenbeck model can be used for
    %the variance terms. The transformation is done assuming that the speed
    %is constant.
        v=x(4);
        sigmaOmega=sqrt(sigmaTurn2/v^2);
        if(~isfinite(sigmaOmega))
            sigmaOmega=0;
        end
        
        %Compute QTurn in terms of a turn parameterized by a turn rate.
        QTurn=QGaussMarkov(T,x([3,5],:),sigmaOmega^2,tauTurn,1);
        
        %Now, transform the appropriate elements of QTurn if a transverse
        %acceleration is used instead of a turn rate.
        QTurn(2,2)=QTurn(2,2)*v^2;
        QTurn(1,2)=QTurn(1,2)*v;
        QTurn(2,1)=QTurn(2,1)*v;
    case 'TurnRate'%The turn is expressed in terms of a turn rate.
        QTurn=QGaussMarkov(T,x([3,5],:),sigmaTurn2,tauTurn,1);
    otherwise
        error('Unknown turn type specified.');
end

switch(length(x))
    case 5%There is no linear acceleration.
%The covariance matrix suggested in Blackman's paper when given position,
%heading, speed, and turn rate.
        Q=[0,   0,  0,                  0,          0;
           0,   0,  0,                  0,          0;
           0,   0,  QTurn(1,1),         0,          QTurn(1,2);
           0,   0,  0,                  T*param4,  0;
           0,   0,  QTurn(2,1),         0,          QTurn(2,2)];
    case 6%There is a linear acceleration component.
        QAccel=QGaussMarkov(T,x([4,6],:),param4,tauLinAccel,1);

%An extension of the concept in Blackman's paper to the case where we are
%given heading, speed, turn rate, and linear acceleration.
        Q=[0,   0,  0,              0,               0,                  0;
           0,   0,  0,              0,               0,                  0;
           0,   0,  QTurn(1,1),     0,               QTurn(1,2),         0;
           0,   0,  0,              QAccel(1,1),     0,                  QAccel(1,2);
           0,   0,  QTurn(2,1),     0,               QTurn(2,2),         0;
           0,   0,  0,              QAccel(2,1),     0,                  QAccel(2,2)];
    otherwise
        error('The dimensionality of the state is neither 5 nor 6.');
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
