function J=JacobPolarLin2D(T,x,qTheta)
%%JACOBPOLARLIN2D Evaluate the Jacobian (matrix of partial derivatives with
%                 respect to the state elements) of a discrete-time state
%                 transition function for a constant-heading 2D dynamic
%                 model where the target state is given in term of
%                 position, the direction of the velocity measured
%                 counterclockwise from the x-axis, and speed. Optionally,
%                 a speed derivative component can be given to model
%                 possible linear target acceleration. This provides the
%                 Jacobian corresponding to the state transition function
%                 fTransPolarLin2D.
%
%INPUTS: x  The 4X1 or 5X1 target state vector in 2D space in the
%           order of [2D position;direction angle;speed; speed derivative]
%           where the heading angle is measured in radians counterclockwise
%           from the x axis and the speed derivative is akin to an
%           acceleration. If a speed derivative is not given (is assumed
%           zero), then x is 4X1.
%        T  The time-duration of the propagation interval in seconds.
%    qTheta If the transition function from a second-order weak
%           Itô-Taylor expansion is desired, then qTheta must be provided.
%           qTheta is the power spectral density of the process noise
%           corrupting the angular component of the state, having units of
%           radians^2/s. If qTheta is omitted, then the transition function
%           obtained by holding the nonlinear components constant over the
%           prediction interval is provided.
%
%OUTPUTS: J The 4X4 or 5X5 Jacobian matrix of the function
%           fTransPolarLin2D, where the partial derivatives are ordered
%           [df/dx(1), df/dx(2),...,df/dx(xDim)]. That is, column i
%           consists of partial derivatives with respect to element i of
%           the x vector.
%
%More comments on the model are given in the function fTransPolarLin2D.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    type='Discretized';
else
    type='ItoTaylorO2Weak';
end

xDim=size(x,1);

theta=x(3,1);%The heading (counterclockwise from the x axis).
v=x(4,1);%The speed

cosT=cos(theta);
sinT=sin(theta);

switch(xDim)
    case 4%If we are given position, heading, and linear velocity.
        switch(type)
            case 'Discretized'
                J=[1,0,-T*v*sinT,T*cosT;
                    0,1,T*v*cosT,T*sinT;
                    0,0,1,0;
                    0,0,0,1];
            case 'ItoTaylorO2Weak'
                J=[1,0,-(T-(T^2/4)*qTheta)*v*sinT,  (T-(T^2/4)*qTheta)*cosT;
                    0,1,(T-(T^2/4)*qTheta)*v*cosT,   (T-(T^2/4)*qTheta)*sinT;
                    0,0,1,                           0;
                    0,0,0,                           1];
            otherwise
                error('Invalid type provided');
        end
        
    case 5%If we are given position, heading, linear velocity, and linear
          %acceleration.
        vDot=x(5,:);%The derivative of the speed.
        switch(type)
            case 'Discretized'
                J=[1,0,-(T*v+(T^2/2)*vDot)*sinT,    T*cosT, (T^2/2)*cosT;
                    0,1,(T*v+(T^2/2)*vDot)*cosT,     T*sinT, (T^2/2)*sinT;
                    0,0,1,                           0,      0;
                    0,0,0,                           1,      T;
                    0,0,0,                           0,      1];
            case 'ItoTaylorO2Weak'
                J=[1,0,-((T-(T^2/4)*qTheta)*v+(T^2/2)*vDot)*sinT,   (T-(T^2/4)*qTheta)*cosT,  (T^2/2)*cosT;
                    0,1,((T-(T^2/4)*qTheta)*v+(T^2/2)*vDot)*cosT,    (T-(T^2/4)*qTheta)*sinT,  (T^2/2)*sinT;
                    0,0,1                                            0,                          0;
                    0,0,0                                            1,                          T;
                    0,0,0                                            0,                          1];
            otherwise
                error('Invalid type provided.')
        end
    otherwise
        error('The state vector has an invalid length.')
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
