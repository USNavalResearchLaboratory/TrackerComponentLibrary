function F=FPolarLin2D(T,x,qTheta)
%%FPOLARLIN2D Obtain the discrete-time state transition matrix for
%             a constant-heading 2D dynamic model where the target state is
%             given in terms of position, the direction of the velocity
%             measured counterclockwise from the x-axis and speed.
%             Optionally, a speed derivative component can be given to
%             model possible linear target acceleration. The matrix can be
%             obtained from a linearization or from a second-order weak
%             Itô-Taylor expansion.
%
%INPUTS: T  The time-duration of the propagation interval in seconds.
%        x  The 4X1 or 5X1 target state vector in 2D space in the order of
%           [2D position;direction angle;speed; speed derivative] where the
%           heading angle is measured in radians counterclockwise from the
%           x axis and the speed derivative is akin to an acceleration. If
%           a speed derivative is not given (is assumed zero), then x is
%           4X1.
%    qTheta If the state transition matrix from a second-order weak
%           Itô-Taylor expansion is desired, then qTheta must be provided.
%           qTheta is the power spectral density of the process noise
%           corrupting the angular component of the state, having units of
%           radians^2/s. If qTheta is omitted, then the discrete-time state
%           transition matrix obtained by holding the nonlinear components
%           constant is provided.
%
%OUTPUTS: F The state transition matrix under a discrete-time linear motion
%           model where the target state is given in terms of position,
%           heading, speed, and optionally, linear acceleration.
%
%The state transition matrix without acceleration is the same as the matrix
%for the turning model as taken from the unnumered transition equation
%after equation 1b in [1], which comes from the initial derivation in [2],
%where a different definition of heading direction is used. The expression
%comes from the limit of a turning model as the turn rates goes to zero.
%
%However, the discretized model can be derived from first principles based
%on the principles in Chapter 6 of [3], where a matrix exponential is used.
%
%The solution using a second-order weak Itô-Taylor expansion is taken from
%[4].
%
%The associated state prediction covariance matrix can be obtained using
%QPolarLin2D. The associated nonlinear transition function, which is useful
%if one wishes to use the Cubature Kalman filter for the prediction step,
%is fTransPolarLin2D.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[2] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%[4] D. Laneuville, "Polar versus Cartesian velocity models for maneuvering
%    target tracking with IMM," in Proceedings of the IEEE Aerospace
%    Conference, Big Sky, MT, 2-9 Mar. 2013.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    type='Discretized';
else
    type='ItoTaylorO2Weak';
end

xDim=size(x,1);

theta=x(3);%The heading (counterclockwise from the x axis).
sinHead=sin(theta);
cosHead=cos(theta);

switch(xDim)
    case 4%If we are given position, heading, and linear velocity.
        switch(type)
            case 'Discretized'
                F=[1,0,0, cosHead*T;
                   0,1,0, sinHead*T;
                   0,0,1,         0;
                   0,0,0,         1];
            case 'ItoTaylorO2Weak'
                F=[1,0,0,T*cosHead-(T^2/4)*qTheta*cosHead;
                   0,1,0,T*sinHead-(T^2/4)*qTheta*sinHead;
                   0,0,1,0;
                   0,0,0,1];
            otherwise
                error('Invalid type provided');
        end

    case 5%If we are given position, heading, linear velocity, and linear
           %acceleration.
        switch(type)
            case 'Discretized'
                F=[1,0,0,cosHead*T, (T^2/2)*cosHead;
                   0,1,0,sinHead*T, (T^2/2)*sinHead;
                   0,0,1,0,         0;
                   0,0,0,1,         T;
                   0,0,0,0,         1];
            case 'ItoTaylorO2Weak'
                F=[1,0,0,T*cosHead-(T^2/4)*qTheta*cosHead,(T^2/2)*cosHead;
                   0,1,0,T*sinHead-(T^2/4)*qTheta*sinHead,(T^2/2)*sinHead;
                   0,0,1,0,0;
                   0,0,0,1,T;
                   0,0,0,0,1];
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
