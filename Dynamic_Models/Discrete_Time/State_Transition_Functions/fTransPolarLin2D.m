function xPred=fTransPolarLin2D(T,x,qTheta)
%%FTRANSPOLARLIN2D Evaluate a discrete-time state transition function for a
%                  constant-heading 2D dynamic model where the target state
%                  is given in term of position, the direction of the
%                  velocity measured counterclockwise from the x-axis, and
%                  speed. Optionally, a speed derivative component can be
%                  given to model possible linear target acceleration. This
%                  function can be useful for state prediction using the
%                  Cubature Kalman filter.
%
%INPUTS: T  The time-duration of the propagation interval in seconds.
%        x  The 4XN or 5XN set of N target state vectors in 2D space in the
%           order of [2D position;direction angle;speed; speed derivative]
%           where the heading angle is measured in radians counterclockwise
%           from the x axis and the speed derivative is akin to an
%           acceleration. If a speed derivative is not given (is assumed
%           zero), then x is 4XN.
%   qTheta  If the transition function from a second-order weak
%           Itô-Taylor expansion is desired, then qTheta must be provided.
%           qTheta is the power spectral density of the process noise
%           corrupting the angular component of the state, having units of
%           radians^2/s. If qTheta is omitted, then the transition function
%           obtained by holding the nonlinear components constant over the
%           prediction interval is provided.
%
%OUTPUTS: xPred The 4XN or 5XN set of N state vectors after being predicted
%               forward a duration of T in time using the method depending
%               on whether qTheta is provided.
%
%The solution using a second-order weak Itô-Taylor expansion is taken from
%[1].
%
%The discretized model is inspired from a discussion in [2], where a
%different definition of heading direction is used. The expression comes
%from the limit of a turning model as the turn rates goes to zero.
%
%However, the discretized model can be derived from first principles based
%on the principles in Chapter 6 of [3], where a matrix exponential is used.
%
%The associated state prediction covariance matrix can be obtained using
%QPolarLin2D. The associated state transition matrix, which can be useful
%if one wants to directly use a Kalman filter, can be obtained from
%FPolarLin2D.
%
%REFERENCES:
%[1] D. Laneuville, "Polar versus Cartesian velocity models for maneuvering
%    target tracking with IMM," in Proceedings of the IEEE Aerospace
%    Conference, Big Sky, MT, 2-9 Mar. 2013.
%[2] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%[3] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    type='Discretized';
else
    type='ItoTaylorO2Weak';
end

numX=size(x,2);
xDim=size(x,1);

xPred=zeros(xDim,numX);
theta=x(3,:);%The heading (counterclockwise from the x axis).
v=x(4,:);%The speed

cosT=cos(theta);
sinT=sin(theta);

switch(xDim)
    case 4%If we are given position, heading, and linear velocity.
        switch(type)
            case 'Discretized'               
                xPred(1,:)=x(1,:)+T*v.*cosT;
                xPred(2,:)=x(2,:)+T*v.*sinT;
                xPred(3,:)=x(3,:);
                xPred(4,:)=x(4,:);
            case 'ItoTaylorO2Weak'
                %This formulation is directly in Laneuville's paper, except
                %he defines the heading angle differently.
                xPred(1,:)=x(1,:)+T*v.*cosT-(T^2/4)*qTheta*v.*cosT;
                xPred(2,:)=x(2,:)+T*v.*sinT-(T^2/4)*qTheta*v.*sinT;
                xPred(3,:)=x(3,:);
                xPred(4,:)=x(4,:);
            otherwise
                error('Invalid type provided');
        end
    case 5%If we are given position, heading, linear velocity, and linear
           %acceleration.
           vDot=x(5,:);%The deriavtive of the speed.
        switch(type)
            case 'Discretized'
                xPred(1,:)=x(1,:)+T*v.*cosT+(T^2/2)*vDot.*cosT;
                xPred(2,:)=x(2,:)+T*v.*sinT+(T^2/2)*vDot.*sinT;
                xPred(3,:)=x(3,:);
                xPred(4,:)=x(4,:)+T*vDot;
                xPred(5,:)=x(5,:);
            case 'ItoTaylorO2Weak'
                xPred(1,:)=x(1,:)+T*v.*cosT-(T^2/4)*v.*qTheta.*cosT+(T^2/2)*vDot.*cosT;
                xPred(2,:)=x(2,:)+T*v.*sinT-(T^2/4)*v.*qTheta.*sinT+(T^2/2)*vDot.*sinT;
                xPred(3,:)=x(3,:);
                xPred(4,:)=x(4,:)+T*vDot;
                xPred(5,:)=x(5,:);
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
