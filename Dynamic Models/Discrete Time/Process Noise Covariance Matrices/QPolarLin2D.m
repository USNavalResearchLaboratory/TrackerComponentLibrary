function Q=QPolarLin2D(T,x,qTheta,param4,type)
%QPOLARLIN2D  Obtain the discrete-time state process noise covariance
%             matrix for a constant-heading 2D dynamic model where the
%             target state is given in terms of position, the direction of
%             the velocity measured counterclockwise from the x-axis and
%             speed. Optionally, a speed derivative component can be given
%             to model possible linear target acceleration. The matrix can
%             be obtained from a linearization or from a second-order weak
%             Itô-Taylor expansion.
%
%INPUTS: T The time-duration of the propagation interval in seconds.
%        x The 4X1 or 5X1 target state vector in 2D space in the order of
%          [2D position;direction angle;speed; speed derivative] where the
%          heading angle is measured in radians counterclockwise from the
%          x axis and the speed derivative is akin to an acceleration. If
%          a speed derivative is not given (is assumed zero), then x is
%          4X1.
%   qTheta The power spectral density of the process noise corrupting the
%          angular component of the state, having units of radians^2/s.
%   param4 Depending on the dimensionality of x, this is either the power
%          spectral density of the process noise corrupting the speed
%          (when x is 4X1), having units of m^2/s^3, or this is the power
%          spectral density of the linear acceleration, having units of
%          m^2/s^3.
%     type An optional parameter specifying the type of linearization used
%          to obtain the process noise covariance matrix. This can be
%          either
%          'Discretized' (The default if omitted), treat the nonlinear
%                       terms in the drift matrix as constant and solve
%                       for the state transition equation. Note that Q is
%                       always singular when using the discretized model.
%          'ItoTaylorO2Weak' Use an order 2.0 weak Itô Taylor expansion.
%                       This provides better results than the simple
%                       linear expansion. Note that Q is singular when
%                       using the weak Ito-Taylor expansion with linear
%                       acceleration in the model. Q is (generally) not
%                       singular when a linear acceleration is not used.
%
%OUTPUT: Q The process noise covariance matrix under a linear dynamic
%          model, possibly with linear acceleration.
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
%The associated state transition matrix can be obtained using
%FPolarLin2D. The associated state transition function, which can be useful
%in some Cubature Kalman filter implementation, is fTransPolarLin2D.
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

if(nargin<5)
   type='Discretized'; 
end

xDim=size(x,1);
theta=x(3);%The heading (counterclockwise from the x axis).
v=x(4);%The speed.
sinT=sin(theta);
cosT=cos(theta);

switch(xDim)
    case 4%If we are given position, heading, and linear velocity.
        qv=param4;
        
        switch(type)
            case 'Discretized'
                a=(T^3/3)*qv*cosT^2;
                b=(T^3/3)*qv*sinT*cosT;
                c=0;
                d=(T^2/2)*qv*cosT;
                e=(T^3/3)*qv*sinT^2;
                f=0;
                g=(T^2/2)*qv*sinT;
                h=T*qTheta;
                i=0;
                j=T*qv;
                
                Q=[a,b,c,d;
                   b,e,f,g;
                   c,f,h,i;
                   d,g,i,j];
            case 'ItoTaylorO2Weak'
                a=(T^3/3)*(qv*cosT^2+qTheta*v^2*sinT^2);
                b=(T^3/3)*cosT*sinT*(qv-qTheta*v^2);
                
                c=-(T^2/2)*qTheta*v*sinT;
                d=(T^2/2)*qv*cosT;
                e=(T^3/3)*(qTheta*v^2*cosT^2+qv*sinT^2);
                f=(T^2/2)*qTheta*v*cosT;
                g=(T^2/2)*qv*sinT;
                h=T*qTheta;
                i=0;
                j=T*qv;

                Q=[a,b,c,d;
                   b,e,f,g;
                   c,f,h,i;
                   d,g,i,j];
            otherwise
                error('Invalid type provided');
        end

    case 5%If we are given position, heading, linear velocity, and linear
           %acceleration.
           qvDot=param4;
           vDot=x(5);
        switch(type)
            case 'Discretized'
                a=(T^5/20)*qvDot*cosT^2;
                b=(T^5/20)*qvDot*sinT*cosT;
                c=0;
                d=(T^4/8)*qvDot*cosT;
                e=(T^3/6)*qvDot*cosT;
                f=(T^5/20)*qvDot*sinT^2;
                g=0;
                h=(T^4/8)*qvDot*sinT;
                i=(T^3/6)*qvDot*sinT;
                j=T*qTheta;
                k=0;
                l=0;
                m=(T^3/3)*qvDot;
                n=(T^2/2)*qvDot;
                o=T*qvDot;
                
                Q=[a,b,c,d,e;
                   b,f,g,h,i;
                   c,g,j,k,l;
                   d,h,k,m,n;
                   e,i,l,n,o];
            case 'ItoTaylorO2Weak'
                a=(T^3/3)*((T^2/4)*qvDot*cosT^2+qTheta*(v+(T/2)*vDot)^2*sinT^2);
                b=(T^3/3)*((T^2/4)*qvDot-qTheta*(v+(T/2)*vDot)^2)*sinT*cosT;
                c=-(T^2/2)*qTheta*(v+(T/2)*vDot)*sinT;
                d=(T^4/6)*qvDot*cosT;
                e=(T^3/4)*qvDot*cosT;
                f=(T^3/3)*((T^2/4)*qvDot*sinT^2+qTheta*(v+(T/2)*vDot)^2*cosT^2);
                g=(T^2/2)*qTheta*(v+(T/2)*vDot)*cosT;
                h=(T^4/6)*qvDot*sinT;
                i=(T^3/4)*qvDot*sinT;
                j=qTheta*T;
                k=0;
                l=0;
                m=(T^3/3)*qvDot;
                n=(T^2/2)*qvDot;
                o=T*qvDot;

                Q=[a,b,c,d,e;
                   b,f,g,h,i;
                   c,g,j,k,l;
                   d,h,k,m,n;
                   e,i,l,n,o];
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
