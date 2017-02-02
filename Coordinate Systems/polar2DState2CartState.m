function xCart=polar2DState2CartState(xPol,systemType)
%%POLAR2DSTATE2CARTSTATE Convert a 2D target state where the velocity had
%               been decomposed into a direction angle (heading) and speed
%               components into Cartesian components. Depending on the
%               system type chosen, the state can have components for a
%               linear acceleration and/ or a turn rate.
%
%INPUTS: xPol  A 4X1, 5X1 or 6X1 target state where the first four
%              components are [position;heading;speed] in 2D. The other
%              components depend on the value of systemType.
%   systemType A string constant specifying the desired type of input and
%              output. In all instances, the heading is measured in terms
%              of radians counterclockwise from the x-axis. Possible values
%              are:
%              'ConstVel'     The target state is [position;heading;speed]
%                             and xCart is [position;velocity]
%              'ConstAccel'   The target state is [position;heading;speed;
%                             speed derivative] and xCart is
%                             [position;velocity;acceleration]
%              'ConstTurn'    The target state is [position;heading;speed;
%                             turn rate] and xCart is
%                             [position;velocity;acceleration]
%              'TurnAndAccel' The target state is [position;heading;speed;
%                             turnrate; speed derivative] and xCart is
%                             [position;velocity;acceleration]
%
%%OUTPUTS: xCart The state converted into 2D Cartesian coordinates
%                consisting of position and velocity and, depending on
%                systemType, possibly acceleration components.
%
%The use of 2D states where the heading and speed have been separated is
%discussed in [1] and [2].
%
%The opposite of this function is Cart2DState2PolarState.
%
%REFERENCES:
%[1] M. Busch and S. Blackman, "Evaluation of IMM filtering for an air
%    defense system application," in Proceedings of SPIE: Signal and Data
%    Processing of Small Targets, vol. 2561, 9 Jul. 1995, pp. 435-447.
%[2] J. L. Gertz, "Multisensor surveillance for improved aircraft
%    tracking," The Lincoln Laboratory Journal, vol. 2, no. 3, pp. 381-396,
%    1989.
%
%July 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get position
x=xPol(1);
y=xPol(2);
%Get heading
theta=xPol(3);
%Get speed
v=xPol(4);

switch(systemType)
    case 'ConstVel'
        %Given position, heading, and speed.
        xCart=[x;
               y;
               v*cos(theta);
               v*sin(theta)];
    case 'ConstAccel'
        %Given position, heading, speed, and linear acceleration.
        
        vDot=xPol(5);%Linear acceleration
        xCart=[x;
               y;
               v*cos(theta);
               v*sin(theta);
               vDot*cos(theta);
               vDot*sin(theta)];
    case 'ConstTurn'
        %Given position, heading, speed, and turn rate.
        
        omega=xPol(5);%Turn rate
        xCart=[x;
               y;
               v*cos(theta);
               v*sin(theta);
               -v*omega*sin(theta);
               v*omega*cos(theta)];
    case 'TurnAndAccel'
        %Given position, heading, speed, turn rate, and linear
        %acceleration.
        
        omega=xPol(5);%Turn rate
        vDot=xPol(6);%Linear acceleration
        xCart=[x;
               y;
               v*cos(theta);
               v*sin(theta);
               vDot*cos(theta)-v*omega*sin(theta);
               vDot*sin(theta)+v*omega*cos(theta)];
    otherwise
        error('Invalid system type given.')
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
