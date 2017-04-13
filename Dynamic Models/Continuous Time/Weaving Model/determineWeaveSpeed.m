function [speed,alphaVal,A]=determineWeaveSpeed(s,tEnd,Nw,betaVal)
%%DETERMINEWEAVESPEED This function returns the speed that a target needs
%       to travel when pursuing a weaving trajectory using the continuous
%       time dynamic model aWeave of duration tEnd seconds to travel a
%       distance of s meters along its initial heading while performing Nw
%       weaves. The trajectory's lack of "straightness" is parameterized by
%       beta. See the comments to the function aWeave for an example of its
%       use.
%
%INPUTS: s The desired distance traveled in meters.
%     tEnd The duration of flight in seconds.
%       Nw The integer number of sinusoidal periods in the maneuver
%          sequence before reaching the destination.
%  betaVal A parameter describing how wavy the trajectory is; 0<betaVal<=1.
%
%OUTPUTS: speed The speed in meters per second that the target must travel
%               along its trajectory.
%             A A term determining the magnitude of the deflection during
%               turns.
%         alpha The positive scalar weave period in radians per second.
%
%The weaving motion model and this function are derived in [1]. This
%function implements Equation 103 in [1]. For convenience, the parameters A
%and alphaVal are also returned to be passed to the function aWeave.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    intFunc=@(t)intArg(t,tEnd);
    intVal=integral(intFunc,0,tEnd);

    speed=s/intVal;
    
    alphaVal=2*pi*Nw/tEnd;
    A=(pi^2*Nw/tEnd)*betaVal;
 
    function val=intArg(t,tEndCur)
        val=cos((pi/2)*betaVal.*sin(2*pi*Nw*t./tEndCur));
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
