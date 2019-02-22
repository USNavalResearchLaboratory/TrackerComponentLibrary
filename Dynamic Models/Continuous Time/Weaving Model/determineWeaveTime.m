function [tEnd,alphaVal,A]=determineWeaveTime(s,speed,Nw,betaVal)
%%DETERMINEWEAVETIME This function returns the amount of time a target
%       needs to travel when pursuing a weaving trajectory using the
%       continuous time dynamic model aWeave to travel a distance of s at a
%       given speed while performing Nw weaves. The trajectory's lack of
%       "straightness" is parameterized by beta.  See the comments to the
%       function aWeave for an example of its use.
%
%INPUTS: s The desired distance traveled in meters.
%     speed The speed in meters per second of the target.
%       Nw The integer number of sinusoidal periods in the maneuver
%          sequence before reaching the destination.
%  betaVal A parameter describing how wavy the trajectory is; 0<betaVal<=1.
%
%OUTPUTS: tEnd The duration of flight in seconds.
%            A A term determining the magnitude of the deflection during
%              turns.
%        alpha The positive scalar weave period in radians per second.
%
%This function implements Equation 102 in [1]. For convenience, the
%parameters A and alphaVal are also returned to be passed to the function
%aWeave.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Solve Equation 102 in [1].
    tEnd=fminbnd(@minFunc,0,3*s/speed,optimset('TolX',1e-10));
    
    %Equation 98 in[1].
    alphaVal=2*pi*Nw/tEnd;
    A=pi*alphaVal/2*betaVal;

    function val=minFunc(tEndCur)
        numVals=length(tEndCur);
        intVals=zeros(numVals,1);
        
        for curVal=1:numVals
            intFunc=@(t)intArg(t,tEndCur(curVal));
            intVals(curVal)=integral(intFunc,0,tEndCur(curVal),'RelTol',1e-12);
        end
        
        val=(s-speed*intVals).^2;
    end

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
