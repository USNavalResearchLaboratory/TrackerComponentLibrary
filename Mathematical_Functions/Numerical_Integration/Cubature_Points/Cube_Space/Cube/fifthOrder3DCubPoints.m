function [xi,w]=fifthOrder3DCubPoints()
%%FIFTHORDER3DCUBPOINTS Generate fifth-order cubature points
%               for integration over a 3-dimensional cube with vertices of 
%               (1,-1,-1), (-1,-1,-1), (-1,1,-1), (1,1,-1), (1,-1,1),
%               (-1,-1,1), (-1,1,1), and (1,1,1).
%
%INPUTS: None
%
%OUTPUTS: xi This is a 3XnumCubPoints set of points for the standard
%            cube.
%          w A 1XnumCubPoints set of cubature weights. This sums to the
%            volume of the standard cube (8).
%
%This function implements the points given in [1] (14 points).
%
%EXAMPLE:
%We compare a 4th-order moment computed using these cubature points
%to one computed using monomialIntCube (a 5th order moment would have
%just been 0). The results are the same within typical finite precision
%limits.
% [xi,w]=fifthOrder3DCubPoints();
% alpha=[2;0;2];
% theMoment=findMomentFromSamp(alpha,xi,w);
% intVal=monomialIntCube(alpha);
% RelErr=(theMoment-intVal)/intVal
%
%REFERENCES:
%[1] F. D. Witherden and P. E. Vincent, "On the identification of symmetric
%    quadrature rules for finite element methods," Computer and Mathematics
%    with Applications, vol. 69, no. 10, pp. 1232-1241, May 2015.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

M=[-0.79582242575422146326454882047613584616,                                          0,                                          0,   0.88642659279778393351800554016620498615;
                                        0,                                          0,   0.79582242575422146326454882047613584616,   0.88642659279778393351800554016620498615;
                                        0,   0.79582242575422146326454882047613584616,                                          0,   0.88642659279778393351800554016620498615;
                                        0,                                          0,  -0.79582242575422146326454882047613584616,   0.88642659279778393351800554016620498615;
 0.79582242575422146326454882047613584616,                                          0,                                          0,   0.88642659279778393351800554016620498615;
                                        0,  -0.79582242575422146326454882047613584616,                                          0,   0.88642659279778393351800554016620498615;
 0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
-0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
-0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
-0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
-0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
 0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
 0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039;
 0.75878691063932814626903427811226742764,  -0.75878691063932814626903427811226742764,   0.75878691063932814626903427811226742764,   0.33518005540166204986149584487534626039];

w=M(:,4);
xi=M(:,1:3)';

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
