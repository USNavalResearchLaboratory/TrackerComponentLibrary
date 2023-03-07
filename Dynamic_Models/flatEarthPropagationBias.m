function vertBias=flatEarthPropagationBias(predDist,r)
%%FLATEARTHPROPAGATIONBIAS This is the vertical offset bias present when
%           predicting a constant-altitude target forward in a local
%           tangent plane versus compensating for the curvature of a
%           spherical Earth. This does not include horizontal biases. This
%           can be useful for determining when flat-Earth biases might
%           become significant.
%
%INPUTS: predDist The distance forward in the local tangent plant that one
%                 wishes to travel without compensating for the Earth's
%                 curvature.
%               r The radius of the spherical Earth. If this parameter is
%                 omitted or an empty matrix is passed, then the default of
%                 Constants.WGS84MeanRadius is used.
%
%OUTPUTS: vertBias The magnitude of the vertical elevation bias. 

%This implements Equation 58 in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, pp. 4-53, Aug. 2014.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(r))
    r=Constants.WGS84MeanRadius; 
end

Lr=predDist/r;
vertBias=sqrt(r^2*(1-cos(Lr))^2+(predDist-r*sin(Lr))^2);

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
