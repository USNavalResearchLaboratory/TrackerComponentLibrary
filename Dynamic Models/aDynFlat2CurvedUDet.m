function xDot=aDynFlat2CurvedUDet(x,t,aDyn,uBasis)
%%ADYNFLAT2CURVEDUDET Adapts a flat-Earth drift function to a curved Earth.
%         The position components of the state are kept in the global
%         coordinate system, whereas the other components are kept in the
%         local (flat Earth) coordinate system and the basis vectors
%         between them are assumed to be deterministically known as a
%         function of the state, and are given by uBasis. In the event that
%         the changes in the basis vectors depend on the changes in the
%         state, the function aDynFlat2CurvedUDyn can be used to map a
%         flat-Earth drift function to a curved Earth when given only the
%         derivatives of the basis vectors.
%
%INPUTS: x The target state such that the first three elements are position
%          in the global coordinate system, the next three elements are
%          velocity in the local, flat-Earth coordinate system, and all
%          other elements can be arbitrarily chosen in the local, flat-
%          earth coordinate system.
%        t The time at which the drift function should be evaluated.
%     aDyn The flat-Earth drift function of the form aDyn(x,t). It is
%          assumed that the first three components returned by aDyn are the
%          local velocity components of the state and that nothing in aDyn
%          is position-dependent.
%   uBasis The function to get the basis vectors. The function takes inputs
%          of the form uBasis(x,t)and returns the 9X9 set of coordinate
%          bases.
%
%OUTPUTS: xDot The derivative of x with respect to time. This maps the
%              drift function aDyn to a curved Earth.
%
%A discussion on mapping flat-Earth models to a curved Earth is given in
%[1].
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Transform the local target velocity into the global coordinate system.
rDot=getGlobalVectors(x(4:6,:),uBasis(x,t));

%Evaluate the flat-Earth drift function. The first three components
%returned will not be used.
xDot=aDyn(x,t);

xDot(1:3,:)=rDot;
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
