function aVal=aDynFlat2CurvedUDyn(xStacked,t,aDyn,uDyn)
%%ADYNFLAT2CURVEDUDYN Adapts a flat-Earth drift function to a curved
%                 Earth by jointly propagating the state with the
%                 derivatives of the basis vectors of a coordinate system
%                 that changes as the target moves. The position
%                 components of the state are kept in the global coordinate
%                 system, whereas the other components are kept in the
%                 local (flat Earth) coordinate system. In the event that
%                 the basis vectors are deterministically known as a
%                 function of the state, then use the function
%                 aDynFlat2CurvedUDet.
%
%INPUTS: xStacked The target state such that the first three elements are
%                 position in the global coordinate system, the next three
%                 elements are velocity in the local, flat-Earth coordinate
%                 system, and all other elements can be arbitrarily chosen
%                 in the local, flat-earth coordinate system. The final 9
%                 elements of xStacked are the three unit basis vectors.
%                 The first corresponds to the local x-axis, the second the
%                 local y- and the third the local z-axis. The global x, y
%                 and z position coordinates are the first three entries of
%                 xStacked.
%               t The time at which the drift function should be evaluated.
%            aDyn The flat-Earth drift function of the form aDyn(x,t),
%                 where is the non-stacked state. It is assumed that
%                 the first three components returned by aDyn are the local
%                 velocity components of the state and that nothing in aDyn
%                 is position-dependent.
%            uDyn The drift function of the basis vectors. The function
%                 takes inputs of the form uDyn(u,x,t), where u is the 9X9
%                 set of coordinate bases and x is the non-stacked state.
%                 The function returns a 9X9 matrix of derivatives with
%                 respect to time. If this input is omitted or an empty
%                 matrix is passed, then uDyn=@(u,x,t)uDotEllipsoid(u,x,t);
%                 will be used. That is the drift function for geodesic
%                 propagation.
%
%OUTPUTS: val The stacked derivatives of the state (global position
%             derivative, local other derivatives) and the local basis
%             vectors.
%
%A discussion on mapping flat-Earth models to a curved Earth is given in
%[1]. Note that when using a fully stochastic propagation model, one might
%need to account for the local/global coordiante system differences in the
%drift function as well. However, if linear models are used, then one can
%typically just use DPoly setting the numAug input to 9 (the number of
%dimensions of uDyn).
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(uDyn))
    uDyn=@(u,x,t)uDotEllipsoid(u,x);
end

%Extract the state and the basis vectors. The basis vectors are the final 9
%elements of xStacked.
xLen=length(xStacked)-9;
x=reshape(xStacked(1:xLen),xLen,1);
u=reshape(xStacked((xLen+1):end),3,3);

%Transform the local tangent velocity into the global coordinate system.
rDot=getGlobalVectors(x(4:6),u);

%Evaluate the flat-Earth drift function. The first three components
%returned will not be used.
xLDot=aDyn(x,t);

%Get the derivatives of the basis vectors.
uDot=uDyn(u,x,t);

aVal=[rDot(:);xLDot(4:end);uDot(:)];
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
