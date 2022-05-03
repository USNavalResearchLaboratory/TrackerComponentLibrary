function azEqPts=Cart2AzEquidistantProj(xCart,latLonRef,a,f)
%%CART2AZEQUIDISTANTPROJ Given points in 3D East centered-Earth fixed
%       coordinates convert them to 3D points in an azimuthal equidistant
%       projection about a specified point.
%
%INPUTS: xCarts A 3XN set of N 3D Cartesian points to convert into azimuthal
%               equidistant coordinates.
%     latLonRef A 2X1 [latitude;longitude] reference point in raidans
%               about which the projection is taken.
%             a The semi-major axis of the reference ellipsoid (in
%               meters). If this argument is omitted or an empty matrix
%               is passed, the value in Constants.WGS84SemiMajorAxis is
%               used.
%             f The flattening factor of the reference ellipsoid. If
%               this argument is omitted or an empty matrix is passed,
%               the value in Constants.WGS84Flattening is used.
%
%OUTPUTS: azEqPoints A 3XN set of the points converted to an azimuthal
%                   equidistant projection about latLonRef.
%
%This function just calls Cart2Ellipse and then ellips2AzEquidistantProj.
%
%EXAMPLE:
%This example just demonstrates that the azEquidistantProj2Cart and
%Cart2AzEquidistantProj are consistent with each other. Converting a
%Cartesian point to azimuthal equidistant coordinates and then back
%produces a value whose relative error from the original point is on the
%order of finite precision limitiations. 
% latLonRef=deg2rad([20.906029;-157.078918]);
% a=Constants.WGS84SemiMajorAxis;
% f=Constants.WGS84Flattening;
% llhPt=[[deg2rad([20.631756;-155.551818]);0],[deg2rad([21.295516;-158.002386]);10e3]];
% xCart=ellips2Cart(llhPt,a,f);
% azEqCoords=Cart2AzEquidistantProj(xCart,latLonRef,a,f);
% xCartBack=azEquidistantProj2Cart(azEqCoords,latLonRef,a,f);
% RelErr=max(max(abs((xCartBack-xCart)./xCart)))
%
%January 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

llhPts=Cart2Ellipse(xCart,[],a,f);
azEqPts=ellips2AzEquidistantProj(llhPts,latLonRef,a,f);

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
