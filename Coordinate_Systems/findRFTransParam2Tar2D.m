function M=findRFTransParam2Tar2D(lRx,lTar)
%%FINDRFTRANSPARAM2TAR2D Find the rotation matrix needed to rotate a
%            Cartesian vector to local radar-facing coordinates such
%            that a specified target is on the boresight of the radar in
%            2D. This makes the local y-axis of the radar point at the
%            target.
%
%INPUTS: lRx The 2X1 Cartesian location of the receiver.
%       lTar The 2X1 Cartesian location of the target.
%
%OUTPUTS: M A 2X2 rotation matrix for the transformation from 2D global 
%           Cartesian coordinates to local radar-facing coordinates. The y-
%           axis represents the pointing direction of the radar. This can
%           be directly fed into the 2D RU coordinate transform functions.
%
%EXAMPLE:
% lRx=[10e3;-15e3];
% lTar=[20e3;5e3];
% M=findRFTransParam2Tar2D(lRx,lTar);
% %One sees that the u component is zero.
% u=getUDirection2D(lTar,lRx,M)
%
%April 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

theta=atan2(lTar(1)-lRx(1),lTar(2)-lRx(2));
M=rotMat2D(theta,false);

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
