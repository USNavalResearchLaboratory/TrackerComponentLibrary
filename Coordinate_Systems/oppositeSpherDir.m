function azEl=oppositeSpherDir(azEl,systemType)
%%OPPOSITESPHERDIR Given spherical angles specifying a direction, adjust
%         them to represent the opposite direction.
%
%INPUTS: azEl A 2XN set of [azimuth;elevation] pais in radians to flip.
%  systemType This specified the type of spherical cooridnate system used.
%             The way to slip the direction is the same for all except for
%             systemType=2. Possible values are the same as in Cart2Sphere
%             and are:
%           0 (The default if omitted) Azimuth is measured
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one desires the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%
%OUTPUTS: azEl A 2XN set of [azimuth;elevation] pairs representing the
%              opposite directions of those on the input.
%
%EXAMPLE:
%This generates 10 random directions. It then converts then into an azimuth
%and elevation. The azimuth and elevation are changed to represent the
%opposite direction and are then converted into a direction vector. Lines
%from the origin in each direction are plotted. The original directions are
%red and the flippsed ones are blue. One can see that the directions have
%been flipped.
% N=10;
% numDim=3;
% uStart=randDirVec(numDim,N);
% systemType=0;
% useHalfRange=true;
% azEl=Cart2Sphere(uStart,systemType,useHalfRange);
% azEl=azEl(2:3,:);
% azElFlip=oppositeSpherDir(azEl,systemType);
% uEnd=spher2Cart(azElFlip,systemType,useHalfRange);
% 
% figure(1)
% clf
% hold on
% for k=1:N
%     plot3([0,uStart(1,k)],[0,uStart(2,k)],[0,uStart(3,k)],'-r','linewidth',2)
%     plot3([0,uEnd(1,k)],[0,uEnd(2,k)],[0,uEnd(3,k)],'-b','linewidth',2)
% end
% axis equal
%
%March 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(systemType))
    systemType=0;
end

azEl(1,:)=wrapRange(azEl(1,:)+pi,-pi,pi);
if(systemType==2)
    azEl(2,:)=pi-azEl(2,:);
else
    azEl(2,:)=-azEl(2,:);
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
