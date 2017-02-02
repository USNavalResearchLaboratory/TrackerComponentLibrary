function vObjRef=relVecAdd(vObsFrame,vObjInFrame)
%%RELVECADD  Special relativistic addition of velocity vectors. In the
%            inertial reference coordinate system, an observer moves with
%            constant velocity vObsFrame. In the coordinate system of the
%            observer, an object moves with constant velocity vObjInFrame.
%            This computes the velocity of the observed object in the
%            inertial reference frame. Under special relativity, it is not
%            just vObsFrame+vObjInFrame as it is in Newtonian mechanics.
%
%INPUTS:  vObsFrame The 3XN set of N velocity vectors in meters per second
%                   of the observer with respect to the inertial reference
%                   coordinate system. The magnitude of the velocity must
%                   be less than the speed of light.
%       vObjInFrame The 3XN set of N velocity vectors in meters per second
%                   of the object with respect to the observer's coordinate
%                   system. The magnitude of the velocity must be less than
%                   or equal to the speed of light.
%
%OUTPUTS: vObjRef  The 3XN set of velocity vectors in vObjInFrame
%                  transformed into the inertial reference coordinate
%                  system.
%
%The formulae for special relativistic velocity addition is derived in
%Chapter 1.4 of [1]. The magnitudes of vObsFrame and vObjInFrame must both
%be less than the speed of light (Constants.speedOfLight).
%
%REFERENCES:
%[1] G. Ludyk, Einstein in Matrix Form: Exact Derivation of the Theory of
%    Special and General Relativity without Tensors. Heidelberg: Springer,
%    2013.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

c=Constants.speedOfLight;

v=vObsFrame;
u=vObjInFrame;

%The magnitudes of all of the v vectors.
vMag=sqrt(sum(v.*v,1));

gamma=1./sqrt(1-vMag.^2/c^2);

%v^T*u for each of the velocity vectors.
uv=sum(v.*u,1);

Num=v+u+bsxfun(@times,(1./gamma-1),(u-bsxfun(@times,(uv./vMag.^2),v)));
Denom=1+uv/c^2;

vObjRef=bsxfun(@rdivide,Num,Denom);

%If vMag=0 for any of the vectors, then NaNs will appear. In such an
%instance, the correct solution is vObjInFrame, because the frame is not
%moving.
colSel=any(~isfinite(vObjRef),1);
vObjRef(:,colSel)=u(:,colSel);
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
