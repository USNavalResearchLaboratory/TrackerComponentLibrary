function [rs,rp]=reflectionCoeffsVec(e1,e2,k,n,u1,u2)
%%REFLECTIONCOEFFSVEC Reflection coefficients for the reflection of an
%           electromagnetic plane wave from a LOSSLESS medium off a LOSSY
%           medium (or another lossless medium) are computed given a vector
%           in the direction of propagation of the light and a vector
%           normal to the surface. Both media must have real
%           permeabilities. For example, one might approximate the
%           troposphere as lossless while accounting for the loss in sea
%           water. Lossy media have complex permittivities. This function
%           is a different parameterization of reflectionCoeffs.
%
%INPUTS: e1 The permittivity (refraction index) of the lossless medium.
%           The incoming ray is traveling through this prior to reflecting
%           off of the surface with permittivity e2. This must be a real
%           quantity (lossless). e1 and e2 can both be either absolute
%           permittivities or relative permittivities. "Relative" means
%           that the permittivity of the medium has been divided by the
%           permittivity of free space and is a dimensionless quantity.
%        e2 The permittivity of the medium against which the ray reflects.
%           This can be complex (lossy).
%         k A 3X1 vector in the direction of propagation of the light
%           approaching the surface from which it will reflect.
%         n A 3X1 normal vector to the surface of reflection. It does not
%           matter whether n is pointing up or down from the surface.
%        u1 The permeability of the lossless medium. This must be
%           a real quantity.  u1 and u2 can both be either absolute
%           permeabilities or relative permeabilities. If omitted or an
%           empty matrix is passed, a value of 1 is used (the permeability
%           equals that of free space). That is a reasonable approximation
%           for air.
%        u2 The permeability of the lossy medium against which the incoming
%           ray reflects. medium. This must be a real quantity. If omitted
%           or an empty matrix is passed, a value of 1 is used. That is a
%           reasonable approximation for water.
%
%OUTPUTS: rs The complex reflection coefficient for s-polarized light, also
%            known as transverse-Electric (TE) polarized light.
%         rp The complex reflection coefficient for p-polarized light, also
%            known as transverse-Magnetic (TM) polarized light and as
%            tangent-plane polarized light.
%
%This function uses angBetweenVecs to get the angle between the incoming
%vector and the normal to the surface (adjusting in case one should have
%used -n instead of n) and then calls reflectionCoeffs.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(u2))
    u2=1;
end

if(nargin<5||isempty(u1))
    u1=1;
end

thetai=angBetweenVecs(k,n);

if(thetai>pi/2)
    %If -n should have been used instead of n.
    thetai=pi-thetai;
end
[rs,rp]=reflectionCoeffs(e1,e2,thetai,u1,u2);

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
