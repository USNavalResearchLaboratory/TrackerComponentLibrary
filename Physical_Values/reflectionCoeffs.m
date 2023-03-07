function [rs,rp]=reflectionCoeffs(e1,e2,thetai,u1,u2)
%%REFLECTIONCOEFFS Reflection coefficients for the reflection of an
%           electromagnetic plane wave from a LOSSLESS medium into a LOSSY
%           medium (or another lossless medium) are computed based on the
%           angle that the incoming ray makes with a normal to the surface.
%           Both media must have real permeabilities. For example, one
%           might approximate the troposphere as lossless while accounting
%           for the loss in sea water. Lossy media have complex
%           permittivities.
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
%    thetai The (acute) incident angle of the incoming ray with respect to
%           a normal to the surface. 0 means that it is coming in 90
%           degrees to the surface and pi/2 means that it is parallel to
%           the surface and never crosses. This should be in the range
%           0<=thetai<=pi/2.
%        u1 The permeability of the lossless medium. This must be
%           a real quantity.  u1 and u2 can both be either absolute
%           permeabilities or relative permeabilities. If omitted or an
%           empty matrix is passed, a value of 1 is used (the permeability
%           equals that of free space). That is a reasonable approximation
%           for air.
%        u2 The permeability of the lossy medium against which the incoming
%           ray reflects. This must be a real quantity. If omitted or an
%           empty matrix is passed, a value of 1 is used. That is a
%           reasonable approximation for water.
%
%OUTPUTS: rs The complex reflection coefficient for s-polarized light, also
%            known as transverse-Electric (TE) polarized light.
%         rp The complex reflection coefficient for p-polarized light, also
%            known as transverse-Magnetic (TM) polarized light and as
%            tangent-plane polarized light.
%
%This function implements Equations 5 and 10 of [1], substituting the
%complex angle expression in 12. Unlike in [1], we expanded epsilon_1 and
%mu_1 in terms of a relative permittivity times the permittivity of free
%space and a relative permeability times the permeability of free space.
%that allows for cancellation with the expanded form of beta_2 when
%expressed in such a manner.
%
%The results of [1] show that the corrected solution in [2] is equivalent
%to a previously ad-hoc approach when using a "complex" angle. The topic
%has been somewhat controversial, with additional discussion in [3] and
%[4].
%
%EXAMPLE:
%We plot the magnitude and phase of the reflection coefficients between air
%and seawater at 13MHz frequency.
% e1=1;%Assume that the relative permittivity of air is 1.
% T=20;
% S=35;%35 parts per thousand.
% v=13e6;
% e2=seawaterRelPermittivity(T,v,S);
% 
% %We take the relative permeabilities to be 1.
% numPts=1000;
% thetai=linspace(0,pi/2,numPts);
% 
% rs=zeros(numPts,1);
% rp=zeros(numPts,1);
% for k=1:numPts
%     [rs(k),rp(k)]=reflectionCoeffs(e1,e2,thetai(k));
% end
% 
% figure(1)
% clf
% hold on
% plot(rad2deg(pi/2-thetai),abs(rs),'-r','linewidth',2');
% plot(rad2deg(pi/2-thetai),abs(rp),'-k','linewidth',2');
% xlabel('Incident Angle Above Plane (\circ)')
% ylabel('Reflection Coefficient Magnitude')
% legend('s-Polarization','p-Polarization','location','southeast')
% 
% figure(2)
% clf
% hold on
% plot(rad2deg(pi/2-thetai),rad2deg(angle(rs)),'-r','linewidth',2');
% plot(rad2deg(pi/2-thetai),rad2deg(angle(rp)),'-k','linewidth',2');
% xlabel('Incident Angle Above Plane (\circ)')
% ylabel('Reflection Coefficient Phase (\circ)')
% legend('s-Polarization','p-Polarization','location','southeast')
%
%REFERENCES:
%[1] I. M. Besieris, "Comment on the "Corrected Fresnel coefficients for
%    lossy materials,"" IEEE Antennas and Propagation Magazine, vol. 53,
%    no. 4, pp. 161-164, Aug. 2011.
%[2] F. X. Canning, "Corrected Fresnel coefficients for lossy materials,"
%    in Proceedings of the IEEE Interntional Symposium on Antennas and
%    Propagation, Spokane, WA, 3-8 Jul. 2011, pp. 2123-2126.
%[3] F. X. Canning, "On the Fresnel coefficients for transmission into a
%    lossy medium," in Proceedings of the International Conference on
%    Electromagnetic in Advanced Applications, Turin, Italy, 7-1 Sep. 2015,
%    pp. 165-168.
%[4] M. Oh, "Complex unit vector for the complex wave constant k in
%    a lossy medium," IEEE Antennas and Propagation Magazine, vol. 63,
%    no. 1, pp. 117-120, Feb. 2021.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(u2))
    u2=1;
end

if(nargin<4||isempty(u1))
    u1=1;
end

cosTheta=cos(thetai);
sinTheta=sin(thetai);

sqrtTerm=sqrt(e2*u2-e1*u1*sinTheta^2);

%p polarization is transverse-Magnetic (TM) polarization.
%s-polarization is transverse-Electric (TE) polarization.
rs=(u2*sqrt(e1*u1)*cosTheta-u1*sqrtTerm)/(u2*sqrt(e1*u1)*cosTheta+u1*sqrtTerm);
rp=(e2*sqrt(e1*u1)*cosTheta-e1*sqrtTerm)/(e2*sqrt(e1*u1)*cosTheta+e1*sqrtTerm);

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
