function resp=stdCircAperBeamPattern(u,R,weighting)
%%STDCIRCAPERBEAMPATTERN Obtain the value of the beampattern for a
%       continuous circular aperture having one of a number of circularly
%       symmetric taperings. As defined in Chapter 2.2 of [1], the
%       beampattern is the frequency-wavenumber response function evaluated
%       versus direction. It describes the complex gain of the aperture to
%       an input plane wave. The values returned by this function are real.
%
%INPUTS: u This can be a 2XN set of N look directions from the aperture
%          plane in terms of direction cosines, in which case sum(uv.^2,1)
%          should all be less than or equal to 1; larger values do not
%          cause a warning or error. Otherwise, this can be a 1XN set of
%          sines of the angle offsets from boresight. Due to the circular
%          symmetry, only the offset from boresight matters.
%        R The radis of the aperture in units of wavelengths.
% weighting A parameter specifying the type of circularly symmetric radial
%          tapering to use across the aperture. Possible values are
%          0 (The default if omitted or an empty matrix is passed). Use a
%            uniform taper of 1.
%          1 Use a taper of (1-(r/R)^2) where r is the distance from the
%            center.
%          2 Use a taper of (1-(r/R)^2)^2 where r is the distance from the
%            center.
%
%OUTPUTS: B The beampattern evaluated at the points in u.
%
%This function implements the formulae in Table 4.4 of [1].
%
%EXAMPLE 1:
%Here, we will plot the beam pattern magnitude power in decibels over the
%viewable region.
% R=4;
% numPoints=100;
% uVals=linspace(-1,1,numPoints);
% [u,v]=meshgrid(uVals,uVals);
% uv=[u(:).';v(:).'];
% B=stdCircAperBeamPattern(uv,R,0);
% %Set the values of the response outside of the valid region to zero.
% sel=sum(uv.^2,1)>1;
% B(sel)=0;
% B=reshape(B,numPoints,numPoints);
% 
% figure(1)
% clf
% surf(u,v,20*log10(abs(B)),'EdgeColor','None')
% axis([-1 1 -1 1 -50 0])
% caxis([-50 0])
% colormap(jet(256))
% view(20,45)
% light()
% h1=xlabel('u');
% h2=ylabel('v');
% h3=zlabel('Array Gain');
% title('Array Power Gain in Decibels')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%This reproduces Figure 4.43 in [1] using a different parameter for the
%horizontal axis.
% R=10;
% numPoints=1000;
% uVals=linspace(0,1,numPoints);
% B=stdCircAperBeamPattern(uVals,R,0);
% 
% figure(1)
% plot(uVals,20*log10(abs(B)))
% axis([0 1 -60 0])
% h1=xlabel('u');
% h2=ylabel('B Power (dB)');
% set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] H. L. Van Trees, Optimum Array Processing. New York:
%    Wiley-Interscience, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(weighting))
   weighting=0; 
end

%If we were given u-v points, turn then directly into the sine of the angle
%offset from boresight.
if(size(u,1)==2)
    u=sqrt(abs(sum(u.^2,1)));
end

%Defined in Equation 4.143 of [1].
psiR=2*pi*R*u;

switch(weighting)
    case 0%Uniform
        resp=2*besselj(1,psiR)./psiR;
    case 1%Radial taper
        resp=8*besselj(2,psiR)./psiR.^2;
    case 2%Radial taper squared
        resp=48*besselj(3,psiR)./psiR.^3;
    otherwise
        error('Unknown weighting specified.')
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
