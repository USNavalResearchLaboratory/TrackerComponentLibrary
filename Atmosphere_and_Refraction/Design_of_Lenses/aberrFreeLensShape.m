function [xb,yb,zb,za]=aberrFreeLensShape(zaFun,xa,ya,n,T,fa,fb)
%%ABERRFREELENSSHAPE Given the shape of one side of a singlet lens,
%          determine the shape of the other side so as to be able to focus
%          light from a given source to a specified point without spherical
%          aberration.
%
%INPUTS: zaFun A handle to a function that takes the two independent
%              parameters, xa and ya, and returns the height of the left
%              side of the lens and the derivative of the height with
%              respect to each of the independent parameters.The function
%              handle is called as [za,dzadx,dzady]=zaFun(xa,ya) and it is
%              required that za=0 for xa=0, ya=0. It is assumed that xa and
%              ya can be same-sized matrix arguments resulting in za,
%              dzadx, and dzady being matrix arguments of the same size.
%        xa,ya A grid of points on parameterizing locations on the given
%              surface of the lens (they can go into zaFun). These matrices
%              must have the same shape.
%            n The index of refraction of the lends. The lens is assumed to
%              have a constant index of refraction and the surrounding
%              medium is assumed to have an index of refraction of 1. Thus,
%              if the surrounding medium does not have an index of
%              refraction of 1, then n should be the relative index of
%              refraction with respect to the medium.
%            T The thickness of the lens at xa=0,ya=0.
%           fa The object is located at (xa=0,,ya=0,fa), so fa is the
%              location of the object from one side of the lens. fa is
%              typically negative.
%           fb The location of the image focal point is (xa=0,ya=0,,T+fb),
%              so fb is the distance from the side of the lens opposite
%              from which fa is measured, to the image. Typically, fb is
%              positive.
%
%OUTPUTS: xb,yb,zb Points on the other side of the lens. These matrices
%                  have the same shape as xa and ya. These can eb used to
%                  plot a surface.
%               za The values returned by zaFun when given xa and ya.
%
%This function implements the explicit solution of [1]. Note that the
%expression for h in Equation 8b appears to be incorrect. Rather, the
%expression used in Appendix A is used in its place.
%
%EXAMPLE 1:
%This is one of the surfaces given in Section 3 of [1].
% n=1.5;
% T=1;
% fa=-5;
% fb=6;
% numPts=50;
% rspan=[-1.25,1.25];
% xySpan=linspace(rspan(1),rspan(2),numPts);
% [xa,ya]=meshgrid(xySpan,xySpan);
% zaFun=@(xa,ya)deal((xa.^2+8*ya.^2)/200,xa/100,2*ya/25);
% [xb,yb,zb,za]=aberrFreeLensShape(zaFun,xa,ya,n,T,fa,fb);
% 
% figure(1)
% clf
% hold on
% surf(xb,yb,zb)
% surf(xa,ya,za)
% view(45,30)
% h1=xlabel('x');
% h2=ylabel('y');
% h3=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%This is the example used in Appendix A of [1].
% n=1.5;
% T=10;
% fa=-50;
% fb=60;
% numPts=50;
% rspan=[-12.5,12.5];
% xySpan=linspace(rspan(1),rspan(2),numPts);
% [xa,ya]=meshgrid(xySpan,xySpan);
% zaFun=@(xa,ya)deal(cos((2/5)*xa).*cos((2/5)*ya),-(2/5)*cos((2*ya)/5).*sin((2*xa)/5),-(2/5)*cos((2*xa)/5).*sin((2*ya)/5));
% [xb,yb,zb,za]=aberrFreeLensShape(zaFun,xa,ya,n,T,fa,fb);
% 
% figure(1)
% clf
% hold on
% surf(xb,yb,zb)
% surf(xa,ya,za)
% view(60,20)
% h1=xlabel('x');
% h2=ylabel('y');
% h3=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] R. G. González-Acuña, H. A. Chaparro-Romo, and J. C. Gutiérrez-Vega,
%    "General formula to design a freeform singlet free of spherical
%    aberration and astigmatism," Applied Optics, vol. 58, no. 4, pp. 1010-
%    1015, 1 Feb. 2019.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[za,dzadx,dzady]=zaFun(xa,ya);

%S is defined in Equation 3.
S=sqrt(1+dzadx.^2+dzady.^2);

%D is defined in Equation 3.
D=sqrt(xa.^2+ya.^2+(za-fa).^2);

%Phi is defined in Equation 5.
Phi=sqrt(1-((ya.*dzadx-xa.*dzady).^2+(xa+dzadx.*(za-fa)).^2+(ya+dzady.*(za-fa)).^2)./(n^2.*D.^2.*S.^2));

%X is defined in Equation 4a
X=(xa.*(1+dzady.^2)-dzadx.*(fa+ya.*dzady-za))./(n.*D.*S.^2)-(dzadx.*Phi)./S;

%Y is defined in Equation 4b.
Y=(ya.*(1+dzadx.^2)-dzady.*(fa+xa.*dzadx-za))./(n.*D.*S.^2)-(dzady.*Phi)./S;

%Z is defined in Equation 4c.
Z=(xa.*dzadx+ya.*dzady+(dzadx.^2+dzady.^2).*(za-fa))./(n*D.*S.^2)+Phi./S;

%p is defined in Equation 8c.
p=-sign(fa).*D+fa-fb;

%q is defined in Equation 8d.
q=n.^2.*T-n*p+xa.*X+ya.*Y;

%g is defined in Equation 8a
g=q.*Z+(n.^2-1).*za+Z.^2.*(za-T-fb);

%h is defined in Equation 8b. However, the expression here is a corrected
%form mirroring more closely what is in Appendix A.
h=(X.^2+Y.^2-n.^2).*za.^2-2.*q.*za.*Z+((T+fb).^2+xa.^2+ya.^2-(p-n.*T).^2).*Z.^2;

%zb is defined in Equation 7.
zb=(g-sqrt(g.^2+(n.^2-1).*h))./(n.^2-1);
%xb is defined in Equation 7.
xb=xa+(X.*(zb-za))./Z;
%yb is defined in Equation 7.
yb=ya+(Y.*(zb-za))./Z;

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
