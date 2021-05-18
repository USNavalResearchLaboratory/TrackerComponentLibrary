function [rb,zb,ra,za,hb,tb]=aberrFreeRadSymLens4Surf(zaFun,hbFun,tbFun,rSpan,n,t,ta,numPoints)
%%ABERRFREERADSYMLENS4SURF Given the shape of one side of a radially-
%          symmetric singlet lens, determine the shape of the other side so
%          as to be able to focus light from a given source to a specified
%          distant surface without spherical aberration. The radial
%          symmetry means that one only needs to consider the shape of a
%          single 2D slice through the center of the lens (and the surface,
%          because the rest of the lens and surface in 3D are given by
%          rotating the slice about an axis coming out of the center of the
%          lens. Geometric optics is assumed.
%
%INPUTS: zaFun A handle to a function that takes the independent parameter
%              (which spans the width of the lens) and returns the height
%              of the left side of the lens and the derivative of the
%              height with respect to the independent parameter. The
%              function handle is called as [za,dza]=zaFun(ra) and it is
%              required that za=0 and dza=0 for ra=0. It is assumed that ra
%              can be a vector argument resulting in za and dza being
%              vector arguments.
%  hbFun,tbFun These functions handles are called as hb=hfFun(ra), and
%              tb=tbFun(ra). The points (hb,tb) are on the surface onto
%              which one wishes to project from the lens, without
%              aberration.
%        rSpan A 2X1 or 1X2 vector giving the span of the values of ra that
%              are to go into zaFun. If a scalar is passed, then rSpan is
%              replaced with [-rSpan,rSpan]. This has the same units (e.g.
%              mm) as za, t, ta, hb, and tb as well as the output
%              parameters.
%            n The index of refraction of the lends. The lens is assumed to
%              have a constant index of refraction and the surrounding
%              medium is assumed to have an index of refraction of 1. Thus,
%              if the surrounding medium does not have an index of
%              refraction of 1, then n should be the relative index of
%              refraction with respect to the medium. n can be negative but
%              not zero.
%            t The thickness of the lens at ra=0; t>0. 
%           ta The object is located at (ra=0,ta), so ta is the location of
%              the object from one side of the lens. ta is typically
%              negative.
%    numPoints The number of points within the span of rSpan over which the
%              surface of the lens are determined. The default if omitted
%              or an empty matrix is passed is 500.
%
%OUTPUTS: rb,zb These are 1XnumPoints vectors giving the Cartesian
%               locations of points of the side of the lens that was
%               computed by this function.
%         ra,za These are the coordinates of the side of the lens specified
%               by zaFun.
%         hb,tb These are the coordinates of the surface specified by hbFun
%               and tbFun.
%
%This function implements the explicit solution of [1]. Note that the
%algorithm in the paper requires the choice of a sign parameter s1. This
%parameter is chosen to be the solution that gives the thickness t at ra=0.
%
%Note that physically irrealizable systems might produce lends deisgns with
%surfaces that intersect themselves.
%
%EXAMPLE 1:
%This is the example in Figure 4b of [1]. From the figure, it is clear that
%it was meant for ta to be -50, not 50.
% n=1.5;
% t=8;
% ta=-50;
% zaFun=@(ra)deal(ra.^2/100,ra/50);
% tbFun=@(ra)60+0.0008*ra.^2;
% hbFun=@(ra)0.04*ra;
% 
% rSpan=15;
% [rb,zb,ra,za,hb,tb]=aberrFreeRadSymLens4Surf(zaFun,hbFun,tbFun,rSpan,n,t,ta);
% figure(1)
% clf
% hold on
% plot(ra,za,'-b','linewidth',2);%The given side of the lens.
% plot(rb,zb,'-r','linewidth',4);%The computed side of the lens.
% plot(hb,t+tb,'-g','linewidth',2);%The surface onto which it is projected.
% h1=xlabel('r');
% h2=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%This is the example in Figure 3a of [1].
% n=1.5;
% t=8;
% ta=-50;
% zaFun=@(ra)deal(cos((2/5)*ra),-(2/5)*sin((2/5)*ra));
% tbFun=@(ra)20+0.4*cos(0.4*ra);
% hbFun=@(ra)ra/2;
% 
% rSpan=15;
% [rb,zb,ra,za,hb,tb]=aberrFreeRadSymLens4Surf(zaFun,hbFun,tbFun,rSpan,n,t,ta);
% figure(1)
% clf
% hold on
% plot(ra,za,'-b','linewidth',2);%The given side of the lens.
% plot(rb,zb,'-r','linewidth',4);%The computed side of the lens.
% plot(hb,t+tb,'-g','linewidth',2);%The surface onto which it is projected.
% h1=xlabel('r');
% h2=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] R. G. González-Acuña, M.Avendaño-Alejo,and J. C. Gutiérrez-Vega,
%    "Singlet lens for generating aberration-free patterns on deformed
%    surfaces," Journal of the Optical Society of America A, vol. 36, no.
%    5, pp. 925-929, May 2019.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If the derivative of za at ra=0 is 0, then rb[ra]=0 for ra=0, so we can
%just choose the sign s1 that leads to the correct answer (The thickness
%must be t).
[za0,dzadra0]=zaFun(0);
tb0=tbFun(0);
hb0=hbFun(0);%This should typically be 0.

[~,zb1]=aberrFreeRadSymLens4SurfExplicit(0,za0,dzadra0,n,t,ta,hb0,tb0,1);
[~,zbm1]=aberrFreeRadSymLens4SurfExplicit(0,za0,dzadra0,n,t,ta,hb0,tb0,-1);

if(abs(zb1-t)<abs(zbm1-t))
    s1=1;
else
    s1=-1;
end

if(nargin<8||isempty(numPoints))
    numPoints=500; 
end

if(isscalar(rSpan))
    rSpan=[-rSpan,rSpan]; 
end

ra=linspace(rSpan(1),rSpan(2),numPoints);

[za,dzadra]=zaFun(ra);
tb=tbFun(ra);
hb=hbFun(ra);

[rb,zb]=aberrFreeRadSymLens4SurfExplicit(ra,za,dzadra,n,t,ta,hb,tb,s1);
end

function [rb,zb]=aberrFreeRadSymLens4SurfExplicit(ra,za,dzadra,n,t,ta,hb,tb,s1)
%%ABERRFREELENSSURFEXPLICIT This function implements Equation 5 in [1] when
%               given ra, za and its derivative, dza, as well as hb and tb.
%
%REFERENCES:
%[1] R. G. González-Acuña, M.Avendaño-Alejo,and J. C. Gutiérrez-Vega,
%    "Singlet lens for generating aberration-free patterns on deformed
%    surfaces," Journal of the Optical Society of America A, vol. 36, no.
%    5, pp. 925-929, May 2019.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

%f is defined in Equation 6a.
f=-ta+n*t+sqrt(hb.^2+tb.^2)+sign(ta)*sqrt(ra.^2+(ta-za).^2);

%Z is defined in Equation 3.
Z=(dzadra.*(ra+(za-ta).*dzadra))./(n*sqrt(ra.^2+(za-ta).^2).*(1+dzadra.^2))+sqrt(1-(ra+(za-ta).*dzadra).^2./(n^2*(ra.^2+(za-ta).^2).*(1+dzadra.^2)))./sqrt(1+dzadra.^2);

%R is defined in Equation 3.
R=(ra+(za-ta).*dzadra)./(n*sqrt(ra.^2+(za-ta).^2).*(1+dzadra.^2))-(dzadra.*sqrt(1-(ra+(za-ta).*dzadra).^2./(n^2*(ra.^2+(za-ta).^2).*(1+dzadra.^2))))./sqrt(1+(dzadra).^2);

%Note that Z.^2+R.^2 should =1.

%zd is given in Equation 5
zb=(1/(n.^2-1)).*((n-R).*(n+R).*za+n*f.*Z+(ra-hb).*R.*Z-(t+tb).*Z.^2+(s1/2)*sqrt(4*(n.^2*za-R.^2.*za+n*f.*Z+(ra-hb).*R.*Z-(t+tb).*Z.^2).^2-4*(1-n.^2).*(-n.^2*za.^2+R.^2.*za.^2-2*n*f.*za.*Z+2*(hb-ra).*R.*za.*Z+((hb-ra).^2+(t+tb).^2-f.^2).*Z.^2)));

%rb is given in Equation 5.
rb=ra+(R.*(zb-za))./Z;
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
