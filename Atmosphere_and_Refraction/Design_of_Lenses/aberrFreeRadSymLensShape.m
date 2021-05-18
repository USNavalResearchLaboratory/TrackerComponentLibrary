function [rb,zb,ra,za]=aberrFreeRadSymLensShape(zaFun,rSpan,n,t,ta,tb,numPoints)
%%ABERRFREERADSYMLENSSHAPE Given the shape of one side of a radially-
%          symmetric singlet lens, determine the shape of the other side so
%          as to be able to focus light from a given source to a specified
%          point without spherical aberration. The radial symmetry means
%          that one only needs to consider the shape of a single 2D slice
%          through the center of the lens, because the rest of the lens in
%          3D is given by rotating the slice about an axis coming out of
%          the center of the lens. Geometric optics is assumed.
%
%INPUTS: zaFun A handle to a function that takes the independent parameter
%              (which spans the width of the lens) and returns the height
%              of the left side of the lens and the derivative of the
%              height with respect to the independent parameter. The
%              function handle is called as [za,dza]=zaFun(ra) and it is
%              required that za=0 and dza=0 for ra=0. It is assumed that ra
%              can be a vector argument resulting in za and dza being
%              vector arguments.
%        rSpan A 2X1 or 1X2 vector giving the span of the values of ra that
%              are to go into zaFun. If a scalar is passed, then rSpan is
%              replaced with [-rSpan,rSpan]. This has the same units (e.g.
%              mm) as za, t, ta, and tb as well as the output parameters.
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
%           tb The location of the image focal point is (ra=0,t+tb), so tb
%              is the distance from the side of the lens opposite from
%              which ta is measured, to the image. Typically, tb is
%              positive.
%    numPoints The number of points within the span of rSpan over which the
%              surface of the lens are determined. The default if omitted
%              or an empty matrix is passed is 500.
%
%OUTPUTS: rb,zb These are 1XnumPoints vectors giving the Cartesian
%               locations of points of the side of the lens that was
%               computed by this function.
%         ra,za These are the coordinates of the side of the lens specified
%               by zaFun.
%
%This function implements the explicit solution of [1]. Note that the
%algorithm in the paper requires the choice of a sign parameter s1. This
%parameter is chosen to be the solution that gives the thickness t at ra=0.
%
%It is possible that rSpan is chosen to be so large that it exceeds the
%region over which a physical lens is possible. In such an instance, the
%lines for the top and bottom of the lens cross. Another unrealisable lens
%possibility arises when the surface for one side intersects itself. Often,
%that problem can be fixed by changing the index of refraction of the
%glass.
%
%EXAMPLE 1:
%This is an example of a lens with a parabolic surface, taken from Section
%3 of [1].
% n=1.5;
% t=8;
% ta=-60;
% tb=70;
% rSpan=16.63665027921;
% zaFun=@(ra)deal(ra.^2/200,2/200*ra);
% [rb,zb,ra,za]=aberrFreeRadSymLensShape(zaFun,rSpan,n,t,ta,tb);
% figure(1)
% clf
% hold on
% plot(ra,za,'-b','linewidth',2);
% plot(rb,zb,'-r','linewidth',4);
% h1=xlabel('r');
% h2=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%EXAMPLE 2:
%This is an example of a lens where one surace is the shape of a cosine.
%Such a surface is also given as an example in [1].
% n=1.5;
% t=8;
% ta=-60;
% tb=70;
% rSpan=16.585123012;
% zaFun=@(ra)deal(cos(ra/2),-(1/2)*sin(ra./2));
% [rb,zb,ra,za]=aberrFreeRadSymLensShape(zaFun,rSpan,n,t,ta,tb);
% figure(1)
% clf
% hold on
% plot(ra,za,'-b');
% plot(rb,zb,'-r');
% h1=xlabel('r');
% h2=ylabel('z');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] R. G. González-Acuña and H. A. Chaparro-Romo, "General formula for bi-
%    aspheric singlet lens design free of spherical aberration," Applied
%    Optics, vol. 57, no. 31, pp. 9341-9345, 1 Nov. 2018.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If the derivative of za at ra=0 is 0, then rb[ra]=0 for ra=0, so we can
%just choose the sign s1 that leads to the correct answer (The thickness
%must be t).
[za0,dza0]=zaFun(0);
[~,zb1]=aberrFreeLensSurfExplicit(0,za0,dza0,n,t,ta,tb,1);
[~,zbm1]=aberrFreeLensSurfExplicit(0,za0,dza0,n,t,ta,tb,-1);

if(abs(zb1-t)<abs(zbm1-t))
    s1=1;
else
    s1=-1;
end

if(nargin<7||isempty(numPoints))
    numPoints=500; 
end

if(isscalar(rSpan))
    rSpan=[-rSpan,rSpan]; 
end

ra=linspace(rSpan(1),rSpan(2),numPoints);
%Returns za and its derivative with respect to ra.
[za,dza]=zaFun(ra);
[rb,zb]=aberrFreeLensSurfExplicit(ra,za,dza,n,t,ta,tb,s1);

end

function [rb,zb]=aberrFreeLensSurfExplicit(ra,za,dza,n,t,ta,tb,s1)
%%ABERRFREELENSSURFEXPLICIT This function implements Equation 5 in [1] when
%               given ra, za and its derivative, dza. 
%
%REFERENCES:
%[1] González-Acuña and H. A. Chaparro-Romo, "General formula for bi-
%    aspheric singlet lens design free of spherical aberration," Applied
%    Optics, vol. 57, no. 31, pp. 9341-9345, 1 Nov. 2018.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

signta=sign(ta);

dza2p1=1+dza.^2;
razatadza=ra+(za-ta).*dza;
ra2zata2=ra.^2+(za-ta).^2;
sqrtra2zata2=sqrt(ra2zata2);

%Phi is Equation 3.
Phi=sqrt(1-(razatadza).^2./(n^2*ra2zata2.*dza2p1))./sqrt(dza2p1);

%ri is in Equation 3.
ri=razatadza./(n*sqrtra2zata2.*dza2p1)-dza.*Phi;

%zi is in Equation 3.
zi=((razatadza).*dza)./(n*sqrtra2zata2.*dza2p1)+Phi;

%fi is in Equation 6.
fi=ta-tb-signta.*sqrtra2zata2;

%h0 is in Equation 6.
h0=n*fi.*zi+ri.^2.*za-ra.*ri.*zi+(t+tb)*zi.^2-n^2*(za+t*zi);

%h1 is in Equation 6.
h1=ra.^2+2*t*ra.*ri+(tb-za).^2+t^2*(ri.^2+(zi-1).^2)-2*t*(tb-za).*(zi-1);

%zb is in Equation 5.
zb=(h0+s1*sqrt(zi.^2.*(fi.^2-2*n*fi.*(ra.*ri+t*ri.^2+zi.*(t*(zi-1)-tb+za))+n^2*h1-(ra.*zi+ri.*(t+tb-za)).^2)))/(1-n^2);

%rb is in Equation 5.
rb=ra+(ri.*(zb-za))./zi;

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
