function [zCent,r]=circumscribedSphere4Tetrahedron(v)
%%CIRCUMSCRIBEDSPHERE4TETRAHEDRON Find the center and radius of a
%       circumscribed sphere for a given tetrahedron. A circumscribed
%       sphere is one such that all of the vertices of the tetrahedron
%       touch the sphere.
%
%INPUTS: v A 3X4 set of the 4 3D vertices of the tetrahedron.
%
%OUTPUTS: zCent The 3X2 center of the circumscribed sphere of the
%               tetrahedron.
%             r The radius of the circumscribed sphere.
%
%The formulae for the center and radius of a sphere circumscribed in a
%tetrahedron are given in [1].
%
%EXAMPLE:
%In this example, a tetrahedron is drawn and the circumscribing sphere is
%also plotted.
% v=100*randn(3,4);
% figure(1)
% clf
% hold on
% [zCent,r]=circumscribedSphere4Tetrahedron(v);
% s=drawEllipse(zCent,eye(3,3),r^2);
% alpha(s{1},0.25);
% scatter3(v(1,:),v(2,:),v(3,:),400,'.b')
% for i1=1:4
%     for i2=(i1+1):4
%         plot3([v(1,i1),v(1,i2)],[v(2,i1),v(2,i2)],[v(3,i1),v(3,i2)],'-b','linewidth',2)
%     end
% end
% view(100,40)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Circumsphere." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/Circumsphere.html
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Equation 3
a=det([v',ones(4,1)]);

sqCol=sum(v.^2,1)';

%Equation 5
Dx=det([sqCol,v(2,:)',v(3,:)',ones(4,1)]);
%Equation 6
Dy=-det([sqCol,v(1,:)',v(3,:)',ones(4,1)]);
%Equation 7
Dz=det([sqCol,v(1,:)',v(2,:)',ones(4,1)]);

c=det([sqCol,v']);

zCent=[Dx;Dy;Dz]/(2*a);
r=sqrt(Dx^2+Dy^2+Dz^2-4*a*c)/(2*abs(a));

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
