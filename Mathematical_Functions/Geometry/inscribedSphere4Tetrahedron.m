function [zCent,r]=inscribedSphere4Tetrahedron(v)
%%INSCRIBEDPHERE4TETRAHEDRON Find the center and radius of an
%       inscribed sphere for a given tetrahedron. An inscribed sphere is
%       the smallest sphere that touches all of the faces of the
%       tetrahedon.
%       
%INPUTS: v A 3X4 set of the 4 3D vertices of the tetrahedron.
%
%OUTPUTS: zCent The 3X2 center of the inscribed sphere of the tetrahedron.
%             r The radius of the inscribed sphere.
%
%While formulae are given in [1], it is not actually necessary to solve a
%linear set of equations to get the center of the inscribed circle. In [2],
%it is noted that the inscribed sphere touches each face of the tetrahedron
%at the centroid of each face. A similar relation is given in 4 of lemma 4
%of [3].
%
%EXAMPLE:
%Here, the sides of a tetrahedron are drawn as well as the inscribed
%sphere. Depending on the dimensions of the tetrahdron, the sphere might be
%small or large.
% v=100*randn(3,4);
% figure(1)
% clf
% hold on
% [zCent,r]=inscribedSphere4Tetrahedron(v);
% drawEllipse(zCent,eye(3,3),r^2);
% scatter3(v(1,:),v(2,:),v(3,:),400,'.b')
% for i1=1:4
%     for i2=(i1+1):4
%         plot3([v(1,i1),v(1,i2)],[v(2,i1),v(2,i2)],[v(3,i1),v(3,i2)],'-b','linewidth',2)
%     end
% end
% view(100,10)
% axis equal
%
%REFERENCES:
%[1] P. P. Klein, "The insphere of a tetrahedron," Applied Mathematics,
%    vol. 11, no. 7, pp. 601-612, Jul. 2020.
%[2] J. Burkardt, "Computational geometry lab: TETRAHEDRONS,"
%    Virginia Tech, Tech. Rep., 23 Dec. 2010. [Online]. Available:
%    https://people.sc.fsu.edu/~jburkardt/classes/cg_2007/cg_lab_tetrahedrons.pdf
%[3] M. Pollicott, "Limit points for tetrahedra and inscribed spheres,"
%    University of Warwick, Tech. Rep., undated, Accessed 12/6/2021.
%    [Online]. Available: https://homepages.warwick.ac.uk/~masdbl/inscribed.pdf
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

a=v(:,1);
b=v(:,2);
c=v(:,3);
d=v(:,4);

alphaVal=det([v',ones(4,1)]);
Nabc=norm(cross(b-a,c-a));
Nabd=norm(cross(b-a,d-a));
Nacd=norm(cross(c-a,d-a));
Nbcd=norm(cross(c-b,d-b));

denom=Nabc+Nabd+Nacd+Nbcd;

zCent=(Nabc*d+Nabd*c+Nacd*b+Nbcd*a)/denom;
r=abs(alphaVal)/denom;

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
