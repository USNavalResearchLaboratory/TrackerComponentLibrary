function points=regularIcosahedronCoords()
%%REGULARICOSAHEDRONCOORDS Obtain the coordinates of the vertices of a
%               regular icosahedron. An icosahedron is a polyhedron having
%               20 faces. The regular icosahedron is a platonic solid
%               having 20 equivalent equilateral triangle faces, where at
%               each verted, five edges meet. The vertices are all places
%               on the surface of the unit sphere.
%
%INPUTS: None
%
%OUTPUTS: points 3X12 set of points that form the 12 vertices of the
%                regular icosahedron.
%
%A formula for the vertices of the regular icosahedron is given in Chapter
%9 of [1]. However, rather than using that, the regular icosahedron
%vertices from formula U3 5-1 in [1] are used as they provide unit-distance
%vertices.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

r=sqrt((5+sqrt(5))/10);
s=sqrt((5-sqrt(5))/10);

points=zeros(3,12);

points(:,1)=[r;s;0];
points(:,2)=[r;-s;0];
points(:,3)=[-r;s;0];
points(:,4)=[-r;-s;0];

points(:,5)=[0;r;s];
points(:,6)=[0;r;-s];
points(:,7)=[0;-r;s];
points(:,8)=[0;-r;-s];

points(:,9)=[s;0;r];
points(:,10)=[s;0;-r];
points(:,11)=[-s;0;r];
points(:,12)=[-s;0;-r];

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
