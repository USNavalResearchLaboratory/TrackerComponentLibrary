function points=regularPGonCoords(p)
%%PGONCOORDS Find the coordinates of a regular polygon with p sides. A
%            polygon is regular if all sides have equal length and if the
%            angles of each vertex are all equal. The polygon returned by
%            this function is centered at the origin and all points are
%            unit distance from the origin. The polygon is unique within a
%            rotation about the origin.
%
%INPUTS: p The number of sides of the polygon. This must be >=3.
%
%OUTPUTS: points A 2Xp set of the vertices of the polygon.
%
%The formula for the vertices of a regular polygon is from Chapter 9 of
%[1].
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

k=0:(p-1);
points=[cos(2*pi*k/p);sin(2*pi*k/p)];
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
