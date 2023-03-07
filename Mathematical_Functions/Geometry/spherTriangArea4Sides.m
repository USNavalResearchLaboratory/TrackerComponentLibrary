function surfArea=spherTriangArea4Sides(a,b,c,r)
%%SPHERTRIANGAREA4SIDES Given the length of the sides of a spherical
%           triangle as well as the radius of the spher,e find the surface
%           area of the triangle. Note that a,b, and c should be consistent
%           with a spherical triangle and using inconsistent values (even
%           if much less than the circumference of the sphere) can result
%           in (invalid) complex solutions). This is similar to the
%           sideAngles2SpherTriangArea function.
%
%INPUTS: a, b, c The lengths of the sides of the spherical triangle
%                in any order. If one wishes to convert multiple triangles
%                at once, these can be matrices.
%              r The scalar radius of the sphere on which the distances are
%                measured. If omitted or an empty matrix is passed, the
%                default of 1 is used.
%
%OUTPUTS: surfArea The surface area of the spherical triangle. The
%                  dimensions of this correspond to those of the inputs.
%
%This function just calls spherTriAngles4Sides and then
%spherTriangArea4Angles.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(r))
    r=1;
end

[A,B,C]=spherTriAngles4Sides(a,b,c,r);
surfArea=spherTriangArea4Angles(A,B,C,r);

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