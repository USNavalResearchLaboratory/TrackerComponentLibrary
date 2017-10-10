function w=projVectorOnPlane(v,n)
%%PROJVECTORONPLANE Consider a vector v and a plane given by the equation
%                   n'*(x-p0)=0, where n is a normal to the plane, p0 is a
%                   point and x is the free variable. This function
%                   computes the projection of the vector onto the plane.
%                   THe quantity p0 is not needed.
%
%INPUTS: v A numDimX1 vector.
%        n The numDimX1 normal to the plane.
%
%OUTPUTS: w The numDimX1 projection of v onto n. This is notmal to n.
%
%The expression for this projection is given in Chapter 12.3 of [1].
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%x'*n+d=0;

u=n/(n'*n);
w=v-dot(v,n)*u;

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
