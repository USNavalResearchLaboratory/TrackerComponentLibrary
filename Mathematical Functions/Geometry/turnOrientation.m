function turnDir=turnOrientation(v1,v2,v3)
%%TURNORIENTATION Given 3 two-dimensional vertices in order, determine
%                 whether they form a left turn (go counterclockwise).
%                 Such a determination is necessary for determining
%                 whether a polygon is convex and plays a role in
%                 finding the convex hull of a set of points.
% 
%INPUTS: v1, v2, v3 A set of 3 2XN matrices, of N vertex sets in the order
%                   [x;y] These are 2-dimensional vertices in order,
%                   whereby one wishes to determine whether the angle
%                   v1-v2-v3 is going counterclockwise (left) or clockwise
%                   (right) for each set of vertices.
%
%OUTPUTS: turnDir  An NX1 vector where the ith element is is 1 if the ith
%                  set of vertices form a counterclockwise angle, -1 if
%                  they are going clockwise and 0 if they are collinear or
%                  if two of them coincide.
% 
%The implementation uses the long double datatype to provide a result with
%lower finite precision errors in intermediate results then if everything
%were done using doubles, since most compiler implement the long double as
%an IEEE-754 double extended precision datatype. A notable exception is
%Microsoft's Visual Studio (as of 2017), which maps long doubles to
%doubles (64 bits rather than the 80-bit size of Intel/ AMD floating point
%registers) as noted in 
%https://docs.microsoft.com/en-us/cpp/c-language/type-long-double
%Additionally, the function exactSignOfSumCPP is used to provide the sign
%of an error-free sum. Unfortunately, since the IEEE-754 double extended
%precision datatype only provides a 64-bit mantissa, which is not much of
%an extension over the 54-bit mantissa of a standard double precision
%floating point number, and is not enough to assure no finite precision
%rounding in the products that go into the sum, one cannot guarantee zero
%finite precision errors in this function. On the other hand, if the
%inputs are single-precision (promoted to doubles, since the function does
%not take singles) floating point numbers, the result is exact.
%
%The cross product rule for 3D vectors a and b says that
%norm(cross(a,b))=norm(a)*norm(b)*sin(theta)
%where theta is the positive angle between the vectors. However, If the
%vectors are 2D, (set the z-components to zero), then the cross product
%only have one nonzero component in the z-direction and that component is
%equal to det([a,b]) (for 2D a and b). The interesting thing now, is that
%the sign of the determinant will tell you whether the subsequent vectors
%are going counterclockwise or clockwise. This function returns -1 if the
%vectors are going counterclockwise, 0 if they are exactly collinear and 1
%if they are clockwise.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%turnDir=turnOrientation(v1,v2,v3);
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
