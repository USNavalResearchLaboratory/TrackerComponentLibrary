function yInterp=cubSplineInterpSimp(xPoints,yPoints,xEval)
%%CUBSPINEINTERPSIMP Given a function y(x) evaluated at certain values of
%                 the scalar parameter x (control points), this function
%                 performs cubic spline interpolation at the values in
%                 xEval using some of the closest points in xPoint for
%                 each point in xEval. One obtains interpolated values of
%                 y(xEval), as well as dy(xEval)/dx and d^2y(xEval)/dx^2.
%
%INPUTS: xPoints An NpX1 or 1XNp vector of values at which the function y
%          and its derivative are given.
%        yPoints A 1XNp or an NpX1 matrix of values of the function y(x)
%          evaluated at values of the points in x.
%    xEval An NX1 or 1XN vector of points at which y(x) and its first two
%          derivatives should be interpolated.
%
%OUTPUTS: yInterp A 3XN matrix where yInterp(1,:) is interpolated values of
%                 y(x) and yInterp(2,:) and yInterp(3,:) are respectively
%                 interpolated values of the first and second derivatives
%                 of y(x).
%
%This function just calls fitCubSpline and then evalCubSpline.
%
%Spine interpolation on a collections of points provides smooth first and
%second derivatives. On the other hand, if piecewise Hermite interpolation
%is used, such as when using the pchip function, which is built into
%Matlab, without specifying derivatives, then discontinuities can exist in
%the interpolated derivatives.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    C=fitCubSpline(xPoints,yPoints);
    yInterp=evalCubSpline(xEval,C,xPoints,2);
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
