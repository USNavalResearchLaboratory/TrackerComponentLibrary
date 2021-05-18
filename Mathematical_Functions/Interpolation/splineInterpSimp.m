function yInterp=splineInterpSimp(xPoints,yPoints,xEval)
%%SPINEINTERPSIMP  Given a function y(x) evaluated at certain values of
%                  the scalar parameter x (control points), this function
%                  performs cubic spline interpolation at the values in
%                  xEval using some of the closest points in xPoint for
%                  each point in xEval. One obtains interpolated values of
%                  y(xEval), as well as dy(xEval)/dx and d^2y(xEval)/dx^2.
%
%INPUTS:    x   An NpX1 or 1XNp vector of values at which the function y
%               and its derivative are given.
%           y   A 1XNp or an NpX1 matrix of values of the function y(x)
%               evaluated at values of the points in x.
%       xEval   An NX1 or 1XN vector of points at which y(x) and its
%               first two derivatives should be interpolated.
%
%OUTPUTS: yInterp A 3XN matrix where yInterp(1,:) is interpolated values of
%                 y(x) and yInterp(2,:) and yInterp(3,:) are respectively
%                 interpolated values of the first and second derivatives
%                 of y(x).
%
%This function just calls Matlab's built-in spline function to get the
%piecewise cubic interpolating polynomials. The polynomials are evaluated
%at the points in xEval to get y(xEval). Differentiation routines built
%into Matlab are then applied to the polynomial coefficients to get the
%polynomial forms of the derivative and the second derivative.
%
%Spine interpolation on a collections of points provides smooth first and
%second derivatives. On the other hand, if piecewise Hermite interpolation
%is used, such as when using the pchip function, which is built into
%Matlab, without specifying derivatives, then discontinuities can exist in
%the interpolated derivatives.
%
%April 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    num2Interp=length(xEval);

    yInterp=zeros(3,num2Interp);

    %Get the spline interpolation data. This provides smooth
    %first and second derivatives
    pp=spline(xPoints,yPoints);
    %Interpolate the polynomial value at the given temperature.
    yInterp(1,:)=ppval(pp,xEval);
    %Break apart the interpolation data so that derivatives can
    %be taken.
    [breaks,coefs,l,k,d] = unmkpp(pp);
    %Make the first derivative polynomial.
    pp1 = mkpp(breaks,repmat(k-1:-1:1,[d*l,1]).*coefs(:,1:k-1),d);
    %Evaluate the interpolated first derivatives.
    yInterp(2,:)=ppval(pp1,xEval);
    %Break apart and differentiate the first derivative data to
    %get the second derivatives.
    [breaks,coefs,l,k,d] = unmkpp(pp1);
    pp2 = mkpp(breaks,repmat(k-1:-1:1,[d*l,1]).*coefs(:,1:k-1),d);
    yInterp(3,:)=ppval(pp2,xEval);
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
