function [zInterp,gradInterp,gradDerivMatInterp]=interpSurf2D(xyPoints,xyPointsRef,zPointsRef,smoothTerm)
%%INTERPSURF2D  Given a set of values zPointsRef at the 2D points
%               xyPointsRef that define a surface, interpolate to the
%               points given in xyPoints and evaluate the interpolated z-
%               value, the gradient, and the partial derivatives of the
%               interpolated gradient at the points. For large numbers of
%               points, b-splines will be more efficient than this
%               approach.
%
%INPUTS: xyPoints A 2XN set of N 2D (x,y) points at which interplation and
%                 gradients are desired.
%     xyPointsRef The 2XnumRef set of numRef (x,y) points in two-dimensions
%                 that define the locations at which the zPointsRef
%                 (surface heights) are provided. There should be no
%                 repeated (or very numerically close) points and the
%                 points should not all be collinear/ nearly collinear.
%      zPointsRef The numRefX1 or 1XnumRef set of heights at the reference
%                 points xyPointsRef.
%      smoothTerm An optional positive term that afects how smooth the
%                 interpolated surface is. Higher values make the surface
%                 smoother. If omitted, a default value of 0 is used.
%
%OUTPUTS: zInterp The NX1 vector of interpolated z value at the points
%                 given in xyPoints.
%      gradInterp The 2XnumRef gradients of the interpolated surface
%                 evaluated at the points in xyPoints.
% gradDerivMatInterp The 2X2XnumRef second derivative matrices of the
%                 interpolated surface at the points of interest.
%
%The 2D surface from which the derivatives are derived is found using the
%hyperbolic multiquadratic method. The algorithm is cited to originate in
%[1]. However, that paper was not consulted in implementing the routine.
%Rather the summaries given in [2] and [3] where the routine is compared to
%a number of other algorithms were used.
%
%None of those methods explicitly provides the derivatives, but analytic
%derivatives for the interpolating function are simple to find.
%The interpolating function is
%fInterp(x,y)=sum_{i=1}^{numRef} c_i sqrt((x-x_i)^2+(y-y_i)^2+smoothTerm)
%Where the x_i and y_i are the refernce points and the c_i terms are found
%by solving numRef linear equations so that the interpolation equation
%matches the points at points in zPointsRef at the reference x and y
%values.
%
%REFERENCES:
%[1] R. L. Hardy, "Multiquadric equations of topography and other irregular
%    surfaces," Journal of Geophysical Research, vol. 76, no. 8, pp.
%    1905-1915, 10 Mar. 1971.
%[2] S. E. Stead, "Estimation of gradients from scattered data," Journal of
%    Mathematics, vol. 14, no. 1, pp. 265-279, Winter 1984.
%[3] R. E. Barnhill, "A survey of representation and design of surfaces,"
%    IEEE Computer Graphics and Applications, vol. 3, no. 7, pp. 9-16, Oct.
%    1984.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(smoothTerm))
    smoothTerm=0;
end

numRefPoints=size(xyPointsRef,2);
numPoints2Interp=size(xyPoints,2);

%Find all pairwise distances between the N points. The loop does it one row
%at a time. The matrix is symmetric and the number of computations could be
%cut in half, but loops in Matlab are slow, so this is faster.
distMat=zeros(numRefPoints,numRefPoints);
for curRefPointIdx=1:numRefPoints
    diffVec=bsxfun(@minus,xyPointsRef(:,curRefPointIdx),xyPointsRef);
    distMat(curRefPointIdx,:)=sqrt(sum(diffVec.*diffVec,1)+smoothTerm);
end

%Next, solve for the constant coefficients in the interpolation function.
cCoeffs=distMat\zPointsRef;

%Allocate space for the return values
zInterp=zeros(numPoints2Interp,1);
gradInterp=zeros(2,numPoints2Interp);
gradDerivMatInterp=zeros(2,2,numPoints2Interp);

%Perform interpolation using the points.
for curPointIdx=1:numPoints2Interp
    curXYPoint=xyPoints(:,curPointIdx);

    %Find the pairwise distances from the current point and the reference
    %points.
    diffVec=bsxfun(@minus,curXYPoint,xyPointsRef);
    distVals=sqrt(sum(diffVec.*diffVec,1)+smoothTerm);
    
    %The point is linearly interpolated.
    zInterp(curPointIdx)=distVals*cCoeffs;
    
    %Compute the interpolated gradient. This comes from differentiating the
    %formula for zInterp. If smoothTerm=0 and the point is located at one
    %of the reference points, then there will be 0/0 division problems.
    %This checks for the results of such problems and replaces those values
    %with zeros if detected. The limit of those values as smoothTerm-> is
    %zero, so this is a reasonable thing to do.
    temp=diffVec(1,:)./distVals;
    sel=~isfinite(temp);
    temp(sel)=0;
    %Derivative of zInterp with respect to x.
    gradInterp(1,curPointIdx)=temp*cCoeffs;
    temp=diffVec(2,:)./distVals;
    sel=~isfinite(temp);
    temp(sel)=0;
    %Derivative of zInterp with respect to y.
    gradInterp(2,curPointIdx)=temp*cCoeffs;
    
    %Compute the interpolated second derivative matrix. The matrix is
    %symmetric. The same method of avoiding 0/0 division is used here.
    denomVals=distVals.^3;
    temp=diffVec(1,:).*diffVec(2,:)./denomVals;
    sel=~isfinite(temp);
    temp(sel)=0;
    %The two diagonal terms, the d^2/dx*dy terms
    gradDerivMatInterp(1,2,curPointIdx)=-temp*cCoeffs;
    gradDerivMatInterp(2,1,curPointIdx)=gradDerivMatInterp(1,2,curPointIdx);
    
    %The d^2/dx^2 term
    temp=(smoothTerm+diffVec(1,:).^2)./denomVals;
    sel=~isfinite(temp);
    temp(sel)=0;
    gradDerivMatInterp(1,1,curPointIdx)=temp*cCoeffs;
    
    %The d^2/dy^2 term
    temp=(smoothTerm+diffVec(2,:).^2)./denomVals;
    sel=~isfinite(temp);
    temp(sel)=0;
    gradDerivMatInterp(1,1,curPointIdx)=temp*cCoeffs;
end

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
